!> summary:  Interface to incompressible single-phase flow problems
!> author:   Joerg Stiller
!> date:     2017/08/09
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module ISP_Flow_Problem
  use Kind_Parameters, only: RNP
  use Constants,       only: ZERO
  use Array_Assignments
  use MPI_Binding

  implicit none
  private

  public :: FlowProblem

  !-----------------------------------------------------------------------------
  !> Abstract type for defining and handling an incompressible flow problem

  type, abstract :: FlowProblem

    ! model type
    logical :: stokes         = .false. !< switch to Stokes
    logical :: active_scalar  = .false. !< switch to active scalar(s)
    logical :: passive_scalar = .false. !< switch to passive scalar(s)

    ! model parameters
    integer :: nc         = 4           !< number of components (flow variables)
    integer :: nc_active  = 0           !< number of active scalars
    integer :: nc_passive = 0           !< number of passive scalars

    ! problem characteristics
    logical :: exact_solution = .false. !< switch availability of exact solution
    logical :: variable_props = .false. !< switch to variable properties

    ! reference values
    real(RNP) :: rho_ref = 1            !< density
    real(RNP) :: v_ref   = 1            !< velocity
    real(RNP), allocatable :: nu_ref(:) !< diffusivity

    ! bounding box
    real(RNP) :: x0(3) = 0              !< position nearest to +(∞,∞,∞)
    real(RNP) :: x1(3) = 1              !< position nearest to -(∞,∞,∞)

    ! boundary conditions
    character, allocatable :: bc(:,:)   !< BC type per boundary and variable

  contains

    procedure :: HasVariableProperties
    procedure :: HasExactSolution
    procedure :: GetExactSolution
    procedure :: GetExactTimeDerivative
    procedure :: GetVariableNames

    procedure(SetProblem),                deferred :: SetProblem
    procedure(GetInitialValues),          deferred :: GetInitialValues
    procedure(GetBoundaryValues),         deferred :: GetBoundaryValues
    procedure(GetBoundaryTimeDerivative), deferred :: GetBoundaryTimeDerivative
    procedure(GetExternalSources),        deferred :: GetExternalSources

  end type FlowProblem

  abstract interface

    !---------------------------------------------------------------------------
    !> Initialization of the flow problem
    !>
    !> This routine is intended to set model type and parameters, reference
    !> values, boundary conditions and, possibly, problem-specific parameters.
    !> For flexibility, parameters can be placed in a file. The MPI communicator
    !> is used for broadcasting of parameters and, hence, must be given in case
    !> of parallel execution.

    subroutine SetProblem(problem, file, comm)
      import
      class(FlowProblem),         intent(inout) :: problem
      character(len=*), optional, intent(in)    :: file    !< input file (*.prm)
      type(MPI_Comm),   optional, intent(in)    :: comm    !< MPI communicator
    end subroutine SetProblem

    !---------------------------------------------------------------------------
    !> Provides the initial values u(x,0) for mesh points x

    subroutine GetInitialValues(problem, x, u)
      import
      class(FlowProblem), intent(in)  :: problem
      real(RNP),          intent(in)  :: x(:,:,:,:,:) !< mesh points
      real(RNP),          intent(out) :: u(:,:,:,:,:) !< flow variables
    end subroutine GetInitialValues

    !---------------------------------------------------------------------------
    !> If known, set `ub = u(xb,t)` on boundary `b`, or zero otherwise

    subroutine GetBoundaryValues(problem, b, xb, t, ub)
      import
      class(FlowProblem), intent(in)  :: problem
      integer,            intent(in)  :: b              !< boundary ID
      real(RNP),          intent(in)  :: xb(:,:,:,:)    !< boundary points
      real(RNP),          intent(in)  :: t              !< time
      real(RNP),          intent(out) :: ub(:,:,:,:)    !< flow variables
    end subroutine GetBoundaryValues

    !---------------------------------------------------------------------------
    !> If known, set `dt_ub = ∂u/∂t(xb,t)` on boundary `b`, or zero otherwise

    subroutine GetBoundaryTimeDerivative(problem, b, xb, t, dt_ub)
      import
      class(FlowProblem), intent(in)  :: problem
      integer,            intent(in)  :: b              !< boundary ID
      real(RNP),          intent(in)  :: xb(:,:,:,:)    !< boundary points
      real(RNP),          intent(in)  :: t              !< time
      real(RNP),          intent(out) :: dt_ub(:,:,:,:) !< ∂u/∂t
    end subroutine GetBoundaryTimeDerivative

    !---------------------------------------------------------------------------
    !> Provides the external sources for all variables at points x and time t

    subroutine GetExternalSources(problem, x, t, f)
      import
      class(FlowProblem), intent(in)  :: problem
      real(RNP),          intent(in)  :: x(:,:,:,:,:) !< mesh points
      real(RNP),          intent(in)  :: t            !< time
      real(RNP),          intent(out) :: f(:,:,:,:,:) !< external sources
    end subroutine GetExternalSources

  end interface

contains

!-------------------------------------------------------------------------------
!> Function for enquiring wether an exact solution is available

logical function HasExactSolution(problem)
  class(FlowProblem), intent(in) :: problem

  HasExactSolution = problem % exact_solution

end function HasExactSolution

!-------------------------------------------------------------------------------
!> Function for enquiring wether fluid properties are variable

logical function HasVariableProperties(problem)
  class(FlowProblem), intent(in) :: problem

  HasVariableProperties = problem % variable_props

end function HasVariableProperties

!-------------------------------------------------------------------------------
!> Dummy procedure for exact solution

subroutine GetExactSolution(problem, x, t, u)
  class(FlowProblem), intent(in)  :: problem
  real(RNP),          intent(in)  :: x(:,:,:,:,:) !< mesh points
  real(RNP),          intent(in)  :: t            !< time
  real(RNP),          intent(out) :: u(:,:,:,:,:) !< solution

  call SetArray(u, ZERO, multi=.true.)

end subroutine GetExactSolution

!-------------------------------------------------------------------------------
!> Dummy procedure for exact time derivative

subroutine GetExactTimeDerivative(problem, x, t, dt_u)
  class(FlowProblem), intent(in)  :: problem
  real(RNP),          intent(in)  :: x(:,:,:,:,:)    !< mesh points
  real(RNP),          intent(in)  :: t               !< time
  real(RNP),          intent(out) :: dt_u(:,:,:,:,:) !< time derivative

  call SetArray(dt_u, ZERO, multi=.true.)

end subroutine GetExactTimeDerivative

!-------------------------------------------------------------------------------
!> Returns the default variable names

subroutine GetVariableNames(problem, names)
  class(FlowProblem),            intent(in)  :: problem
  character(len=:), allocatable, intent(out) :: names(:)  !< variable names

  integer :: i, k, l

  l = 5 + ceiling(log10(real( max(1, problem%nc_active, problem%nc_passive) )))

  allocate(character(len=l) :: names(problem % nc))

  names(1) = 'v_1'
  names(2) = 'v_2'
  names(3) = 'v_3'
  names(4) = 'p'

  k = 4

  do i = 1, problem % nc_active
    k = k + 1
    write(names(k), '(A,I0)') 'Phi_', i
  end do

  do i = 1, problem % nc_passive
    k = k + 1
    write(names(k), '(A,I0)') 'Psi_', i
  end do

end subroutine GetVariableNames

!===============================================================================

end module ISP_Flow_Problem
