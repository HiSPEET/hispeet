!> summary:  Interface to incompressible Navier-Stokes problems
!> author:   Joerg Stiller
!> date:     2022/01/06
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module INS__Problem__3D
  use Kind_Parameters, only: RNP
  use Constants,       only: ZERO
  use Execution_Control
  use Array_Assignments
  use MPI_Binding

  implicit none
  private

  public :: INS_Problem_3D

  !-----------------------------------------------------------------------------
  !> Abstract type for defining and handling an incompressible flow problem

  type, abstract :: INS_Problem_3D

    ! model type
    logical :: stokes = .false.         !< switch to Stokes

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
    real(RNP) :: nu_ref  = 1            !< kinematic viscosity

    ! bounding box
    real(RNP) :: x0(3) = 0              !< position nearest to +(∞,∞,∞)
    real(RNP) :: x1(3) = 1              !< position nearest to -(∞,∞,∞)

    ! boundary conditions
    character, allocatable :: bc_v(:)   !< velocity BC type per boundary
    character, allocatable :: bc_p(:)   !< pressure BC type per boundary

  contains

    procedure :: HasVariableProperties
    procedure :: HasExactSolution
    procedure :: GetVariableNames
    procedure :: GetExternalSources
    procedure :: GetViscosity

    procedure :: GetExactSolution
    procedure :: GetExactTimeDerivative
    procedure :: GetExactConvectiveTerm
    procedure :: GetExactPressureTerm
    procedure :: GetExactDiffusiveTerm

    procedure :: SetPressureBC

    procedure(SetProblem),         deferred :: SetProblem
    procedure(GetInitialValues),   deferred :: GetInitialValues
    procedure(GetBoundaryValues),  deferred :: GetBoundaryValues

  end type INS_Problem_3D

  abstract interface

    !---------------------------------------------------------------------------
    !> Initialization of the flow problem
    !>
    !> This routine is intended to set model type and parameters, reference
    !> values, boundary conditions and, possibly, problem-specific parameters.
    !> For flexibility, parameters can be placed in a file. The MPI communicator
    !> is used for broadcasting of parameters and, hence, must be given in case
    !> of parallel execution.

    subroutine SetProblem(problem, nb, file, comm)
      import :: INS_Problem_3D, MPI_Comm
      class(INS_Problem_3D),      intent(inout) :: problem
      integer,                    intent(in)    :: nb    !< number of boundaries
      character(len=*), optional, intent(in)    :: file  !< input file (*.prm)
      type(MPI_Comm),   optional, intent(in)    :: comm  !< MPI communicator
    end subroutine SetProblem

    !---------------------------------------------------------------------------
    !> Provides the initial values u(x,0) for mesh points x

    subroutine GetInitialValues(problem, x, u)
      import :: INS_Problem_3D, RNP
      class(INS_Problem_3D), intent(in) :: problem
      real(RNP), intent(in)  :: x(:,:,:,:,:) !< mesh points
      real(RNP), intent(out) :: u(:,:,:,:,:) !< flow variables
    end subroutine GetInitialValues

    !---------------------------------------------------------------------------
    !> If known, set `ub = u(xb,t)` in boundary points `xb`, or zero otherwise

    subroutine GetBoundaryValues(problem, b, xb, t, ub)
      import :: INS_Problem_3D, RNP
      class(INS_Problem_3D), intent(in) :: problem
      integer,   intent(in)  :: b           !< boundary identifier
      real(RNP), intent(in)  :: xb(:,:,:,:) !< boundary points
      real(RNP), intent(in)  :: t           !< time
      real(RNP), intent(out) :: ub(:,:,:,:) !< flow variables
    end subroutine GetBoundaryValues

  end interface

contains

  !-----------------------------------------------------------------------------
  !> Function for enquiring wether an exact solution is available

  logical function HasExactSolution(problem)
    class(INS_Problem_3D), intent(in) :: problem

    HasExactSolution = problem % exact_solution

  end function HasExactSolution

  !-----------------------------------------------------------------------------
  !> Function for enquiring wether fluid properties are variable

  logical function HasVariableProperties(problem)
    class(INS_Problem_3D), intent(in) :: problem

    HasVariableProperties = problem % variable_props

  end function HasVariableProperties

  !-----------------------------------------------------------------------------
  !> Returns the default variable names

  subroutine GetVariableNames(problem, names)
    class(INS_Problem_3D), intent(in) :: problem
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

  !-----------------------------------------------------------------------------
  !> Automatic generation of pressure boundary conditions

  subroutine SetPressureBC(problem)
    class(INS_Problem_3D), intent(inout) :: problem

    integer :: b

    if (allocated(problem%bc_p)) deallocate(problem%bc_p)
    allocate(problem%bc_p( size(problem%bc_v) ))

    do b = 1, size(problem % bc_v)
      select case(problem % bc_v(b))
      case('D') ! Dirichlet conditions for velocity
        problem % bc_p(b) = 'N'
      case('O') ! Outflow conditions
        problem % bc_p(b) = 'D'
      case('P') ! periodic conditions
        problem % bc_p(b) = 'P'
      case default
        call Error('SetPressureBC', 'invalid velocity BC', 'INS__Problem__3D')
      end select
    end do

  end subroutine SetPressureBC

  !-----------------------------------------------------------------------------
  !> Dummy procedure for external sources

  subroutine GetExternalSources(problem, x, t, F_s)
    class(INS_Problem_3D), intent(in)  :: problem
    real(RNP), intent(in)  :: x(:,:,:,:,:)   !< mesh points
    real(RNP), intent(in)  :: t              !< time
    real(RNP), intent(out) :: F_s(:,:,:,:,:) !< external sources

    call SetArray(F_s, ZERO, multi=.true.)

  end subroutine GetExternalSources

  !-----------------------------------------------------------------------------
  !> Dummy procedure for variable viscosity

  subroutine GetViscosity(problem, x, t, u, nu)
    class(INS_Problem_3D), intent(in) :: problem
    real(RNP), intent(in)  :: x (:,:,:,:,:) !< mesh points
    real(RNP), intent(in)  :: t             !< time
    real(RNP), intent(in)  :: u (:,:,:,:,:) !< flow variables
    real(RNP), intent(out) :: nu(:,:,:,:)   !< viscosity

    call SetArray(nu, problem % nu_ref, multi=.true.)

  end subroutine GetViscosity

  !-----------------------------------------------------------------------------
  !> Dummy procedure for exact solution

  subroutine GetExactSolution(problem, x, t, u)
    class(INS_Problem_3D), intent(in) :: problem
    real(RNP), intent(in)  :: x(:,:,:,:,:) !< mesh points
    real(RNP), intent(in)  :: t            !< time
    real(RNP), intent(out) :: u(:,:,:,:,:) !< solution

    call SetArray(u, ZERO, multi=.true.)

  end subroutine GetExactSolution

  !-----------------------------------------------------------------------------
  !> Dummy procedure for exact time derivative

  subroutine GetExactTimeDerivative(problem, x, t, dt_u)
    class(INS_Problem_3D), intent(in) :: problem
    real(RNP), intent(in)  :: x(:,:,:,:,:)    !< mesh points
    real(RNP), intent(in)  :: t               !< time
    real(RNP), intent(out) :: dt_u(:,:,:,:,:) !< time derivative ∂u/∂t

    call SetArray(dt_u, ZERO, multi=.true.)

  end subroutine GetExactTimeDerivative

  !---------------------------------------------------------------------------
  !> Dummy procedure for exact convection term, F_c = -[∇⋅vv,0]

  subroutine GetExactConvectiveTerm(problem, x, t, F_c)
    class(INS_Problem_3D), intent(in) :: problem
    real(RNP), intent(in)  :: x(:,:,:,:,:)   !< mesh points
    real(RNP), intent(in)  :: t              !< time
    real(RNP), intent(out) :: F_c(:,:,:,:,:) !< convection term

    call SetArray(F_c, ZERO, multi=.true.)

  end subroutine GetExactConvectiveTerm

  !---------------------------------------------------------------------------
  !> Dummy procedure for exact pressure term, F_p = [-∇p,0]

  subroutine GetExactPressureTerm(problem, x, t, F_p)
    class(INS_Problem_3D), intent(in) :: problem
    real(RNP), intent(in)  :: x(:,:,:,:,:)   !< mesh points
    real(RNP), intent(in)  :: t              !< time
    real(RNP), intent(out) :: F_p(:,:,:,:,:) !< pressure term

    call SetArray(F_p, ZERO, multi=.true.)

  end subroutine GetExactPressureTerm

  !---------------------------------------------------------------------------
  !> Dummy procedure for exact diffusion term, F_d = [∇⋅τ,0]

  subroutine GetExactDiffusiveTerm(problem, x, t, F_d)
    class(INS_Problem_3D), intent(in) :: problem
    real(RNP), intent(in)  :: x(:,:,:,:,:)   !< mesh points
    real(RNP), intent(in)  :: t              !< time
    real(RNP), intent(out) :: F_d(:,:,:,:,:) !< diffusion term

    call SetArray(F_d, ZERO, multi=.true.)

  end subroutine GetExactDiffusiveTerm

  !=============================================================================

end module INS__Problem__3D
