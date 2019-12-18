!> summary:  Interface to 1D conservation problems discretized with DG-SEM
!> author:   Joerg Stiller
!> date:     2019/11/08
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Conservation_Problem_1D

  use Kind_Parameters, only: RNP
  use IP_Element_Operators_1D
  use DG_Utilities_1D

  implicit none
  private

  public :: ConservationProblem1D

  !-----------------------------------------------------------------------------
  !> Abstract type for defining and handling 1D conservation problems

  type, abstract :: ConservationProblem1D

    ! problem parameters .......................................................

    integer                     :: nc = -1 !< number of conservation variables
    real(RNP)                   :: xb1     !< position of 1st boundary (left)
    real(RNP)                   :: xb2     !< position of 2nd boundary (right)
    character, allocatable      :: bc(:,:) !< BC type per boundary and variable

    ! space discretization .....................................................

    integer                     :: ne = -1 !< number of elements
    real(RNP)                   :: dx      !< element size
    type(IP_ElementOperators1D) :: eop     !< IP-DG element operators
    real(RNP), allocatable      :: x(:,:)  !< mesh points

    ! private control parameters ...............................................

    logical, private :: exact_solution     = .false.
    logical, private :: implicit_diffusion = .false.

  contains

    procedure :: SetSpaceDiscretization
    procedure :: HasExactSolution
    procedure :: ExactSolution

    procedure(SetProblem),     deferred :: SetProblem
    procedure(InitialValues),  deferred :: InitialValues
    procedure(ConvectiveFlux), deferred :: ConvectiveFlux
    procedure(DiffusiveFlux),  deferred :: DiffusiveFlux
    procedure(RHS_Convection), deferred :: RHS_Convection
    procedure(RHS_Diffusion),  deferred :: RHS_Diffusion

  end type ConservationProblem1D

  abstract interface

    !---------------------------------------------------------------------------
    !> Initialization of the flow problem

    subroutine SetProblem(problem, file)
      import
      class(ConservationProblem1D), intent(inout) :: problem
      character(len=*), optional,   intent(in)    :: file !< (*.prm)
    end subroutine SetProblem

    !---------------------------------------------------------------------------
    !> Provides the initial values u(x,0) for mesh points x

    function InitialValues(problem) result(u)
      import
      class(ConservationProblem1D), intent(in) :: problem
      real(RNP) :: u(0:problem%eop%po, problem%ne, problem%nc)
    end function InitialValues

    !---------------------------------------------------------------------------
    !> Convective flux of the given conservation problem

    function ConvectiveFlux(problem, u) result(fc)
      import
      class(ConservationProblem1D), intent(in) :: problem
      real(RNP), contiguous,        intent(in) :: u(0:,:,:) !< u(x,t)

      real(RNP) :: fc(0:ubound(u,1), size(u,2), problem%nc)

    end function ConvectiveFlux

    !---------------------------------------------------------------------------
    !> Diffusive flux of the given conservation problem

    function DiffusiveFlux(problem, u) result(fd)
      import
      class(ConservationProblem1D), intent(in) :: problem
      real(RNP), contiguous,        intent(in) :: u   (0:,:,:) !< u(x,t)

      real(RNP) :: fd(0:ubound(u,1), size(u,2), problem%nc)

    end function DiffusiveFlux

    !---------------------------------------------------------------------------
    !> Convective contribution to RHS of DG-SEM formulation

    function RHS_Convection(problem, t, u) result(rc)
      import
      class(ConservationProblem1D), intent(in) :: problem
      real(RNP),                    intent(in) :: t         !< time
      real(RNP), contiguous,        intent(in) :: u(0:,:,:) !< u(x,t)

      real(RNP) :: rc(0:ubound(u,1), size(u,2), problem%nc)

    end function RHS_Convection

    !---------------------------------------------------------------------------
    !> Diffusive contribution to RHS of DG-SEM formulation

    function RHS_Diffusion(problem, t, u) result(rd)
      import
      class(ConservationProblem1D), intent(in) :: problem
      real(RNP),                    intent(in) :: t         !< time
      real(RNP), contiguous,        intent(in) :: u(0:,:,:) !< u(x,t)

      real(RNP) :: rd(0:ubound(u,1), size(u,2), problem%nc)

    end function RHS_Diffusion

  end interface

contains

  !-----------------------------------------------------------------------------
  !> Initialization of space-discretization

  subroutine SetSpaceDiscretization(problem, opt, ne)
    class(ConservationProblem1D), intent(inout) :: problem
    class(IP_ElementOptions1D),   intent(in)    :: opt !< IP/DG-SEM options
    integer,                      intent(in)    :: ne  !< number of elements

    problem % eop = IP_ElementOperators1D(opt)
    problem % ne  = ne

    if (allocated(problem % x)) then
      deallocate(problem % x)
    end if
    allocate(problem % x(0:problem%eop%po, ne))

    call GetMeshPoints( problem % eop  &
                      , problem % xb1  &
                      , problem % xb2  &
                      , problem % dx   &
                      , problem % x    )

  end subroutine SetSpaceDiscretization

  !-----------------------------------------------------------------------------
  !> Function for enquiring wether an exact solution is available

  logical function HasExactSolution(problem)
    class(ConservationProblem1D), intent(in) :: problem

    HasExactSolution = problem % exact_solution

  end function HasExactSolution

  !-----------------------------------------------------------------------------
  !> Dummy procedure for the exact solution u(x,t)

  function ExactSolution(problem, x, t) result(u)
    class(ConservationProblem1D), intent(in) :: problem
    real(RNP), contiguous,        intent(in) :: x(0:,:)
    real(RNP),                    intent(in) :: t
    real(RNP) :: u(0:ubound(x,1), size(x,2), problem%nc)

    u = 0

    ! just to avoid compiler warnings ;)
    if (abs(t) >= 0) return

  end function ExactSolution

  !=============================================================================

end module Conservation_Problem_1D
