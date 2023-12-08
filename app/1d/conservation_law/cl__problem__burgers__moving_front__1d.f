!> summary:  Moving front problem for Burgers equation
!> author:   Joerg Stiller
!> date:     2023/05/04
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> This problem is based on an exact solution of Burgers equation in the domain
!> [-1,1]:
!>
!>       u(x,t) = 1 - tanh[(x + 0.5 - t) / (2ν)],
!>
!> given in
!> J.S. Hesthaven & T. Warburton, Nodal Discontinuous Galerkin Methods,
!> Springer 2008, on page 255.
!>
!> In the inviscid limit, the solution takes the form
!>
!>       u(x,t) = 2   for x < xs = t - 0.5 and
!>       u(x,t) = 0   for x > xs.
!>
!===============================================================================

module CL__Problem__Burgers__Moving_Front__1D

  use Kind_Parameters, only: RNP
  use Constants,       only: ONE, HALF, ZERO
  use Execution_Control
  use CL__Problem__Burgers__1D
  use CL__Operator__1D

  implicit none
  private

  public :: CL_Problem_Burgers_MovingFront_1D

  !-----------------------------------------------------------------------------
  !> Type for defining and handling 1D Burgers moving front problem

  type, extends(CL_Problem_Burgers_1D) :: CL_Problem_Burgers_MovingFront_1D
  contains

    procedure :: SetProblem
    procedure :: GetInitialValues
    procedure :: GetExactSolution
    procedure :: GetBoundaryValues

  end type CL_Problem_Burgers_MovingFront_1D

  !=============================================================================

contains

  !-----------------------------------------------------------------------------
  !> Initialization of the Burgers moving front problem

  subroutine SetProblem(this, file)
    class(CL_Problem_Burgers_MovingFront_1D), intent(inout) :: this
    character(len=*), optional, intent(in) :: file !< (*.prm)

    real(RNP) :: nu = 0.1
    character :: bc(2) = ['D','N']
    integer   :: nu_sd_filter = -1

    namelist /burgers_moving_front_prm/ nu, bc, nu_sd_filter

    logical :: exists, opened
    integer :: prm

    ! check for input file .....................................................

    if (present(file)) then

      inquire(file=trim(file)//'.prm', exist=exists, opened=opened, number=prm)

      if (exists .and. .not. opened) then
        open(newunit=prm, file=trim(file)//'.prm', action='READ')
      else if (opened) then
        rewind(prm)
      else
        call Error( 'SetProblem'                                  &
                  , 'Input file "'//trim(file)//'.prm" not found' &
                  , 'CL__Problem__Burgers__Moving_Front__1D'      )
      end if

      ! read parameters ........................................................

      read(prm, nml=burgers_moving_front_prm)
      if (.not. opened) close(prm)

    end if

    ! set parameters ...........................................................

    this % nc  =  1  ! number of conservation variables
    this % xb1 = -1  ! position of left boundary
    this % xb2 =  1  ! position of right boundary
    this % bc  =  bc ! BC types at left and right boundaries
    this % nu  =  nu ! viscosity

    this % nu_sd_filter = nu_sd_filter
    this % has_exact_solution = .true.

  end subroutine SetProblem

  !---------------------------------------------------------------------------
  !> Provides the initial values u(x,0)

  subroutine GetInitialValues(this, cl_operator, u)
    class(CL_Problem_Burgers_MovingFront_1D), intent(in) :: this
    class(CL_Operator_1D), intent(in)  :: cl_operator
    real(RNP), contiguous, intent(out) :: u(0:,:,:) !< u⁰(0:po,1:n1,1:nc)

    u(:,:,1) = ExactSolution(this%nu, cl_operator%x, ZERO)

  end subroutine GetInitialValues

  !-----------------------------------------------------------------------------
  !> Exact solution u(x,t)

  subroutine GetExactSolution(this, cl_operator, t, u)
    class(CL_Problem_Burgers_MovingFront_1D), intent(in) :: this
    class(CL_Operator_1D), intent(in)  :: cl_operator
    real(RNP),             intent(in)  :: t
    real(RNP), contiguous, intent(out) :: u(0:,:,:) !< u(x,t) at mesh points

    u(:,:,1) = ExactSolution(this%nu, cl_operator%x, t)

  end subroutine GetExactSolution

  !---------------------------------------------------------------------------
  !> Provides the left and right boundary values for time t

  subroutine GetBoundaryValues(this, t, bv)
    class(CL_Problem_Burgers_MovingFront_1D), intent(in) :: this
    real(RNP), intent(in)  :: t       !< time
    real(RNP), intent(out) :: bv(:,:) !< boundary values

    ! left boundary
    select case(this % bc(1))
    case('D')
      bv(1,1) = ExactSolution(this%nu, this%xb1, t)
    case('N')
      bv(1,1) = ExactDerivative(this%nu, this%xb1, t)
    case default
      bv(1,1) = 0
    end select

    ! right boundary
    select case(this % bc(2))
    case('D')
      bv(1,2) = ExactSolution(this%nu, this%xb2, t)
    case('N')
      bv(1,2) = ExactDerivative(this%nu, this%xb2, t)
    case default
      bv(1,2) = 0
    end select

  end subroutine GetBoundaryValues

  !-----------------------------------------------------------------------------
  !> Exact solution u(x,t)

  elemental function ExactSolution(nu, x, t) result(u)
    real(RNP), intent(in) :: nu
    real(RNP), intent(in) :: x
    real(RNP), intent(in) :: t
    real(RNP) :: u

    if (nu > 0) then
      u = ONE - tanh((x + HALF - t) / (2*nu))
    else
      if (x < t - HALF) then
        u = 2
      else
        u = 0
      end if
    end if

  end function ExactSolution

  !-----------------------------------------------------------------------------
  !> Exact viscous flux q = ν ∂u/∂x(x,t)

  elemental function ExactDerivative(nu, x, t) result(q)
    real(RNP), intent(in) :: nu
    real(RNP), intent(in) :: x
    real(RNP), intent(in) :: t
    real(RNP) :: q

    if (nu > 0) then
      q = ONE/(2*nu) * (tanh((x + HALF - t) / (2*nu))**2 - ONE)
    else
      q = 0
    end if

  end function ExactDerivative

  !=============================================================================

end module CL__Problem__Burgers__Moving_Front__1D

