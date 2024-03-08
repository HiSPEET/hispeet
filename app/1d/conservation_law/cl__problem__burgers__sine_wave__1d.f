!> summary:  Sine wave problem for Burgers equation
!> author:   Joerg Stiller
!> date:     2024/01/29
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> This problem considers the solution of Burgers equation in the domain [0,1]
!> with initial conditions
!>
!>       u₀(x) = c₀ + sin(2πx)
!>
!===============================================================================

module CL__Problem__Burgers__Sine_Wave__1D

  use Kind_Parameters, only: RNP
  use Constants,       only: HALF, PI
  use Execution_Control
  use CL__Problem__Burgers__1D
  use CL__Operator__1D

  implicit none
  private

  public :: CL_Problem_Burgers_SineWave_1D
  public :: CL_Problem_Burgers_SineWave_Options_1D

  !-----------------------------------------------------------------------------
  !> Type for defining and handling 1D Burgers moving front problem

  type, extends(CL_Problem_Burgers_1D) :: CL_Problem_Burgers_SineWave_1D

    real(RNP) :: c0

  contains

    procedure :: SetProblem
    procedure :: GetInitialValues
    procedure :: GetBoundaryValues

  end type CL_Problem_Burgers_SineWave_1D

  !-----------------------------------------------------------------------------
  !> 1D Burgers moving front options

  type, extends(CL_Problem_Burgers_Options_1D) :: &
      CL_Problem_Burgers_SineWave_Options_1D

    real(RNP) :: c0 = HALF

  end type CL_Problem_Burgers_SineWave_Options_1D

  !=============================================================================

contains

  !-----------------------------------------------------------------------------
  !> Initialization of the Burgers moving front problem

  subroutine SetProblem(this, file)
    class(CL_Problem_Burgers_SineWave_1D), intent(inout) :: this
    character(len=*), optional, intent(in) :: file !< (*.prm)

    type(CL_Problem_Burgers_SineWave_Options_1D) :: opt
    namelist /burgers_moving_front_prm/ opt

    logical :: exists, opened
    integer :: prm

    ! preset options ...........................................................

    opt % xb1 =  0
    opt % xb2 =  1
    opt % nu  =  1e-3
    opt % bc  = 'P'

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
                  , 'CL__Problem__Burgers__Sine_Wave__1D'         )
      end if

      ! read parameters ........................................................

      read(prm, nml=burgers_moving_front_prm)
      if (.not. opened) close(prm)

    end if

    ! set parameters ...........................................................

    ! initialize base type
    call this % Init_CL_Problem_Burgers_1D(opt)

    this % c0 = opt % c0

  end subroutine SetProblem

  !---------------------------------------------------------------------------
  !> Provides the initial values u(x,0)

  subroutine GetInitialValues(this, cl_operator, u)
    class(CL_Problem_Burgers_SineWave_1D), intent(in) :: this
    class(CL_Operator_1D), intent(in)  :: cl_operator
    real(RNP), contiguous, intent(out) :: u(0:,:,:) !< u⁰(0:po,1:n1,1:nc)

    u(:,:,1) = this%c0 + sin(2 * PI * cl_operator%x)

  end subroutine GetInitialValues

  !---------------------------------------------------------------------------
  !> Provides the left and right boundary values for time t

  subroutine GetBoundaryValues(this, t, bv)
    class(CL_Problem_Burgers_SineWave_1D), intent(in) :: this
    real(RNP), intent(in)  :: t       !< time
    real(RNP), intent(out) :: bv(:,:) !< boundary values

    if (any(this%bc /= 'P')) then
      call Error( 'GetBoundaryValues'                       &
                , 'Problem supports only periodic boundary' &
                , 'CL__Problem__Burgers__Sine_Wave__1D'     )
    end if

  end subroutine GetBoundaryValues

  !=============================================================================

end module CL__Problem__Burgers__Sine_Wave__1D

