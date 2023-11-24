module CL__Time_Integrator__ISD__1D

  use, intrinsic :: ISO_Fortran_Env, only: OUTPUT_UNIT

  use Kind_Parameters, only: RNP
  use Constants,       only: ZERO, ONE, THREE, HALF

  use CL__Problem__Scalar__1D
  use CL__Time_Integrator__1D

  implicit none
  private

  public :: CL_TimeIntegrator_ISD_1D
  public :: CL_TimeIntegrator_Options_ISD_1D

  !-----------------------------------------------------------------------------
  !> ...

  type, extends(CL_TimeIntegrator_1D) :: CL_TimeIntegrator_ISD_1D
    integer :: order  !< theoretical order of convergence {1,2}
    integer :: method !< 1: MP, 2: TR, default: TR/MP, for order 2 only
  contains
    procedure :: Init_CL_TimeIntegrator_Options_ISD_1D
    procedure :: Show => Show_CL_TimeIntegrator_Options_ISD_1D
    procedure :: TimeStep
  end type CL_TimeIntegrator_ISD_1D

  ! overloading the constructor
  interface CL_TimeIntegrator_ISD_1D
    module procedure New_CL_TimeIntegrator_Options_ISD_1D
  end interface

  !-----------------------------------------------------------------------------
  !> Type for providing ISD time-integrator options (none, so far)

  type, extends(CL_TimeIntegrator_Options_1D) :: &
      CL_TimeIntegrator_Options_ISD_1D
    integer :: order  = 1 !< theoretical order of convergence {1,2}
    integer :: method = 0 !< 1: MP, 2: TR, default: TR/MP, for order 2 only
  end type CL_TimeIntegrator_Options_ISD_1D

contains

  !-----------------------------------------------------------------------------
  !> Constructor for objects of type CL_TimeIntegrator_ISD_1D with options

  function New_CL_TimeIntegrator_Options_ISD_1D(opt) result(this)
    class(CL_TimeIntegrator_Options_ISD_1D), optional, intent(in) :: opt
    type(CL_TimeIntegrator_ISD_1D) :: this

    call Init_CL_TimeIntegrator_Options_ISD_1D(this, opt)

  end function New_CL_TimeIntegrator_Options_ISD_1D

  !-----------------------------------------------------------------------------
  !> Initialization of a Init_CL_TimeIntegrator_Options_ISD_1D object

  subroutine Init_CL_TimeIntegrator_Options_ISD_1D(this, opt)
    class(CL_TimeIntegrator_ISD_1D),                   intent(inout) :: this
    class(CL_TimeIntegrator_Options_ISD_1D), optional, intent(in)    :: opt

    ! intialize parent type
    call this % Init_CL_TimeIntegrator_1D(opt)

    this % order  = opt % order
    this % method = opt % method

    select case(this % order)
    case(1)
      this % name = 'ISD IMEX method of order 1'
    case(2)
      select case(this % method)
      case(1)
        this % name = 'ISD IMEX midpoint rule'
      case(2)
        this % name = 'ISD IMEX trapezoidal rule'
      case default
        this % name = 'ISD IMEX trapezoidal/midpoint rule'
      end select
    end select

  end subroutine Init_CL_TimeIntegrator_Options_ISD_1D

  !-----------------------------------------------------------------------------
  !> Output of CL_TimeIntegrator_ISD_1D settings

  subroutine Show_CL_TimeIntegrator_Options_ISD_1D(this, unit)
    class(CL_TimeIntegrator_ISD_1D), intent(in) :: this
    integer, optional, intent(in) :: unit  !< output unit

    integer :: io

    if (present(unit)) then
      io = unit
    else
      io = OUTPUT_UNIT
    end if

    ! show parent settings
    call this % Show_CL_TimeIntegrator_1D(unit)

    write(io,'(2X,A,T15,G0)')  'name:' , trim(this % name)
    write(io,'(2X,A,T15,G0)')  'order:', this % impl

  end subroutine Show_CL_TimeIntegrator_Options_ISD_1D

  !-----------------------------------------------------------------------------
  !> Performs an IMEX ISD step

  subroutine TimeStep(this, problem, t, dt, u)
    class(CL_TimeIntegrator_ISD_1D), intent(inout) :: this
    class(CL_Problem_Scalar_1D),     intent(in)    :: problem
    real(RNP),    intent(inout) :: t
    real(RNP),    intent(in)    :: dt       !< step size ∆t
    complex(RNP), intent(in)    :: lambda  !< ...
    complex(RNP), intent(inout) :: u       !< u(t) → u(t+ ∆t)

    complex(RNP), parameter :: i = (ZERO, ONE)
    complex(RNP) :: u1, u2, u3
    real(RNP) :: hdt

    select case(this%order)

    case(1)

      u = u + dt * (ZERO, ONE) * lambda%im * u
      u = u / (ONE - dt * lambda%re + (dt * lambda%im)**2)

    case(2)

      hdt = HALF * dt

      select case(this%method)

      case(1)

        ! midpoint .............................................................

        u1 = u

        u2 = u  +  hdt * i * lambda%im * u1
        u2 = u2 / (ONE  - hdt * lambda%re + (hdt * lambda%im)**2)

        u3 = u  +  dt * lambda * u2

        u = u3

      case(2)

        ! trapezoid ............................................................

        u1 = u

        u2 = u  +  dt * i * lambda%im * u1
        u2 = u2 / (ONE - dt * lambda%re + (dt * lambda%im)**2)

        u3 = u  +  hdt * (lambda * u1 + i * lambda%im * u2 )
        u3 = u3 / (ONE - hdt * lambda%re)

        u = u3

      case default

        ! trapezoid/midpoint ...................................................

        u1 = u

        u2 = u  +  hdt * i * lambda%im * u1
        u2 = u2 / (ONE  - hdt * lambda%re + (hdt * lambda%im)**2)

        u3 = u  +  hdt * lambda%re * u1  +  dt * i * lambda%im * u2
        u3 = u3 / (ONE - hdt * lambda%re)

        u = u3

      end select

    end select

  end subroutine TimeStep

  !=============================================================================

end module CL__Time_Integrator__ISD__1D
