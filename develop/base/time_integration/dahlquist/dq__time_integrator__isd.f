module DQ__Time_Integrator__ISD
  use Kind_Parameters, only: RNP
  use Constants,       only: ZERO, ONE, THREE, HALF
  use DQ__Time_Integrator

  implicit none
  private

  public :: DQ_TimeIntegrator_ISD
  public :: DQ_TimeIntegrator_Options_ISD

  !-----------------------------------------------------------------------------
  !> IMEX trapezoidal rule for Dahlquist equation

  type, extends(DQ_TimeIntegrator) :: DQ_TimeIntegrator_ISD
    integer :: order  !< theoretical order of convergence {1,2}
    integer :: method !< 1: MP, 2: TR, default: TR/MP, for order 2 only
  contains
    procedure :: Init_DQ_TimeIntegrator_ISD
    procedure :: Show => Show_DQ_TimeIntegrator_ISD
    procedure :: TimeStep
  end type DQ_TimeIntegrator_ISD

  ! overloading the constructor
  interface DQ_TimeIntegrator_ISD
    module procedure New_DQ_TimeIntegrator_ISD
  end interface

  !-----------------------------------------------------------------------------
  !> Type for providing ISD time-integrator options (none, so far)

  type, extends(DQ_TimeIntegratorOptions) :: DQ_TimeIntegrator_Options_ISD
    integer :: order  = 1 !< theoretical order of convergence {1,2}
    integer :: method = 0 !< 1: MP, 2: TR, default: TR/MP, for order 2 only
  end type DQ_TimeIntegrator_Options_ISD

contains

  !-----------------------------------------------------------------------------
  !> Constructor for objects of type DQ_TimeIntegrator_ISD with options

  function New_DQ_TimeIntegrator_ISD(opt) result(this)
    class(DQ_TimeIntegrator_Options_ISD), optional, intent(in) :: opt
    type(DQ_TimeIntegrator_ISD) :: this

    call Init_DQ_TimeIntegrator_ISD(this, opt)

  end function New_DQ_TimeIntegrator_ISD

  !-----------------------------------------------------------------------------
  !> Initialization of a Init_DQ_TimeIntegrator_ISD object

  subroutine Init_DQ_TimeIntegrator_ISD(this, opt)
    class(DQ_TimeIntegrator_ISD),                   intent(inout) :: this
    class(DQ_TimeIntegrator_Options_ISD), optional, intent(in)    :: opt

    ! intialize parent type
    call this % Init_DQ_TimeIntegrator(opt)

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

  end subroutine Init_DQ_TimeIntegrator_ISD

  !-----------------------------------------------------------------------------
  !> Output of DQ_TimeIntegrator_ISD settings

  subroutine Show_DQ_TimeIntegrator_ISD(this, unit)
    class(DQ_TimeIntegrator_ISD), intent(in) :: this
    integer, optional, intent(in) :: unit  !< output unit

    ! show parent settings
    call this % Show_DQ_TimeIntegrator(unit)

  end subroutine Show_DQ_TimeIntegrator_ISD

  !-----------------------------------------------------------------------------
  !> Performs an IMEX ISD step

  subroutine TimeStep(this, lambda, dt, u)
    class(DQ_TimeIntegrator_ISD), intent(inout) :: this
    real(RNP),    intent(in)    :: dt      !< step size ∆t
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

end module DQ__Time_Integrator__ISD
