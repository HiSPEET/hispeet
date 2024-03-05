module DQ__Time_Integrator__ISD

  use, intrinsic :: ISO_Fortran_Env, only: OUTPUT_UNIT

  use Kind_Parameters, only: RNP
  use Constants,       only: ZERO, ONE, THREE, HALF
  use DQ__Time_Integrator

  implicit none
  private

  public :: DQ_TimeIntegrator_ISD
  public :: DQ_TimeIntegrator_Options_ISD

  !-----------------------------------------------------------------------------
  !> ISD methods for Dahlquist equation

  type, extends(DQ_TimeIntegrator) :: DQ_TimeIntegrator_ISD
    integer   :: order   !< theoretical order of convergence {1,2}
    integer   :: n_stage !< number of stages {1,2}, ignored with order 2
    real(RNP) :: sigma   !< ISD scaling factor
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
    integer   :: order   = 1    !< theoretical order of convergence
    integer   :: n_stage = 1    !< number of stages, ignored with order 2
    real(RNP) :: sigma   = HALF !< ISD scaling factor
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

    this % order   = opt % order
    this % n_stage = opt % n_stage
    this % sigma   = opt % sigma

    select case(this % order)
    case(1)
      select case(this % n_stage)
      case(2)
        this % name = 'ISD method of order 1 with two stages'
      case default
        this % name = 'ISD method of order 1 with one stage'
      end select
    case(2)
      this % name = 'ISD method of order 2 with two stages'
    end select

  end subroutine Init_DQ_TimeIntegrator_ISD

  !-----------------------------------------------------------------------------
  !> Output of DQ_TimeIntegrator_ISD settings

  subroutine Show_DQ_TimeIntegrator_ISD(this, unit)
    class(DQ_TimeIntegrator_ISD), intent(in) :: this
    integer, optional, intent(in) :: unit  !< output unit

    integer :: io

    if (present(unit)) then
      io = unit
    else
      io = OUTPUT_UNIT
    end if

    ! show parent settings
    call this % Show_DQ_TimeIntegrator(unit)

    write(io,'(2X,A,T15,G0)')     'name:' , trim(this % name)
    write(io,'(2X,A,T14,ES12.5)') 'sigma:', this % sigma

  end subroutine Show_DQ_TimeIntegrator_ISD

  !-----------------------------------------------------------------------------
  !> Performs an IMEX ISD step

  subroutine TimeStep(this, lambda, dt, u)
    class(DQ_TimeIntegrator_ISD), intent(inout) :: this
    real(RNP),    intent(in)    :: dt      !< step size ∆t
    complex(RNP), intent(in)    :: lambda  !< ...
    complex(RNP), intent(inout) :: u       !< u(t) → u(t+ ∆t)

    complex(RNP), parameter :: i = (ZERO, ONE)
    complex(RNP) :: u0, u1, u2, u3
    real(RNP) :: hdt

    associate(sigma => this % sigma)

      hdt = HALF * dt

      select case(this%order)

      case(1)

        select case(this%n_stage)

        case(2)

          ! two-stage ............................................................

          u0 = u

          u1 = u0 + dt * (ZERO, ONE) * lambda%im * u0
          u1 = u1 / (ONE - dt * (lambda%re - sigma * dt * lambda%im**2))

          u2 = u0 + dt * (ZERO, ONE) * lambda%im * u1
          u2 = u2 / (ONE - dt * (lambda%re - sigma * dt * lambda%im**2))

          u = u2

        case default

          ! one-stage ............................................................

          u = u + dt * (ZERO, ONE) * lambda%im * u
          u = u / (ONE - dt * (lambda%re - sigma * dt * lambda%im**2))

        end select

      case(2)

        u0 = u

        u1 = u0 +  hdt * i * lambda%im * u0
        u1 = u1 / (ONE - sigma * dt * lambda%re + (hdt * lambda%im)**2)

        u2 = u0 +  hdt * i * lambda%im * u1
        u2 = u2 / (ONE - sigma * dt * lambda%re + (hdt * lambda%im)**2)

        u = u0 +  dt * lambda * u2

      end select

    end associate

  end subroutine TimeStep

  !=============================================================================

end module DQ__Time_Integrator__ISD
