module DQ__Time_Integrator__TR
  use Kind_Parameters, only: RNP
  use Constants,       only: ZERO, ONE, THREE, HALF
  use DQ__Time_Integrator

  implicit none
  private

  public :: DQ_TimeIntegrator_TR
  public :: DQ_TimeIntegrator_Options_TR

  !-----------------------------------------------------------------------------
  !> IMEX trapezoidal rule for Dahlquist equation

  type, extends(DQ_TimeIntegrator) :: DQ_TimeIntegrator_TR
  contains
    procedure :: Init_DQ_TimeIntegrator_TR
    procedure :: Show => Show_DQ_TimeIntegrator_TR
    procedure :: TimeStep
  end type DQ_TimeIntegrator_TR

  ! overloading the constructor
  interface DQ_TimeIntegrator_TR
    module procedure New_DQ_TimeIntegrator_TR
  end interface

  !-----------------------------------------------------------------------------
  !> Type for providing TR time-integrator options (none, so far)

  type, extends(DQ_TimeIntegratorOptions) :: DQ_TimeIntegrator_Options_TR
  end type DQ_TimeIntegrator_Options_TR

contains

  !-----------------------------------------------------------------------------
  !> Constructor for objects of type DQ_TimeIntegrator_TR with options

  function New_DQ_TimeIntegrator_TR(opt) result(this)
    class(DQ_TimeIntegrator_Options_TR), optional, intent(in) :: opt
    type(DQ_TimeIntegrator_TR) :: this

    call Init_DQ_TimeIntegrator_TR(this, opt)

  end function New_DQ_TimeIntegrator_TR

  !-----------------------------------------------------------------------------
  !> Initialization of a Init_DQ_TimeIntegrator_TR object

  subroutine Init_DQ_TimeIntegrator_TR(this, opt)
    class(DQ_TimeIntegrator_TR),                   intent(inout) :: this
    class(DQ_TimeIntegrator_Options_TR), optional, intent(in)    :: opt

    ! intialize parent type
    call this % Init_DQ_TimeIntegrator(opt)
    this % name = 'IMEX TR method'

  end subroutine Init_DQ_TimeIntegrator_TR

  !-----------------------------------------------------------------------------
  !> Output of DQ_TimeIntegrator_TR settings

  subroutine Show_DQ_TimeIntegrator_TR(this, unit)
    class(DQ_TimeIntegrator_TR), intent(in) :: this
    integer, optional, intent(in) :: unit  !< output unit

    ! show parent settings
    call this % Show_DQ_TimeIntegrator(unit)

  end subroutine Show_DQ_TimeIntegrator_TR

  !-----------------------------------------------------------------------------
  !> Performs an IMEX TR step

  subroutine TimeStep(this, lambda, dt, z)
    class(DQ_TimeIntegrator_TR), intent(inout) :: this
    real(RNP),    intent(in)    :: dt      !< step size ∆t
    complex(RNP), intent(in)    :: lambda  !< ...
    complex(RNP), intent(inout) :: z       !< z(t) → z(t+ ∆t)

    complex(RNP), parameter :: i = (ZERO, ONE)
    complex(RNP):: z1, z2, z3

    select case(this%impl)

    case(0)

      ! explicit TR ............................................................

      z1 = z
      z2 = z + dt * lambda * z1
      z  = z + dt * HALF * lambda * (z1 + z2)

    case(2)

      ! implicit TR ............................................................

      z = (z + dt * HALF * lambda * z) / (ONE - dt * HALF * lambda)

    case default

      ! IMEX TR ................................................................

      z1 = z

      z2 = z  + dt * i * lambda%im * z1
      z2 = z2 / (ONE - dt * lambda%re)

      z3 = z  + dt * HALF * (lambda * z1 + i * lambda%im * z2 )
      z3 = z3 / (ONE - dt * HALF * lambda%re)

      z = z3

    end select

  end subroutine TimeStep

  !=============================================================================

end module DQ__Time_Integrator__TR
