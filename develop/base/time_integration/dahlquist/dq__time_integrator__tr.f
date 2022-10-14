module DQ__Time_Integrator__TR

  use, intrinsic :: ISO_Fortran_Env, only: OUTPUT_UNIT

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

    integer :: io

    if (present(unit)) then
      io = unit
    else
      io = OUTPUT_UNIT
    end if

    ! show parent settings
    call this % Show_DQ_TimeIntegrator(unit)

    write(io,'(2X,A,T15,G0)')  'name:', trim(this % name)
    write(io,'(2X,A,T15,G0)')  'impl:', this % impl

  end subroutine Show_DQ_TimeIntegrator_TR

  !-----------------------------------------------------------------------------
  !> Performs an IMEX TR step

  subroutine TimeStep(this, lambda, dt, u)
    class(DQ_TimeIntegrator_TR), intent(inout) :: this
    real(RNP),    intent(in)    :: dt      !< step size ∆t
    complex(RNP), intent(in)    :: lambda  !< ...
    complex(RNP), intent(inout) :: u       !< u(t) → u(t+ ∆t)

    complex(RNP), parameter :: i = (ZERO, ONE)
    complex(RNP):: u1, u2, u3

    select case(this%impl)

    case(0)

      ! explicit TR ............................................................

      u1 = u
      u2 = u + dt * lambda * u1
      u  = u + dt * HALF * lambda * (u1 + u2)

    case(1)

      ! implicit TR ............................................................

      u = (u + dt * HALF * lambda * u) / (ONE - dt * HALF * lambda)

    case default

      ! IMEX TR ................................................................

      u1 = u

      u2 = u  + dt * i * lambda%im * u1
      u2 = u2 / (ONE - dt * lambda%re)

      u3 = u  + dt * HALF * (lambda * u1 + i * lambda%im * u2 )
      u3 = u3 / (ONE - dt * HALF * lambda%re)

      u = u3

    end select

  end subroutine TimeStep

  !=============================================================================

end module DQ__Time_Integrator__TR
