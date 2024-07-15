module DQ__Time_Integrator

  use, intrinsic :: ISO_Fortran_Env, only: OUTPUT_UNIT

  use Kind_Parameters, only: RNP
  implicit none
  private

  public :: DQ_TimeIntegrator
  public :: DQ_TimeIntegrator_Options

  !-----------------------------------------------------------------------------
  !> Abstract type of a one-step time integrator for Dahlquist equation

  type, abstract :: DQ_TimeIntegrator

    character(len=80) :: name = ''  !< time-integrator name
    integer :: impl !< switch to explicit/implicit/IMEX corrector (0/1/default)

  contains

    procedure, non_overridable :: Init_DQ_TimeIntegrator
    procedure, non_overridable :: Show_DQ_TimeIntegrator
    procedure :: Show => Show_DQ_TimeIntegrator
    procedure(TimeStep), deferred :: TimeStep
  end type DQ_TimeIntegrator

  abstract interface

    !---------------------------------------------------------------------------
    !> Execution of a single time step

    subroutine TimeStep(this, lambda, dt, u)
      import
      class(DQ_TimeIntegrator), intent(inout) :: this
      real(RNP),    intent(in)    :: dt      !< step size ∆t
      complex(RNP), intent(in)    :: lambda  !< ...
      complex(RNP), intent(inout) :: u       !< u(t) → u(t+ ∆t)
    end subroutine TimeStep

  end interface

  !-----------------------------------------------------------------------------
  !> Base type for providing time integrator options

  type DQ_TimeIntegrator_Options
    integer :: impl = 2  !< 0: explicit, 1: implicit, default:  IMEX
  end type DQ_TimeIntegrator_Options

contains

  !=============================================================================
  ! DQ_TimeIntegrator: type-bound procedures

  !-----------------------------------------------------------------------------
  !> Initialization of DQ_TimeIntegrator object

  subroutine Init_DQ_TimeIntegrator(this, opt)
    class(DQ_TimeIntegrator),                  intent(inout) :: this
    class(DQ_TimeIntegrator_Options), optional, intent(in)    :: opt !< options

    if (present(opt)) then
      this % impl = opt % impl
    end if

  end subroutine Init_DQ_TimeIntegrator

  !-----------------------------------------------------------------------------
  !> Output of DQ_TimeIntegrator settings

  subroutine Show_DQ_TimeIntegrator(this, unit)
    class(DQ_TimeIntegrator), intent(in) :: this
    integer, optional, intent(in) :: unit  !< output unit

    integer :: io

    if (present(unit)) then
      io = unit
    else
      io = OUTPUT_UNIT
    end if

    write(io,'(/,A)') 'DQ_TimeIntegrator settings'
    write(io,'(A,/)') repeat('=',80)

    ! append further settings in corresponding routines of derived types

  end subroutine Show_DQ_TimeIntegrator

  !=============================================================================

end module DQ__Time_Integrator
