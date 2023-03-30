module CL__Time_Integrator__1D

  use, intrinsic :: ISO_Fortran_Env, only: OUTPUT_UNIT

  use Kind_Parameters, only: RNP
  use CL__Problem__Scalar__1D

  implicit none
  private

  public :: CL_TimeIntegrator_1D
  public :: CL_TimeIntegrator_Options_1D

  !-----------------------------------------------------------------------------
  !> Abstract type of a one-step time integrator for Dahlquist equation

  type, abstract :: CL_TimeIntegrator_1D

    character(len=80) :: name = ''  !< time-integrator name
    integer :: impl !< switch to explicit/implicit/IMEX corrector (0/1/2)

  contains

    procedure, non_overridable :: Init_CL_TimeIntegrator_1D
    procedure, non_overridable :: Show_CL_TimeIntegrator_1D
    procedure :: Show => Show_CL_TimeIntegrator_1D
    procedure(TimeStep), deferred :: TimeStep
  end type CL_TimeIntegrator_1D

  abstract interface

    !---------------------------------------------------------------------------
    !> Execution of a single time step

    subroutine TimeStep(this, problem, t, dt, u, M_inv)
      import
      class(CL_TimeIntegrator_1D), intent(inout) :: this
      class(CL_Problem_Scalar_1D), intent(in)     :: problem
      real(RNP),    intent(inout) :: t
      real(RNP),    intent(in)    :: dt              !< step size ∆t
      real(RNP),    intent(inout) :: u(0:,:,:)       !< u(t) → u(t+ ∆t)
      real(RNP),    intent(in)    :: M_inv(:,:,:)
    end subroutine TimeStep

  end interface

  !-----------------------------------------------------------------------------
  !> Base type for providing time integrator options

  type CL_TimeIntegrator_Options_1D
    integer :: impl = 0  !< 0=default:  explicit, 1: implicit, 2: IMEX,
  end type CL_TimeIntegrator_Options_1D

contains

  !=============================================================================
  ! CL_TimeIntegrator_1D: type-bound procedures

  !-----------------------------------------------------------------------------
  !> Initialization of CL_TimeIntegrator_1D object

  subroutine Init_CL_TimeIntegrator_1D(this, opt)
    class(CL_TimeIntegrator_1D),                   intent(inout) :: this
    class(CL_TimeIntegrator_Options_1D), optional, intent(in)    :: opt

    if (present(opt)) then
      this % impl = opt % impl
    end if

  end subroutine Init_CL_TimeIntegrator_1D

  !-----------------------------------------------------------------------------
  !> Output of CL_TimeIntegrator_1D settings

  subroutine Show_CL_TimeIntegrator_1D(this, unit)
    class(CL_TimeIntegrator_1D), intent(in) :: this
    integer, optional, intent(in) :: unit  !< output unit

    integer :: io

    if (present(unit)) then
      io = unit
    else
      io = OUTPUT_UNIT
    end if

    write(io,'(/,A)') 'CL_TimeIntegrator_1D settings'
    write(io,'(A,/)') repeat('=',80)

    ! append further settings in corresponding routines of derived types

  end subroutine Show_CL_TimeIntegrator_1D

  !=============================================================================

end module CL__Time_Integrator__1D
