module CL__Time_Integrator__TR__1D

  use, intrinsic :: ISO_Fortran_Env, only: OUTPUT_UNIT

  use Kind_Parameters, only: RNP
  use Constants,       only: ZERO, ONE, THREE, HALF

  use CL__Problem__Scalar__1D
  use CL__Time_Integrator__1D

  implicit none
  private

  public :: CL_TimeIntegrator_TR_1D
  public :: CL_TimeIntegrator_Options_TR_1D

  !-----------------------------------------------------------------------------
  !> IMEX trapezoidal rule for Dahlquist equation

  type, extends(CL_TimeIntegrator_1D) :: CL_TimeIntegrator_TR_1D
  contains
    procedure :: Init_CL_TimeIntegrator_TR_1D
    procedure :: Show => Show_CL_TimeIntegrator_TR_1D
    procedure :: TimeStep
  end type CL_TimeIntegrator_TR_1D

  ! overloading the constructor
  interface CL_TimeIntegrator_TR_1D
    module procedure New_CL_TimeIntegrator_TR_1D
  end interface

  !-----------------------------------------------------------------------------
  !> Type for providing TR time-integrator options (none, so far)

  type, extends(CL_TimeIntegrator_Options_1D) :: CL_TimeIntegrator_Options_TR_1D
  end type CL_TimeIntegrator_Options_TR_1D

contains

  !-----------------------------------------------------------------------------
  !> Constructor for objects of type CL_TimeIntegrator_TR_1D with options

  function New_CL_TimeIntegrator_TR_1D(opt) result(this)
    class(CL_TimeIntegrator_Options_TR_1D), optional, intent(in) :: opt
    type(CL_TimeIntegrator_TR_1D) :: this

    call Init_CL_TimeIntegrator_TR_1D(this, opt)

  end function New_CL_TimeIntegrator_TR_1D

  !-----------------------------------------------------------------------------
  !> Initialization of a Init_CL_TimeIntegrator_TR_1D object

  subroutine Init_CL_TimeIntegrator_TR_1D(this, opt)
    class(CL_TimeIntegrator_TR_1D),                   intent(inout) :: this
    class(CL_TimeIntegrator_Options_TR_1D), optional, intent(in)    :: opt

    ! intialize parent type
    call this % Init_CL_TimeIntegrator_1D(opt)
    this % name = 'IMEX TR method'

  end subroutine Init_CL_TimeIntegrator_TR_1D

  !-----------------------------------------------------------------------------
  !> Output of CL_TimeIntegrator_TR_1D settings

  subroutine Show_CL_TimeIntegrator_TR_1D(this, unit)
    class(CL_TimeIntegrator_TR_1D), intent(in) :: this
    integer, optional, intent(in) :: unit  !< output unit

    integer :: io

    if (present(unit)) then
      io = unit
    else
      io = OUTPUT_UNIT
    end if

    ! show parent settings
    call this % Show_CL_TimeIntegrator_1D(unit)

    write(io,'(2X,A,T15,G0)')  'name:', trim(this % name)
    write(io,'(2X,A,T15,G0)')  'impl:', this % impl

  end subroutine Show_CL_TimeIntegrator_TR_1D

  !-----------------------------------------------------------------------------
  !> Performs an IMEX TR step

  subroutine TimeStep(this, problem, t, dt, u, M_inv)
    class(CL_TimeIntegrator_TR_1D), intent(inout) :: this
    class(CL_Problem_Scalar_1D),    intent(in)    :: problem
    real(RNP), intent(inout) :: t
    real(RNP), intent(in)    :: dt              !< step size ∆t
    real(RNP), intent(inout) :: u(0:,:,:)       !< u(t) → u(t+ ∆t)
    real(RNP), intent(in)    :: M_inv(:,:,:)

    real(RNP), allocatable, dimension(:,:,:), save :: u1, u2, u3
    real(RNP), allocatable, dimension(:,:,:), save :: f

    associate( po => problem % eop % po &
             , ne => problem % ne       &
             , nc => problem % nc       )

      ! workspace
      if (.not. allocated(f)) then
        allocate(f, u1, u2, u3, mold = u)
      end if

      select case(this%impl)

      case(0)

        ! explicit TR ............................................................

        u1 = u

        f  = problem % RHS_Convection(t, u)  &
           + problem % RHS_Diffusion (t, u)

        u2 = u + dt * M_inv * f

        f  = problem % RHS_Convection(t+dt, u2)  &
           + problem % RHS_Diffusion (t+dt, u2)

        u  = u + dt * HALF * M_inv * f

      case(1)

        ! implicit TR ............................................................

        !u = (u + dt * HALF * lambda * u) / (ONE - dt * HALF * lambda)

        print *, "Actually there is no way to use an implicit method at this moment!"
        stop

      case default

        ! IMEX TR ................................................................

        !u1 = u

        !u2 = u  + dt * i * lambda%im * u1
        !u2 = u2 / (ONE - dt * lambda%re)

        !u3 = u  + dt * HALF * (lambda * u1 + i * lambda%im * u2 )
        !u3 = u3 / (ONE - dt * HALF * lambda%re)

        !u = u3

        print *, "Actually there is no way to use an implicit method at this moment!"
        stop

      end select

    end associate

  end subroutine TimeStep

  !=============================================================================

end module CL__Time_Integrator__TR__1D
