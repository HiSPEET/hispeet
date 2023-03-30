module CL__Time_Integrator__Euler__1D

  use, intrinsic :: ISO_Fortran_Env, only: OUTPUT_UNIT

  use Kind_Parameters, only: RNP
  use Constants,       only: ZERO, ONE

  use CL__Problem__Scalar__1D
  use CL__Time_Integrator__1D

  implicit none
  private

  public :: CL_TimeIntegrator_Euler_1D
  public :: CL_TimeIntegrator_Options_Euler_1D

  !-----------------------------------------------------------------------------
  !> IMEX Euler method for Dahlquist equation

  type, extends(CL_TimeIntegrator_1D) :: CL_TimeIntegrator_Euler_1D
  contains
    procedure :: Init_CL_TimeIntegrator_Euler_1D
    procedure :: Show => Show_CL_TimeIntegrator_Euler_1D
    procedure :: TimeStep
  end type CL_TimeIntegrator_Euler_1D

  ! overloading the constructor
  interface CL_TimeIntegrator_Euler_1D
    module procedure New_CL_TimeIntegrator_Euler_1D
  end interface

  !-----------------------------------------------------------------------------
  !> Type for providing Euler time-integrator options (none, so far)

  type, extends(CL_TimeIntegrator_Options_1D) :: &
    CL_TimeIntegrator_Options_Euler_1D
  end type CL_TimeIntegrator_Options_Euler_1D

contains

  !-----------------------------------------------------------------------------
  !> Constructor for objects of type CL_TimeIntegrator_Euler_1D with options

  function New_CL_TimeIntegrator_Euler_1D(opt) result(this)
    class(CL_TimeIntegrator_Options_Euler_1D), optional, intent(in) :: opt
    type(CL_TimeIntegrator_Euler_1D) :: this

    call Init_CL_TimeIntegrator_Euler_1D(this, opt)

  end function New_CL_TimeIntegrator_Euler_1D

  !-----------------------------------------------------------------------------
  !> Initialization of a Init_CL_TimeIntegrator_Euler_1D object

  subroutine Init_CL_TimeIntegrator_Euler_1D(this, opt)
    class(CL_TimeIntegrator_Euler_1D),                   intent(inout) :: this
    class(CL_TimeIntegrator_Options_Euler_1D), optional, intent(in)    :: opt

    ! intialize parent type
    call this % Init_CL_TimeIntegrator_1D(opt)
    this % name = 'IMEX Euler method'

  end subroutine Init_CL_TimeIntegrator_Euler_1D

  !-----------------------------------------------------------------------------
  !> Output of CL_TimeIntegrator_Euler_1D settings

  subroutine Show_CL_TimeIntegrator_Euler_1D(this, unit)
    class(CL_TimeIntegrator_Euler_1D), intent(in) :: this
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

  end subroutine Show_CL_TimeIntegrator_Euler_1D

  !-----------------------------------------------------------------------------
  !> Performs an IMEX Euler step: u₁ = u₀ + ∆t (iλᵢ u₀ + λᵣ u₁)

  subroutine TimeStep(this, problem, t, dt, u, M_inv)
    class(CL_TimeIntegrator_Euler_1D), intent(inout) :: this
    class(CL_Problem_Scalar_1D),       intent(in)    :: problem
    real(RNP), intent(inout) :: t
    real(RNP), intent(in)    :: dt              !< step size ∆t
    real(RNP), intent(inout) :: u(0:,:,:)       !< u(t) → u(t+ ∆t)
    real(RNP), intent(in)    :: M_inv(:,:,:)

    real(RNP), allocatable, dimension(:,:,:), save :: f

    associate( po => problem % eop % po &
             , ne => problem % ne       &
             , nc => problem % nc       )

      if (.not. allocated(f)) then
        allocate(f, mold = u)
      end if

      select case(this%impl)
        case(0)
          ! explicit Euler step: u₁ = u₀ + ∆t λ u₀
          f  = problem % RHS_Convection(t, u)  &
             + problem % RHS_Diffusion (t, u)
          u = u + dt * M_inv * f
        case(1)

          ! implicit Euler step: u₁ = u₀ + ∆t λ u₀

          !u = u / (ONE - dt * lambda)
          print *, "Actually there is no way to use an implicit method at this moment!"
          stop

        case default

          ! IMEX Euler step: u₁ = u₀ + ∆t λ u₀

          !u = u + dt * (ZERO, ONE) * lambda%im * u  ! explicit convection with i λᵢ
          !u = u / (ONE - dt * lambda%re)            ! implicit diffusion with λᵣ
          print *, "Actually there is no way to use an implicit method at this moment!"
          stop

      end select

    end associate

  end subroutine TimeStep

  !=============================================================================

end module CL__Time_Integrator__Euler__1D
