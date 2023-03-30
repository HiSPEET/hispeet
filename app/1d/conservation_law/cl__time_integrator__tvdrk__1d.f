module CL__Time_Integrator__TVDRK__1D

  use, intrinsic :: ISO_Fortran_Env, only: OUTPUT_UNIT

  use Kind_Parameters, only: RNP
  use Constants,       only: TWO, THIRD
  use CL__Time_Integrator__1D

  use CL__Problem__Scalar__1D           ! For datatype CL_Problem_Scalar_1D

  implicit none
  private

  public :: CL_TimeIntegrator_TVDRK_1D
  public :: CL_TimeIntegrator_Options_TVDRK_1D

  !-----------------------------------------------------------------------------
  !> IMEX tvdrk method for Dahlquist equation

  type, extends(CL_TimeIntegrator_1D) :: CL_TimeIntegrator_TVDRK_1D
  contains
    procedure :: Init_CL_TimeIntegrator_TVDRK_1D
    procedure :: Show => Show_CL_TimeIntegrator_TVDRK_1D
    procedure :: TimeStep
  end type CL_TimeIntegrator_TVDRK_1D

  ! overloading the constructor
  interface CL_TimeIntegrator_TVDRK_1D
    module procedure New_CL_TimeIntegrator_TVDRK_1D
  end interface

  !-----------------------------------------------------------------------------
  !> Type for providing tvdrk time-integrator options (none, so far)

  type, extends(CL_TimeIntegrator_Options_1D) :: &
    CL_TimeIntegrator_Options_TVDRK_1D
  end type CL_TimeIntegrator_Options_TVDRK_1D

contains

  !-----------------------------------------------------------------------------
  !> Constructor for objects of type CL_TimeIntegrator_TVDRK_1D with options

  function New_CL_TimeIntegrator_TVDRK_1D(opt) result(this)
    class(CL_TimeIntegrator_Options_TVDRK_1D), optional, intent(in) :: opt
    type(CL_TimeIntegrator_TVDRK_1D) :: this

    call Init_CL_TimeIntegrator_TVDRK_1D(this, opt)

  end function New_CL_TimeIntegrator_TVDRK_1D

  !-----------------------------------------------------------------------------
  !> Initialization of a Init_CL_TimeIntegrator_TVDRK_1D object

  subroutine Init_CL_TimeIntegrator_TVDRK_1D(this, opt)
    class(CL_TimeIntegrator_TVDRK_1D),                   intent(inout) :: this
    class(CL_TimeIntegrator_Options_TVDRK_1D), optional, intent(in)    :: opt

    ! intialize parent type
    call this % Init_CL_TimeIntegrator_1D(opt)
    this % name = 'Explicit tvdrk method'

  end subroutine Init_CL_TimeIntegrator_TVDRK_1D

  !-----------------------------------------------------------------------------
  !> Output of CL_TimeIntegrator_TVDRK_1D settings

  subroutine Show_CL_TimeIntegrator_TVDRK_1D(this, unit)
    class(CL_TimeIntegrator_TVDRK_1D), intent(in) :: this
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

  end subroutine Show_CL_TimeIntegrator_TVDRK_1D

  !-----------------------------------------------------------------------------
  !> 3rd order TVD Runge-Kutta step

  subroutine TimeStep(this, problem, t, dt, u, M_inv)
    class(CL_TimeIntegrator_TVDRK_1D), intent(inout) :: this
    class(CL_Problem_Scalar_1D),       intent(in)    :: problem
    real(RNP), intent(inout) :: t
    real(RNP), intent(in)    :: dt              !< step size ∆t
    real(RNP), intent(inout) :: u(0:,:,:)       !< u(t) → u(t+ ∆t)
    real(RNP), intent(in)    :: M_inv(:,:,:)

    real(RNP), allocatable, dimension(:,:,:), save :: f, u0, u1, u2
    real(RNP) :: c0, c1, c2, ct
    integer   :: e, k

    associate( po => problem % eop % po &
             , ne => problem % ne       &
             , nc => problem % nc       )

      ! workspace
      if (.not. allocated(f)) then
        allocate(u0, u1, u2, f, mold = u)
      end if

      ! initial condition
      u0 = u

      ! step 1
      ct = dt
      f  = problem % RHS_Convection(t, u0)  & ! t is the wrong time here, but
         + problem % RHS_Diffusion (t, u0)    ! does not matter with periodic BC
      u1 = u0 + ct * M_inv * f

      ! step 2
      c0 = 0.75_RNP
      c1 = 0.25_RNP
      ct = 0.25_RNP * dt
      f  = problem % RHS_Convection(t, u1)  & ! t is the wrong time here, but
         + problem % RHS_Diffusion (t, u1)    ! does not matter with periodic BC
      u2 = c0 * u0 + c1 * u1 + ct * M_inv * f

      ! step 3
      c0 = THIRD
      c2 = TWO * THIRD
      ct = TWO * THIRD * dt
      f  = problem % RHS_Convection(t, u2)  & ! t is the wrong time here, but
         + problem % RHS_Diffusion (t, u2)    ! does not matter with periodic BC
      u  = c0 * u0 + c2 * u2 + ct * M_inv * f

    end associate

  end subroutine TimeStep

  !=============================================================================

end module CL__Time_Integrator__TVDRK__1D
