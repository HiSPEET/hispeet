module CL__Time_Integrator__RK__1D

  use, intrinsic :: ISO_Fortran_Env, only: OUTPUT_UNIT

  use Kind_Parameters, only: RNP
  use Constants,       only: ZERO, ONE
  use IMEX_Runge_Kutta_Method

  use CL__Problem__Scalar__1D
  use CL__Time_Integrator__1D

  implicit none
  private

  public :: CL_TimeIntegrator_RK_1D
  public :: CL_TimeIntegrator_Options_RK_1D

  !-----------------------------------------------------------------------------
  !> IMEX Runge-Kutta method for Dahlquist equation

  type, extends(CL_TimeIntegrator_1D) :: CL_TimeIntegrator_RK_1D
    type(IMEX_RK_Method) :: imex_rk  !< IMEX Runge-Kutta method
  contains
    procedure :: Init_CL_TimeIntegrator_RK_1D
    procedure :: Show => Show_CL_TimeIntegrator_RK_1D
    procedure :: TimeStep
  end type CL_TimeIntegrator_RK_1D

  ! overloading the constructor
  interface CL_TimeIntegrator_RK_1D
    module procedure New_CL_TimeIntegrator_RK_1D
  end interface

  !-----------------------------------------------------------------------------
  !> Type for providing RK time-integrator options

  type, extends(CL_TimeIntegrator_Options_1D) :: CL_TimeIntegrator_Options_RK_1D
    integer :: n_stage = 4  !< number of stages
    integer :: method  = 1  !< RK method selector, if more than one exist
  end type CL_TimeIntegrator_Options_RK_1D

contains

  !-----------------------------------------------------------------------------
  !> Constructor for objects of type CL_TimeIntegrator_RK_1D with options

  function New_CL_TimeIntegrator_RK_1D(opt) result(this)
    class(CL_TimeIntegrator_Options_RK_1D), optional, intent(in) :: opt
    type(CL_TimeIntegrator_RK_1D) :: this

    call Init_CL_TimeIntegrator_RK_1D(this, opt)

  end function New_CL_TimeIntegrator_RK_1D

  !-----------------------------------------------------------------------------
  !> Initialization of a Init_CL_TimeIntegrator_RK_1D object

  subroutine Init_CL_TimeIntegrator_RK_1D(this, opt)
    class(CL_TimeIntegrator_RK_1D),                   intent(inout) :: this
    class(CL_TimeIntegrator_Options_RK_1D), optional, intent(in)    :: opt

    ! intialize parent type
    call this % Init_CL_TimeIntegrator_1D(opt)
    this % name = 'IMEX Runge-Kutta method'

    ! initialize RK method
    call this % imex_rk % Init_IMEX_RK_Method(opt % n_stage, opt % method)

  end subroutine Init_CL_TimeIntegrator_RK_1D

  !-----------------------------------------------------------------------------
  !> Output of CL_TimeIntegrator_RK_1D settings

  subroutine Show_CL_TimeIntegrator_RK_1D(this, unit)
    class(CL_TimeIntegrator_RK_1D), intent(in) :: this
    integer, optional, intent(in) :: unit  !< output unit

    integer :: io

    if (present(unit)) then
      io = unit
    else
      io = OUTPUT_UNIT
    end if

    ! show parent settings
    call this % Show_CL_TimeIntegrator_1D(unit)

    ! show IMEX RK settings
    call this % imex_rk % Show(unit)

    write(io,'(2X,A,T15,G0)') 'impl:', this % impl

  end subroutine Show_CL_TimeIntegrator_RK_1D

  !-----------------------------------------------------------------------------
  !> Performs an RK time step

  subroutine TimeStep(this, problem, t, dt, u, M_inv)
    class(CL_TimeIntegrator_RK_1D), intent(inout) :: this
    class(CL_Problem_Scalar_1D),    intent(in)    :: problem
    real(RNP),    intent(inout) :: t
    real(RNP),    intent(in)    :: dt              !< step size ∆t
    real(RNP),    intent(inout) :: u(0:,:,:)       !< u(t) → u(t+ ∆t)
    real(RNP),    intent(in)    :: M_inv(:,:,:)

    select case(this%impl)
    case(0)
      !print *, 'Proceding explicit timestep'
      call TimeStep_EX(this, problem, t, dt, u, M_inv)

    case(1)
      print *, 'Implicit integrator doesnt exist yet!'
      !call TimeStep_IM(this, problem, t, dt, u, M_inv)
    case default
      print *, 'IMEX integrator doesnt exist yet!'
      !call TimeStep_IMEX(this, problem, t, dt, u, M_inv)
    end select

  end subroutine TimeStep

  !-----------------------------------------------------------------------------
  !> Performs an IMEX RK step

  subroutine TimeStep_IMEX(this, problem, t, dt, u, M_inv)
    class(CL_TimeIntegrator_RK_1D), intent(inout) :: this
    class(CL_Problem_Scalar_1D),    intent(in)    :: problem
    real(RNP), intent(inout) :: t
    real(RNP), intent(in)    :: dt              !< step size ∆t
    real(RNP), intent(inout) :: u(0:,:,:)       !< u(t) → u(t+ ∆t)
    real(RNP), intent(in)    :: M_inv(:,:,:)

    complex(RNP), allocatable :: u_s(:,:,:,:)  ! stage solutions
    integer :: i, j, a

    associate( a_im    => this % imex_rk % a_im    &
             , a_ex    => this % imex_rk % a_ex    &
             , b_im    => this % imex_rk % b_im    &
             , b_ex    => this % imex_rk % b_ex    &
             , ns      => this % imex_rk % n_stage &
             , po      => problem % eop % po       &
             , ne      => problem % ne             &
             , nc      => problem % nc             )

      ! initialization .........................................................

      !allocate(u_s(ns))

      allocate(u_s(0:po,ne,nc,ns))

      ! stage 1 ................................................................

      u_s(:,:,:,1) = u

      ! stages 2:ns ............................................................

!      do i = 2, ns
!        u_s(i) = u
!        do j = 1, i-1
!          u_s(i) = u_s(i)                                             &
!                 + dt * a_ex(i,j) * (ZERO, ONE) * lambda%im * u_s(j)  &
!                 + dt * a_im(i,j) *               lambda%re * u_s(j)
!        end do
!        u_s(i) = u_s(i) / (ONE - dt * a_im(i,i) * lambda%re)
!      end do

      ! assembly ...............................................................

!      do i = 1, ns
!        u = u + dt * b_ex(i) * (ZERO, ONE) * lambda%im * u_s(i)  &
!              + dt * b_im(i) *               lambda%re * u_s(i)
!      end do

    end associate

  end subroutine TimeStep_IMEX

  !-----------------------------------------------------------------------------
  !> Performs an implicit RK step

  subroutine TimeStep_IM(this, problem, t, dt, u, M_inv)
    class(CL_TimeIntegrator_RK_1D), intent(inout) :: this
    class(CL_Problem_Scalar_1D),    intent(in)    :: problem
    real(RNP), intent(inout) :: t
    real(RNP), intent(in)    :: dt              !< step size ∆t
    real(RNP), intent(inout) :: u(0:,:,:)       !< u(t) → u(t+ ∆t)
    real(RNP), intent(in)    :: M_inv(:,:,:)

    real(RNP), allocatable :: u_s(:,:,:,:)   ! stage solutions
    integer :: i, j, e, k, a

    associate( a_im  => this % imex_rk % a_im    &
             , b_im  => this % imex_rk % b_im    &
             , ns    => this % imex_rk % n_stage &
             , po => problem % eop % po          &
             , ne => problem % ne                &
             , nc => problem % nc                )

      ! initialization .........................................................

      allocate(u_s(0:po,ne,nc,ns))

      ! stage 1 ................................................................

!      u_s(1) = u

      ! stages 2:ns ............................................................

!      do i = 2, ns
!        u_s(i) = u
!        do j = 1, i-1
!          u_s(i) = u_s(i) + dt * a_im(i,j) * lambda * u_s(j)
!        end do
!        u_s(i) = u_s(i) / (ONE - dt * a_im(i,i) * lambda)
!      end do

      ! assembly ...............................................................

!      do i = 1, ns
!        u = u + dt * b_im(i) * lambda * u_s(i)
!      end do

    end associate

  end subroutine TimeStep_IM

  !-----------------------------------------------------------------------------
  !> Performs a single explicit RK step

  subroutine TimeStep_EX(this, problem, t, dt, u, M_inv)
    class(CL_TimeIntegrator_RK_1D), intent(inout) :: this
    class(CL_Problem_Scalar_1D),    intent(in)    :: problem
    real(RNP), intent(inout) :: t
    real(RNP), intent(in)    :: dt              !< step size ∆t
    real(RNP), intent(inout) :: u(0:,:,:)       !< u(t) → u(t+ ∆t)
    real(RNP), intent(in)    :: M_inv(:,:,:)

    real(RNP), allocatable :: u_s(:,:,:,:)         ! stage solutions
    real(RNP), allocatable :: f(:,:,:)
    real(RNP) :: ct
    integer   :: i, j

    associate( a_ex  => this % imex_rk % a_ex    &
             , b_ex  => this % imex_rk % b_ex    &
             , c     => this % imex_rk % c    &
             , ns    => this % imex_rk % n_stage &
             , po => problem % eop % po &
             , ne => problem % ne       &
             , nc => problem % nc       )

      ! workspace ..............................................................

      allocate(u_s(0:po,ne,nc,ns))

      if (.not. allocated(f)) then
        allocate(f, mold = u)
      end if

      ! stage 1 ................................................................

      u_s(:,:,:,1) = u

      ! stages 2:ns ............................................................

      do i = 2, ns

        ct = c(i)*dt !intermediate timesteps

        u_s(:,:,:,i) = u ! initialise with u_t

        do j = 1, i-1
          u_s(:,:,:,i) = u_s(:,:,:,i) + dt*a_ex(i,j)*u_s(:,:,:,j)
        end do
        f  = problem % RHS_Convection( t+ct , u_s(:,:,:,i) )  & ! from tvdrk: t is the wrong time here, but with periodic bc doesn't matter t+ct should be correct
           + problem % RHS_Diffusion ( t+ct , u_s(:,:,:,i) )  ! k_i = M_inv * f ( t+c_i*dt , u_t + dt*sum_j=1^(i-1) ( aij*u_j )
        u_s(:,:,:,i) = M_inv * f
      end do

      ! assembly ...............................................................

      do i = 1, ns
        u = u + dt * b_ex(i) * u_s(:,:,:,i) ! assembly u_(t+1) + u_t + dt*sum_i=1^ns ( b_i * k_i )
      end do

    end associate

  end subroutine TimeStep_EX

  !=============================================================================

end module CL__Time_Integrator__RK__1D
