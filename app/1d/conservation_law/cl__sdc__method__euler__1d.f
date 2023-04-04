module CL__SDC__Method__Euler__1D

  use, intrinsic :: ISO_Fortran_Env, only: OUTPUT_UNIT

  use Kind_Parameters, only: RNP
  use Constants,       only: ONE, ZERO

  use CL__Problem__Scalar__1D
  use CL__Time_Integrator__1D
  use CL__SDC__Method__1D

  implicit none
  private

  public :: CL_SDC_Method_Euler_1D
  public :: CL_SDC_Options_Euler_1D

  !-----------------------------------------------------------------------------
  !> IMEX Euler SDC ...

  type, extends(CL_SDC_Method_1D) :: CL_SDC_Method_Euler_1D
  contains
    procedure :: Init_CL_SDC_Method_Euler_1D
    procedure :: Show => Show_CL_SDC_Method_Euler_1D
    procedure :: CorrectorRHS
    procedure :: CorrectorStep
  end type CL_SDC_Method_Euler_1D

  ! overloading the constructor
  interface CL_SDC_Method_Euler_1D
    module procedure New_CL_SDC_Method_Euler_1D
  end interface

  !-----------------------------------------------------------------------------
  !> Type for providing SDC-Euler options

  type, extends(CL_SDC_Options_1D) :: CL_SDC_Options_Euler_1D
  end type CL_SDC_Options_Euler_1D

contains

  !=============================================================================
  ! SDC_Corrector_Euler: type-bound procedures

  !-----------------------------------------------------------------------------
  !> Constructor for objects of type SDC_Corrector_Euler

  function New_CL_SDC_Method_Euler_1D(pre_opt, sdc_opt) result(this)
    class(CL_TimeIntegrator_Options_1D), intent(in) :: pre_opt !< predictor opts
    class(CL_SDC_Options_Euler_1D),      intent(in) :: sdc_opt !< SDC options
    type(CL_SDC_Method_Euler_1D) :: this

    call Init_CL_SDC_Method_Euler_1D(this, pre_opt, sdc_opt)

  end function New_CL_SDC_Method_Euler_1D

  !-----------------------------------------------------------------------------
  !> Initialization of a SDC_Corrector_Euler object

  subroutine Init_CL_SDC_Method_Euler_1D(this, pre_opt, sdc_opt)
    class(CL_SDC_Method_Euler_1D),       intent(inout) :: this
    class(CL_TimeIntegrator_Options_1D), intent(in) :: pre_opt !< predictor opts
    class(CL_SDC_Options_Euler_1D),      intent(in) :: sdc_opt !< SDC options

    ! intialize parent type
    call this % Init_CL_SDC_Method_1D(pre_opt, sdc_opt)

    this % corrector_name = 'IMEX Euler method'
  end subroutine Init_CL_SDC_Method_Euler_1D

  !-----------------------------------------------------------------------------
  !> Output of SDC_Corrector_Euler settings

  subroutine Show_CL_SDC_Method_Euler_1D(this, unit)
    class(CL_SDC_Method_Euler_1D), intent(in) :: this
    integer, optional, intent(in) :: unit  !< output unit

    integer :: io

    if (present(unit)) then
      io = unit
    else
      io = OUTPUT_UNIT
    end if

    call this % Show_CL_SDC_Method_1D(unit)

    write(io,'(2X,A,T15,G0)') 'name:', this % corrector_name

  end subroutine Show_CL_SDC_Method_Euler_1D

  !-----------------------------------------------------------------------------
  !> Computes F_ex and F_im as defined in the corrector for a single matrix u

  subroutine CorrectorRHS(this, problem, t, dt, u, F_ex, F_im)

    class(CL_SDC_Method_Euler_1D), intent(in) :: this
    class(CL_Problem_Scalar_1D),   intent(in) :: problem
    real(RNP), intent(in)    :: t
    real(RNP), intent(in)    :: dt           !< step size ∆t
    real(RNP), intent(inout) :: u(:,:,:)     !< u(t) → u(t+ ∆t)
    real(RNP), intent(out)   :: F_ex(:,:,:)  !< explicit RHS for corrector
    real(RNP), intent(out)   :: F_im(:,:,:)  !< implicit RHS for corrector

    select case (this % impl)
    case(0) ! explicit
      F_im = 0
      F_ex(:,:,:) = problem % RHS_Convection(t, u(:,:,:)) &
                  + problem % RHS_Diffusion (t, u(:,:,:))

    case(1) ! implicit
!      F_im = lambda * u
!      F_ex = 0

    case(2) ! IMEX
!      F_im =     lambda % re * u
!      F_ex = i * lambda % im * u

    end select

    if (dt > 0) return  ! just to avoid compiler warnings

  end subroutine CorrectorRHS

  !-----------------------------------------------------------------------------
  !> Execution of a single correction step

  subroutine CorrectorStep( this, problem, m, t, u, M_inv, F &
                          , F_ex, F_im, F_ex_new, F_im_new   )

    class(CL_SDC_Method_Euler_1D), intent(inout) :: this
    class(CL_Problem_Scalar_1D), intent(in)       :: problem
    integer  , intent(in)    :: m                   !< current SDC interval index
    real(RNP), intent(in)    :: t(0:)               !< SDC time nodes
    real(RNP), intent(inout) :: u(0:,:,:,0:)        !< uᵏ⁺¹(:m-1),uᵏ→uᵏ⁺¹(m),uᵏ(m+1:)
    real(RNP), intent(in)    :: M_inv(:,:,:)
    real(RNP), intent(in)    :: F(0:,:,:,0:)        !< Fᵏ
    real(RNP), intent(in)    :: F_ex(0:,:,:,0:)     !< F_exᵏ
    real(RNP), intent(in)    :: F_im(0:,:,:,0:)     !< F_imᵏ
    real(RNP), intent(inout) :: F_ex_new(0:,:,:,0:) !< F_exᵏ⁺¹(0:m-1) → F_exᵏ⁺¹(0:m)
    real(RNP), intent(inout) :: F_im_new(0:,:,:,0:) !< F_imᵏ⁺¹(0:m-1) → F_imᵏ⁺¹(0:m)

    ! auxiliary variables ......................................................

    real(RNP), allocatable   :: S(:,:,:)
    real(RNP), allocatable   :: u1(:,:,:)
    real(RNP), allocatable   :: u2(:,:,:)
    real(RNP)    :: t0, t1, dt_sub
    real(RNP)    :: dt
    integer      :: j

    associate( po   => problem % eop % po &
             , ne   => problem % ne       &
             , nc   => problem % nc       &
             , Mat  => problem % mm       )

      ! workspace ..............................................................

      allocate(S  (0:po, ne, nc))
      allocate(u1 (0:po, ne, nc))
      allocate(u2 (0:po, ne, nc))

      ! initialization .........................................................

      t0     = t(m-1)     ! SDC subinterval start
      t1     = t(m)       ! SDC subinterval end
      dt_sub = t1 - t0    ! SDC subinterval length

      ! SDC quadrature .........................................................

      associate(n_sub => this % n_sub, w_sub => this % w_sub)

        dt = t(n_sub) - t(0) ! original t(n_sub) - t(0)
        S = 0
        do j = 0, n_sub
          S = S + dt * M_inv * F(:,:,:,j) * w_sub(j,m) !Why is here a dt_swub. Isn't that in w_sub = int_(t_m-1)^(t_m) \alpha(t)dt incoporated? Seems like not because solution explodes if dt_sub missing!
        end do
        u1 = u(:,:,:,m-1) + S  ! in christlieb paper 2009a:     u(:,:,:,m-1) would be u(:,:,:,m-1)^k  in Timeste one uses i instead of k as index.
                               ! if one uses m instead of m-1 the solution is usable. With m-1 the corrected solution is worse than the predicted one.

      end associate

      ! correction .............................................................

      select case(this % impl)
      case(0) ! explicit
        u2 = u1 + dt_sub * M_inv * (F_ex_new(:,:,:,m-1) - F_ex(:,:,:,m-1)) !F_new_(...,0) wo wird das denn definiert?
      case(1) ! implicit
         print *, 'Implicit euler-based SDC doesnt exist yet!'
  !      u2 = (u1 - dt_sub * F_im(m)) / (ONE - dt_sub * lambda)
      case(2) ! IMEX
         print *, 'IMEX Euler-based SDC doesnt exist yet!'
  !      u1 = u1 + dt_sub * (F_ex_new(m-1) - F_ex(m-1) - F_im(m))
  !      u2 = u1 / (ONE - dt_sub * lambda%re)
      case(3) ! IMEX partitioned
         print *, 'Partitioned IMEX Euler-based SDC doesnt exist yet!'
  !      u2 = u1 + dt_sub * (F_ex_new(m-1) - F_ex(m-1) + F_im_new(m-1) - F_im(m-1))
  !      u2 = u1 + dt_sub * (i * lambda%im * u2 - F_ex(m) - F_im(m))
  !      u2 = u2 / (ONE - dt_sub * lambda%re)
      end select

      u(:,:,:,m) = u2 ! if uncommented one can see just the predictor-solution which is much better than the one after the corrector step. that is be u^k+1_m

      ! update RHS F_ex_new = f( t_m, u^k+1_m )
      call this % CorrectorRHS(problem, t(m), dt_sub, u(:,:,:,m), F_ex_new(:,:,:,m), F_im_new(:,:,:,m)) ! -> F_ex_new = f( t_(m-1), u^k_m ) but should be f( t_(m), u^k_m ) for computation of u^k_m+1

    end associate

  end subroutine CorrectorStep

  !=============================================================================

end module CL__SDC__Method__Euler__1D
