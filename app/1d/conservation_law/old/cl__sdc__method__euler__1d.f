!> @note  *** WORK IN PROGRESS ***

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
  !> IMEX Euler SDC

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
      F_im = problem % RHS_Diffusion (t, u(:,:,:))
      F_ex = problem % RHS_Convection(t, u(:,:,:))

    end select

    if (dt > 0) return  ! just to avoid compiler warnings

  end subroutine CorrectorRHS

  !-----------------------------------------------------------------------------
  !> Execution of a single correction step

  subroutine CorrectorStep( this, problem, m, t, u, F      &
                          , F_ex, F_im, F_ex_new, F_im_new )

    class(CL_SDC_Method_Euler_1D), intent(inout) :: this
    class(CL_Problem_Scalar_1D), intent(in)       :: problem
    integer  , intent(in)    :: m                   !< current SDC interval index
    real(RNP), intent(in)    :: t(0:)               !< SDC time nodes
    real(RNP), intent(inout) :: u(0:,:,:,0:)        !< uᵏ⁺¹(:m-1),uᵏ→uᵏ⁺¹(m),uᵏ(m+1:)
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

    real(RNP), allocatable, save :: G(:,:,:)
    real(RNP), allocatable, save :: mm_inv(:,:)
    real(RNP) :: bv(2) = 0

    associate( po  => problem % eop % po &
             , ne  => problem % ne       &
             , nc  => problem % nc       &
             , mm  => problem % mm       )

      ! workspace ..............................................................

      allocate(S  (0:po, ne, nc))
      allocate(u1 (0:po, ne, nc))
      allocate(u2 (0:po, ne, nc))

      if (.not. allocated(G)) then
        allocate(G(0:po, ne, nc))
        allocate(mm_inv, source = 1/mm)
      end if

      ! initialization .........................................................

      t0     = t(m-1)     ! SDC subinterval start
      t1     = t(m)       ! SDC subinterval end
      dt_sub = t1 - t0    ! SDC subinterval length

      ! SDC quadrature .........................................................

      associate(n_sub => this % n_sub, w_sub => this % w_sub)

        dt = t(n_sub) - t(0)
        S = 0
        do j = 0, n_sub
          S(:,:,1) = S(:,:,1) + dt * mm_inv * F(:,:,1,j) * w_sub(j,m)
        end do
        u1 = u(:,:,:,m-1) + S

      end associate

      ! correction .............................................................

      select case(this % impl)
      case(0) ! explicit
        u2(:,:,1) = u1(:,:,1) + dt_sub * mm_inv * (F_ex_new(:,:,1,m-1) - F_ex(:,:,1,m-1))
      case(1) ! implicit
         print *, 'Implicit euler-based SDC doesnt exist yet!'
  !      u2 = (u1 - dt_sub * F_im(m)) / (ONE - dt_sub * lambda)
      case(2) ! IMEX
         print *, 'IMEX Euler-based SDC doesnt exist yet!'
  !      u1 = u1 + dt_sub * (F_ex_new(m-1) - F_ex(m-1) - F_im(m))
  !      u2 = u1 / (ONE - dt_sub * lambda%re)
        u1(:,:,1) = 1/dt_sub * mm * u1(:,:,1)
        u1(:,:,1) = u1(:,:,1) + ( F_ex_new(:,:,1,m-1) - F_ex(:,:,1,m-1) - F_im(:,:,1,m) )!( F_im_new(:,:,1,m-1) - F_ex_new(:,:,1,m)- F_im(:,:,1,m) )
        !                       ( F_ex(u^k+1_m)       - F_ex(u^k_m)     - F_im(u^k_m+1) )
        ! Das scheint soweit richtig zu sein, vgl. TeX

        !u1(:,:,1) = 1/dt_sub * mm * ( F_ex_new(:,:,1,m-1) - F_ex(:,:,1,m-1)- F_im(:,:,1,m) )
        !G(:,:,:) = u1 + problem % RHS_Convection(t(m), u1(:,:,:))
        ! t(m) oder t(m-1) oder loop ueber t und u(:,:,1,m-1), u(:,:,1,m) oder u1 oder u2

        call problem % elliptic_op(1) &
                     % HybridSolver( dx      =  problem % dx   &
                                   , lambda  =  ONE/dt_sub     &
                                   , nu      =  problem % nu_c &
                                   , f       =  u1(:,:,1)      &
                                   , bv      =  bv             &
                                   , u       =  u2(:,:,1)      &
                                   , standby = .false.         ) ! Do not change this into .true. or you will get numerical-analysts nightmare!

      case(3) ! IMEX partitioned
         print *, 'Partitioned IMEX Euler-based SDC doesnt exist yet!'
  !      u2 = u1 + dt_sub * (F_ex_new(m-1) - F_ex(m-1) + F_im_new(m-1) - F_im(m-1))
  !      u2 = u1 + dt_sub * (i * lambda%im * u2 - F_ex(m) - F_im(m))
  !      u2 = u2 / (ONE - dt_sub * lambda%re)
      end select

      u(:,:,:,m) = u2

      ! update RHS F_ex_new = f( t_m, u^k+1_m )
      call this % CorrectorRHS(problem, t(m), dt_sub, u(:,:,:,m), F_ex_new(:,:,:,m), F_im_new(:,:,:,m))

    end associate

  end subroutine CorrectorStep

  !=============================================================================

end module CL__SDC__Method__Euler__1D
