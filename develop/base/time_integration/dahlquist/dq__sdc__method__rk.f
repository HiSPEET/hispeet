module DQ__SDC__Method__RK

  use, intrinsic :: ISO_Fortran_Env, only: OUTPUT_UNIT

  use Kind_Parameters, only: RNP
  use Constants,       only: ONE, ZERO, HALF
  use DQ__Time_Integrator
  use DQ__SDC__Method

  use Lagrange_Interpolation
  use IMEX_Runge_Kutta_Method
  use Spectral_Deferred_Correction

  implicit none
  private

  public :: DQ_SDC_Method_RK
  public :: DQ_SDC_Options_RK

  !-----------------------------------------------------------------------------
  !> IMEX Runge-Kutta SDC corrector

  type, extends(DQ_SDC_Method) :: DQ_SDC_Method_RK
    type(IMEX_RK_Method) :: imex_rk
    real(RNP), allocatable :: w_rk(:,:,:) !< quadrature weights at RK nodes
    real(RNP), allocatable :: l_rk(:,:,:) !< interpolation weights at RK nodes
  contains
    procedure :: Init_DQ_SDC_Method_RK
    procedure :: Show => Show_DQ_SDC_Method_RK
    procedure :: CorrectorRHS
    procedure :: CorrectorStep
  end type DQ_SDC_Method_RK

  ! overloading the constructor
  interface DQ_SDC_Method_RK
    module procedure New_DQ_SDC_Method_RK
  end interface

  !-----------------------------------------------------------------------------
  !> Type for providing IMEX Runge-Kutta SDC-corrector options

  type, extends(DQ_SDC_Options) :: DQ_SDC_Options_RK
    integer :: n_stage = 5  !< number of stages of the predictor
    integer :: method  = 1  !< RK method selector, if more than one exist
  end type DQ_SDC_Options_RK

contains

  !=============================================================================
  ! SDC_Corrector_RK: type-bound procedures

  !-----------------------------------------------------------------------------
  !> Constructor for objects of type SDC_Corrector_RK

  function New_DQ_SDC_Method_RK(pre_opt, sdc_opt) result(this)
    class(DQ_TimeIntegratorOptions), intent(in) :: pre_opt !< predictor options
    class(DQ_SDC_Options_RK),        intent(in) :: sdc_opt !< SDC options
    type(DQ_SDC_Method_RK) :: this

    call Init_DQ_SDC_Method_RK(this, pre_opt, sdc_opt)

  end function New_DQ_SDC_Method_RK

  !-----------------------------------------------------------------------------
  !> Initialization of a DQ_SDC_Method object

  subroutine Init_DQ_SDC_Method_RK(this, pre_opt, sdc_opt)
    class(DQ_SDC_Method_RK),         intent(inout) :: this
    class(DQ_TimeIntegratorOptions), intent(in) :: pre_opt !< predictor options
    class(DQ_SDC_Options_RK),        intent(in) :: sdc_opt !< SDC options

    real(RNP) :: t0, t1, ti
    integer   :: i, j, m

    ! intialize parent type
    call this % Init_DQ_SDC_Method(pre_opt, sdc_opt)

    ! initialize RK method
    call this % imex_rk % Init_IMEX_RK_Method(sdc_opt % n_stage, sdc_opt % method)

    this % corrector_name = 'IMEX RK method'

    ! auxiliary components .....................................................

    if (allocated(this % w_rk)) deallocate(this % w_rk)
    if (allocated(this % l_rk)) deallocate(this % l_rk)

    associate(n_sub => this % n_sub, imex_rk => this % imex_rk)

      allocate(this % w_rk(0:n_sub, 1:imex_rk%n_stage, 1:n_sub), source = ZERO)
      allocate(this % l_rk(0:n_sub, 1:imex_rk%n_stage, 1:n_sub), source = ZERO)

      subintervals: do m = 1, n_sub

        ! starting and end point of the SDC interval
        t0 = this % t(m-1)
        t1 = this % t(m)

        rk_points: do i = 1, imex_rk % n_stage

          ! position the current RK node
          ti = t0 + (t1 - t0) * imex_rk % c(i)

          ! weights for numerical integration over interval (t0,ti), where
          ! t0 is the start and ti is the time of RK stage i in subinterval m
          this % w_rk(:,i,m) = this % SubintervalWeights(t0, ti)

          ! Lagrange polynomial to SDC point j at RK node i in subinterval m
          do j = 0, n_sub
            this % l_rk(j,i,m) = LagrangePolynomial(j, this % t, ti)
          end do

        end do rk_points
      end do subintervals

    end associate

  end subroutine Init_DQ_SDC_Method_RK

  !-----------------------------------------------------------------------------
  !> Output of SDC_Corrector_RK settings

  subroutine Show_DQ_SDC_Method_RK(this, unit)
    class(DQ_SDC_Method_RK), intent(in) :: this
    integer, optional, intent(in) :: unit  !< output unit

    integer :: io

    if (present(unit)) then
      io = unit
    else
      io = OUTPUT_UNIT
    end if

    ! show parent settings
    call this % Show_DQ_SDC_Method(unit)

    ! show IMEX RK settings
    call this % imex_rk % Show(unit)

  end subroutine Show_DQ_SDC_Method_RK

  !-----------------------------------------------------------------------------
  !> Computes F_ex and F_im as defined in the correctors

  elemental subroutine CorrectorRHS(this, lambda, dt, u, F_ex, F_im)
    class(DQ_SDC_Method_RK), intent(in) :: this
    complex(RNP), intent(in)  :: lambda !< λ
    real   (RNP), intent(in)  :: dt     !< step size, used with ISD only
    complex(RNP), intent(in)  :: u      !< u
    complex(RNP), intent(out) :: F_ex   !< explicit RHS for corrector
    complex(RNP), intent(out) :: F_im   !< implicit RHS for corrector

    complex(RNP), parameter :: i = (ZERO, ONE)

    select case (this % impl)

    case(0) ! explicit
      F_im = 0
      F_ex = lambda * u

    case(1) ! implicit
      F_im = lambda * u
      F_ex = 0

    case default ! IMEX
      F_im =     lambda % re * u
      F_ex = i * lambda % im * u

    end select

    if (dt > 0) return  ! just to avoid compiler warnings !

  end subroutine CorrectorRHS

  !-----------------------------------------------------------------------------
  !> Execution of a single correction step

  subroutine CorrectorStep( this, lambda, m, t, u , F      &
                          , F_ex, F_im, F_ex_new, F_im_new )

    class(DQ_SDC_Method_RK), intent(inout) :: this
    complex(RNP), intent(in)    :: lambda       !< λ
    integer     , intent(in)    :: m            !< current SDC interval index
    real   (RNP), intent(in)    :: t(0:)        !< SDC time nodes
    complex(RNP), intent(inout) :: u(0:)        !< uᵏ⁺¹(:m-1),uᵏ→uᵏ⁺¹(m),uᵏ(m+1:)
    complex(RNP), intent(in)    :: F(0:)        !< Fᵏ
    complex(RNP), intent(in)    :: F_ex(0:)     !< F_exᵏ
    complex(RNP), intent(in)    :: F_im(0:)     !< F_imᵏ
    complex(RNP), intent(inout) :: F_ex_new(0:) !< F_exᵏ⁺¹(0:m-1) → F_exᵏ⁺¹(0:m)
    complex(RNP), intent(inout) :: F_im_new(0:) !< F_imᵏ⁺¹(0:m-1) → F_imᵏ⁺¹(0:m)

    complex(RNP), allocatable :: G_im(:)
    complex(RNP), allocatable :: G_ex(:)
    complex(RNP) :: S, S_rk, F_ex_rk, F_im_rk, u_rk
    real   (RNP) :: t0, t1, ti, dt, dt_sub

    integer :: i, j

    associate( a_im    => this % imex_rk % a_im    &
             , a_ex    => this % imex_rk % a_ex    &
             , b_im    => this % imex_rk % b_im    &
             , b_ex    => this % imex_rk % b_ex    &
             , c       => this % imex_rk % c       &
             , n_stage => this % imex_rk % n_stage &
             , w_rk    => this % w_rk              &
             , l_rk    => this % l_rk              &
             , impl    => this % impl              &
             , n_sub   => this % n_sub             &
             , w_sub   => this % w_sub             )

      ! initialization .........................................................

      allocate(G_im(n_stage))
      allocate(G_ex(n_stage))

      t0     = t(m-1)       ! SDC subinterval start
      t1     = t(m)         ! SDC subinterval end
      dt_sub = t1 - t0      ! SDC subinterval length

      dt = t(n_sub) - t(0)  ! full time step length


      ! stage 1 ................................................................

      G_im(1) = F_im_new(m-1) - F_im(m-1)
      G_ex(1) = F_ex_new(m-1) - F_ex(m-1)

      ! stages 2:n_stage .......................................................

      do i = 2, n_stage

        ti = t0 + dt_sub * c(i)

        u_rk    = u(m-1)
        S_rk    = 0
        G_ex(i) = 0
        G_im(i) = 0

        do j = 0, this % n_sub
          S_rk    = S_rk + dt * F(j) * w_rk(j,i,m)
          G_ex(i) = G_ex(i) - l_rk(j,i,m) * F_ex(j)
          G_im(i) = G_im(i) - l_rk(j,i,m) * F_im(j)
        end do

        do j = 1, i-1 ! explicit correction
          u_rk = u_rk + dt_sub * ( a_ex(i,j) * G_ex(j) )
        end do

        do j = 1, i   ! implicit correction
          u_rk = u_rk + dt_sub * ( a_im(i,j) * G_im(j) )
        end do

        u_rk = u_rk + S_rk

        select case(impl) ! solve implicit equation

        case(0) ! explicit (identity)
          u_rk = u_rk

        case(1) ! implicit
          u_rk = u_rk / (ONE - dt_sub * a_im(i,i) * lambda)

        case default ! IMEX
          u_rk = u_rk / (ONE - dt_sub * a_im(i,i) * lambda%re)

        end select

        call this % CorrectorRHS(lambda, dt_sub,  u_rk, F_ex_rk, F_im_rk)
        G_ex(i) = G_ex(i) + F_ex_rk
        G_im(i) = G_im(i) + F_im_rk

      end do

      ! correction .............................................................

      if (all(b_im == a_im(n_stage,:)) .and. all(b_ex == a_ex(n_stage,:))) then

        ! globally stiffly accurate method: already done
        u(m) = u_rk

      else

        ! SDC subinterval integral
        S = 0
        do i = 0, n_sub
          S = S + dt * F(i) * w_sub(i,m)
        end do

        ! assembly
        u_rk = u_rk - S_rk
        do i = 1, n_stage
          u_rk = u_rk + dt_sub * ( (b_ex(i) - a_ex(n_stage,i)) * G_ex(i) &
                                 + (b_im(i) - a_im(n_stage,i)) * G_im(i) )
        end do

        u(m) = u_rk + S

      end if

      ! update RHS .............................................................

      call this % CorrectorRHS(lambda, dt, u(m), F_ex_new(m), F_im_new(m))

    end associate

  end subroutine CorrectorStep

  !=============================================================================

end module DQ__SDC__Method__RK
