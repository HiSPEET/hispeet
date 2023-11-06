module CL__SDC__Method__RK__1D

  use, intrinsic :: ISO_Fortran_Env, only: OUTPUT_UNIT

  use Kind_Parameters, only: RNP
  use Constants,       only: ONE, ZERO, HALF, THIRD, TWO

  use Lagrange_Interpolation
  use IMEX_Runge_Kutta_Method
  use Spectral_Deferred_Correction

  use CL__Problem__Scalar__1D
  use CL__Time_Integrator__1D
  use CL__SDC__Method__1D

  implicit none
  private

  public :: CL_SDC_Method_RK_1D
  public :: CL_SDC_Options_RK_1D

  !-----------------------------------------------------------------------------
  !> IMEX Runge-Kutta SDC corrector

  type, extends(CL_SDC_Method_1D) :: CL_SDC_Method_RK_1D
    type(IMEX_RK_Method)   :: imex_rk
    real(RNP), allocatable :: w_rk(:,:,:) !< quadrature weights at RK nodes
    real(RNP), allocatable :: l_rk(:,:,:) !< interpolation weights at RK nodes
  contains
    procedure :: Init_CL_SDC_Method_RK_1D
    procedure :: Show => Show_CL_SDC_Method_RK_1D
    procedure :: CorrectorRHS
    procedure :: CorrectorStep
  end type CL_SDC_Method_RK_1D

  ! overloading the constructor
  interface CL_SDC_Method_RK_1D
    module procedure New_CL_SDC_Method_RK_1D
  end interface

  !-----------------------------------------------------------------------------
  !> Type for providing IMEX Runge-Kutta SDC-corrector options

  type, extends(CL_SDC_Options_1D) :: CL_SDC_Options_RK_1D
    integer :: n_stage = 5  !< number of stages of the predictor
    integer :: method  = 1  !< RK method selector, if more than one exist
  end type CL_SDC_Options_RK_1D

contains

  !=============================================================================
  ! SDC_Corrector_RK: type-bound procedures

  !-----------------------------------------------------------------------------
  !> Constructor for objects of type SDC_Corrector_RK

  function New_CL_SDC_Method_RK_1D(pre_opt, sdc_opt) result(this)
    class(CL_TimeIntegrator_Options_1D), intent(in) :: pre_opt !< predictor options
    class(CL_SDC_Options_RK_1D),        intent(in) :: sdc_opt !< SDC options
    type(CL_SDC_Method_RK_1D) :: this

    call Init_CL_SDC_Method_RK_1D(this, pre_opt, sdc_opt)

  end function New_CL_SDC_Method_RK_1D

  !-----------------------------------------------------------------------------
  !> Initialization of a CL_SDC_Method_1D object

  subroutine Init_CL_SDC_Method_RK_1D(this, pre_opt, sdc_opt)
    class(CL_SDC_Method_RK_1D),         intent(inout) :: this
    class(CL_TimeIntegrator_Options_1D), intent(in) :: pre_opt !< predictor options
    class(CL_SDC_Options_RK_1D),        intent(in) :: sdc_opt !< SDC options

    real(RNP) :: t0, t1, ti
    integer   :: i, j, m

    ! intialize parent type
    call this % Init_CL_SDC_Method_1D(pre_opt, sdc_opt)

    ! initialize RK method
    call this % imex_rk % Init_IMEX_RK_Method(sdc_opt % n_stage, sdc_opt % method)

    this % corrector_name = 'IMEX RK method'

    ! auxiliary components .....................................................

    if (allocated(this % w_rk)) deallocate(this % w_rk)
    if (allocated(this % l_rk)) deallocate(this % l_rk)

    associate(n_sub => sdc_opt % n_sub, imex_rk => this % imex_rk)

      allocate(this % w_rk(0:n_sub, 1:imex_rk%n_stage, 1:n_sub), source = ZERO)
      allocate(this % l_rk(0:n_sub, 1:imex_rk%n_stage, 1:n_sub), source = ZERO)

      subintervals: do m = 1, n_sub

        ! starting and end point of the SDC interval
        t0 = this % t(m-1) !eventually time( ... )
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

  end subroutine Init_CL_SDC_Method_RK_1D

  !-----------------------------------------------------------------------------
  !> Output of SDC_Corrector_RK settings

  subroutine Show_CL_SDC_Method_RK_1D(this, unit)
    class(CL_SDC_Method_RK_1D), intent(in) :: this
    integer, optional, intent(in) :: unit  !< output unit

    integer :: io

    if (present(unit)) then
      io = unit
    else
      io = OUTPUT_UNIT
    end if

    ! show parent settings
    call this % Show_CL_SDC_Method_1D(unit)

    ! show IMEX RK settings
    call this % imex_rk % Show(unit)

  end subroutine Show_CL_SDC_Method_RK_1D

  !-----------------------------------------------------------------------------
  !> Computes F_ex and F_im as defined in the correctors

  subroutine CorrectorRHS(this, problem, t, dt, u, F_ex, F_im)

    class(CL_SDC_Method_RK_1D), intent(in) :: this
    class(CL_Problem_Scalar_1D), intent(in) :: problem
    real(RNP), intent(in)    :: t
    real(RNP), intent(in)    :: dt            !< step size, used with ISD only
    real(RNP), intent(inout) :: u(:,:,:)      !< u
    real(RNP), intent(out)   :: F_ex(:,:,:)   !< explicit RHS for corrector
    real(RNP), intent(out)   :: F_im(:,:,:)   !< implicit RHS for corrector

    select case (this % impl)
    case(0) ! explicit
      F_im = 0
      F_ex(:,:,:) = problem % RHS_Convection(t, u(:,:,:))  & ! t is the wrong time here, but
                  + problem % RHS_Diffusion (t, u(:,:,:))    ! does not matter with periodic BC

    case(1) ! implicit
!      F_im = lambda * u
!      F_ex = 0

    case default ! IMEX
!      F_im =     lambda % re * u
!      F_ex = i * lambda % im * u

    end select

    if (dt > 0) return  ! just to avoid compiler warnings

  end subroutine CorrectorRHS

  !-----------------------------------------------------------------------------
  !> Execution of a single correction step

  subroutine CorrectorStep( this, problem, m, t, u ,M_inv, F      &
                          , F_ex, F_im, F_ex_new, F_im_new )

    class(CL_SDC_Method_RK_1D), intent(inout) :: this
    class(CL_Problem_Scalar_1D), intent(in)    :: problem
    integer  , intent(in)    :: m                    !< current SDC interval index
    real(RNP), intent(in)    :: t(0:)                !< SDC time nodes
    real(RNP), intent(inout) :: u(0:,:,:,0:)         !< uᵏ⁺¹(:m-1),uᵏ→uᵏ⁺¹(m),uᵏ(m+1:)
    real(RNP), intent(in)    :: M_inv(:,:,:)
    real(RNP), intent(in)    :: F(0:,:,:,0:)         !< Fᵏ
    real(RNP), intent(in)    :: F_ex(0:,:,:,0:)      !< F_exᵏ
    real(RNP), intent(in)    :: F_im(0:,:,:,0:)      !< F_imᵏ
    real(RNP), intent(inout) :: F_ex_new(0:,:,:,0:)  !< F_exᵏ⁺¹(0:m-1) → F_exᵏ⁺¹(0:m)
    real(RNP), intent(inout) :: F_im_new(0:,:,:,0:)  !< F_imᵏ⁺¹(0:m-1) → F_imᵏ⁺¹(0:m)

    real(RNP), allocatable :: G_im(:,:,:,:)
    real(RNP), allocatable :: G_ex(:,:,:,:)
    real(RNP), allocatable :: u_rk(:,:,:)
    real(RNP), allocatable :: S(:,:,:)
    real(RNP), allocatable :: S_rk(:,:,:)
    real(RNP), allocatable :: F_ex_rk(:,:,:)
    real(RNP), allocatable :: F_im_rk(:,:,:)
    real(RNP) :: t0, t1, ti, dt, dt_sub

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
             , w_sub   => this % w_sub             &
             , po      => problem % eop % po       &
             , ne      => problem % ne             &
             , nc      => problem % nc             &
             , Mat     => problem % mm             )

      ! workspace ..............................................................

      allocate( G_im(0:po, ne, nc, n_stage) )
      allocate( G_ex(0:po, ne, nc, n_stage) )

      allocate( F_ex_rk(0:po, ne, nc) )
      allocate( F_im_rk(0:po, ne, nc) ) !wasn't there a nice method to allocate a matrix of the same size as u in the time integrator modules?
                                        !that one:    allocate(f, M_inv, mold = u)
      allocate( S(0:po,ne,nc) )
      allocate( S_rk(0:po,ne,nc) )      !is this in fortran a zero-matrix per default or should one multiply it for safety reasons

      ! initialization .........................................................

      t0     = t(m-1)       ! SDC subinterval start
      t1     = t(m)         ! SDC subinterval end
      dt_sub = t1 - t0      ! SDC subinterval length

      dt = t(n_sub) - t(0)  ! full time step length

      ! stage 1 ................................................................

      G_im(:,:,:,1) = F_im_new(:,:,:,m-1) - F_im(:,:,:,m-1)
      G_ex(:,:,:,1) = F_ex_new(:,:,:,m-1) - F_ex(:,:,:,m-1)


      ! stages 2:n_stage .......................................................

      do i = 2, n_stage

        ti = t0 + dt_sub * c(i)

        u_rk    = u(:,:,:,m-1)

        S_rk    = 0

        G_ex(:,:,:,i) = 0
        G_im(:,:,:,i) = 0

        ! SDC quadrature .......................................................

        do j = 0, this % n_sub
          S_rk    = S_rk + dt * M_inv *  F(:,:,:,j) * w_rk(j,i,m)
          G_ex(:,:,:,i) = G_ex(:,:,:,i) - l_rk(j,i,m) * F_ex(:,:,:,j) ! lagrangian interpolation
          G_im(:,:,:,i) = G_im(:,:,:,i) - l_rk(j,i,m) * F_im(:,:,:,j)
        end do

        ! correction ...........................................................

        do j = 1, i-1 ! explicit correction
          u_rk = u_rk + dt_sub * ( a_ex(i,j) * M_inv * G_ex(:,:,:,j) )
        end do

        do j = 1, i   ! implicit correction
          u_rk = u_rk + dt_sub * ( a_im(i,j) * M_inv * G_im(:,:,:,j) )
        end do

        u_rk = u_rk + S_rk

        select case(impl) ! solve implicit equation
        case(0) ! explicit (identity)
          u_rk = u_rk
        case(1) ! implicit
           print *, 'Implicit RK-based SDC doesnt exist yet!'
!          u_rk = u_rk / (ONE - dt_sub * a_im(i,i) * lambda)
        case default ! IMEX
           print *, 'IMEX RK-based SDC doesnt exist yet!'
!          u_rk = u_rk / (ONE - dt_sub * a_im(i,i) * lambda%re)
        end select

        call this % CorrectorRHS(problem, ti, dt_sub, u_rk, F_ex_rk, F_im_rk)


        G_ex(:,:,:,i) = G_ex(:,:,:,i) + F_ex_rk
        G_im(:,:,:,i) = G_im(:,:,:,i) + F_im_rk

      end do

      ! correction .............................................................

      if (all(b_im == a_im(n_stage,:)) .and. all(b_ex == a_ex(n_stage,:))) then

        ! globally stiffly accurate method: already done
        u(:,:,:,m) = u_rk

      else

        S = 0

        do i = 0, n_sub
          S = S + dt * M_inv * F(:,:,:,i) * w_sub(i,m) ! Is that used somewhere?
        end do

        ! assembly .............................................................

        u_rk = u_rk - S_rk
        do i = 1, n_stage
          u_rk = u_rk + dt_sub * ( (b_ex(i) - a_ex(n_stage,i)) * M_inv *G_ex(:,:,:,i) &
                                 + (b_im(i) - a_im(n_stage,i)) * M_inv *G_im(:,:,:,i) )
        end do

      end if

      ! update RHS .............................................................

      call this % CorrectorRHS(problem, t(m), dt, u(:,:,:,m), F_ex_new(:,:,:,m), F_im_new(:,:,:,m))

    end associate

  end subroutine CorrectorStep

  !=============================================================================

end module CL__SDC__Method__RK__1D
