module CL__SDC__Method__ISD1__1D

  use, intrinsic :: ISO_Fortran_Env, only: OUTPUT_UNIT

  use Kind_Parameters, only: RNP
  use Constants,       only: ONE, ZERO
  use Array_Assignments

  use CL__Problem__1D
  use CL__Operator__1D
  use CL__Time_Integrator__1D
  use CL__SDC__Method__1D

  implicit none
  private

  public :: CL_SDC_Method_ISD1_1D
  public :: CL_SDC_Options_ISD1_1D

  !-----------------------------------------------------------------------------
  !> SDC method based on ISD1

  type, extends(CL_SDC_Method_1D) :: CL_SDC_Method_ISD1_1D
    integer :: stab_method !< stabilization method ∊ {0,1}
    integer :: stab_base   !< basis for exponential stabilization
    integer :: stab_power  !< power for polynomial stabilization
    integer :: stab_cutoff !< stabilization cutoff sweep
  contains
    procedure :: Init_CL_SDC_Method_ISD1_1D
    procedure :: Show => Show_CL_SDC_Method_ISD1_1D
    procedure :: GetHighOrderRHS
    procedure :: GetCorrectorRHS
    procedure :: CorrectorStep
  end type CL_SDC_Method_ISD1_1D

  ! overloading the constructor
  interface CL_SDC_Method_ISD1_1D
    module procedure New_CL_SDC_Method_ISD1_1D
  end interface

  !-----------------------------------------------------------------------------
  !> Type for providing SDC ISD1 options

  type, extends(CL_SDC_Options_1D) :: CL_SDC_Options_ISD1_1D
    integer :: stab_method = 0 !< 0/1/3 no/exponential/polynomial stabilization
    integer :: stab_base   = 2 !< basis for exponential stabilization
    integer :: stab_power  = 2 !< power for polynomial stabilization
    integer :: stab_cutoff = huge(1) !< stabilization cutoff sweep
  end type CL_SDC_Options_ISD1_1D

contains

  !=============================================================================
  ! SDC_Corrector_ISD: type-bound procedures

  !-----------------------------------------------------------------------------
  !> Constructor for objects of type CL_SDC_Method_ISD1_1D

  function New_CL_SDC_Method_ISD1_1D(pre_opt, sdc_opt) result(this)
    class(CL_TimeIntegrator_Options_1D), intent(in) :: pre_opt !< predictor opts
    class(CL_SDC_Options_ISD1_1D),       intent(in) :: sdc_opt !< SDC options
    type(CL_SDC_Method_ISD1_1D) :: this

    call Init_CL_SDC_Method_ISD1_1D(this, pre_opt, sdc_opt)

  end function New_CL_SDC_Method_ISD1_1D

  !-----------------------------------------------------------------------------
  !> Initialization of a CL_SDC_Method_ISD1_1D object

  subroutine Init_CL_SDC_Method_ISD1_1D(this, pre_opt, sdc_opt)
    class(CL_SDC_Method_ISD1_1D),        intent(inout) :: this
    class(CL_TimeIntegrator_Options_1D), intent(in) :: pre_opt !< predictor opts
    class(CL_SDC_Options_ISD1_1D),       intent(in) :: sdc_opt !< SDC options

    ! intialize parent type
    call this % Init_CL_SDC_Method_1D(pre_opt, sdc_opt)

    this % corrector_name = 'ISD method of order 1'
    this % stab_method    = sdc_opt % stab_method
    this % stab_base      = sdc_opt % stab_base
    this % stab_power     = sdc_opt % stab_power
    this % stab_cutoff    = sdc_opt % stab_cutoff

  end subroutine Init_CL_SDC_Method_ISD1_1D

  !-----------------------------------------------------------------------------
  !> Output of SDC_Corrector_ISD settings

  subroutine Show_CL_SDC_Method_ISD1_1D(this, unit)
    class(CL_SDC_Method_ISD1_1D), intent(in) :: this
    integer, optional, intent(in) :: unit  !< output unit

    integer :: io

    if (present(unit)) then
      io = unit
    else
      io = OUTPUT_UNIT
    end if

    call this % Show_CL_SDC_Method_1D(unit)

    write(io,'(/)')
    write(io,'(2X,A,T22,I0)') 'stab_method:', this % stab_method
    write(io,'(2X,A,T22,I0)') 'stab_base:'  , this % stab_base
    write(io,'(2X,A,T22,I0)') 'stab_power:' , this % stab_power
    write(io,'(2X,A,T22,I0)') 'stab_cutoff:', this % stab_cutoff

  end subroutine Show_CL_SDC_Method_ISD1_1D

  !-----------------------------------------------------------------------------
  !> RHS for high-order quadrature

  subroutine GetHighOrderRHS(this, cl_problem, cl_operator, t, dt, k, u, F)
    class(CL_SDC_Method_ISD1_1D), intent(in) :: this
    class(CL_Problem_1D),    intent(in)  :: cl_problem
    class(CL_Operator_1D),   intent(in)  :: cl_operator
    real(RNP),               intent(in)  :: t        !< time
    real(RNP),               intent(in)  :: dt       !< time step size
    integer,                 intent(in)  :: k        !< sweep counter
    real(RNP), contiguous,   intent(in)  :: u(:,:,:) !< u(x,t)
    real(RNP), contiguous,   intent(out) :: F(:,:,:) !< RHS

    real(RNP), allocatable, save :: r_c(:,:,:)
    real(RNP), allocatable, save :: r_d(:,:,:)
    real(RNP), allocatable, save :: r_sd(:,:,:)
    real(RNP), allocatable, save :: f_s(:,:,:)
    real(RNP), allocatable, save :: bv(:,:)

    real(RNP), allocatable :: Me_inv(:)
    real(RNP) :: theta

    integer :: e, j

    associate( eop => cl_operator % eop      &
             , dx  => cl_operator % dx       &
             , po  => cl_operator % eop % po &
             , ne  => cl_operator % ne       &
             , nc  => cl_problem  % nc       )

      allocate(r_c , mold = u)
      allocate(r_d , mold = u)
      allocate(r_sd, mold = u)
      allocate(f_s , mold = u)
      allocate(bv(nc,2))
      allocate(Me_inv(0:po), source = ONE/(dx/2 * eop%w))

      call cl_problem % GetBoundaryValues (t, bv)
      call cl_problem % GetConvectionTerm (cl_operator, bv, u, r_c)
      call cl_problem % GetDiffusionTerm  (cl_operator, bv, u, r_d)
      call cl_problem % GetSources        (cl_operator, t , u, f_s)

      if (k < this % stab_cutoff .and. k >= 0) then
        select case(this % stab_method)
        case(1)
          theta = 0.001 * dt / this % stab_base ** k
        case(2)
          theta = dt * (ONE - real(k,RNP)/this%stab_cutoff) ** this%stab_power
        case default
          theta = 0
        end select
      else
        theta = 0
      end if

      if (theta > 0) then
        call cl_problem % GetSDTerm(cl_operator, theta, bv, u, u, r_sd)
      else
        call SetArray(r_sd, ZERO, multi=.true.)
      end if

      do j = 1, nc
      do e = 1, ne
        F(:,e,j) = Me_inv * (r_c(:,e,j) + r_d(:,e,j) + r_sd(:,e,j))+ f_s(:,e,j)
      end do
      end do

      deallocate(r_c, r_d, r_sd, f_s, bv)

    end associate

  end subroutine GetHighOrderRHS

  !---------------------------------------------------------------------------
  !> Computes F_ex and F_im as defined in the corrector
  !>
  !> The initial values `u₀` are used for evaluating nonlinear coefficients,
  !> such as streamline diffusivity.

  subroutine GetCorrectorRHS( this, cl_problem, cl_operator &
                            , t, dt, u_0, u, F_ex, F_im     )

    class(CL_SDC_Method_ISD1_1D), intent(in) :: this
    class(CL_Problem_1D),  intent(in)  :: cl_problem
    class(CL_Operator_1D), intent(in)  :: cl_operator
    real(RNP),             intent(in)  :: t           !< time
    real(RNP),             intent(in)  :: dt          !< time step size
    real(RNP), contiguous, intent(in)  :: u_0 (:,:,:) !< u₀(x,t)
    real(RNP), contiguous, intent(in)  :: u   (:,:,:) !< u(x,t)
    real(RNP), contiguous, intent(out) :: F_ex(:,:,:) !< explicit RHS part
    real(RNP), contiguous, intent(out) :: F_im(:,:,:) !< implicit RHS part

    real(RNP), allocatable, save :: r_c(:,:,:)
    real(RNP), allocatable, save :: r_d(:,:,:)
    real(RNP), allocatable, save :: r_sd(:,:,:)
    real(RNP), allocatable, save :: bv(:,:)

    real(RNP), allocatable :: Me_inv(:)
    integer :: e, k

    associate( nc       => cl_problem  % nc       &
             , eop      => cl_operator % eop      &
             , dx       => cl_operator % dx       &
             , po       => cl_operator % eop % po &
             , ne       => cl_operator % ne       &
             , activity => cl_operator % activity )

      allocate(r_c , mold = u)
      allocate(r_d , mold = u)
      allocate(r_sd, mold = u)
      allocate(bv(nc,2))

      allocate(Me_inv(0:po), source = ONE / (dx/2 * eop%w))

      call cl_problem % GetBoundaryValues (t, bv)
      call cl_problem % GetConvectionTerm (cl_operator, bv, u, r_c)
      call cl_problem % GetDiffusionTerm  (cl_operator, bv, u, r_d)
      call cl_problem % GetSDTerm         (cl_operator, dt, bv, u_0, u, r_sd)

      do k = 1, nc
      do e = 1, ne
        if (activity(e) > 0) then
          F_ex(:,e,k) = Me_inv * r_c(:,e,k)
          F_im(:,e,k) = Me_inv * (r_d(:,e,k) + r_sd(:,e,k))
        else
          F_ex(:,e,k) = 0
          F_im(:,e,k) = 0
        end if
      end do
      end do

      deallocate(r_c, r_d, r_sd, bv)

    end associate

    ! silence compiler warnings
    if (this % n_sub > 0) return

  end subroutine GetCorrectorRHS

  !-----------------------------------------------------------------------------
  !> Execution of a single correction step

  subroutine CorrectorStep( this, cl_problem, cl_operator, m, t, u &
                          , F, F_ex, F_im, F_ex_new, F_im_new, G   )

    class(CL_SDC_Method_ISD1_1D), intent(in) :: this
    class(CL_Problem_1D),         intent(in) :: cl_problem
    class(CL_Operator_1D),        intent(in) :: cl_operator

    real(RNP), intent(in) :: t(0:)
      !< SDC time nodes
    integer, intent(in) :: m
      !< current SDC interval
    real(RNP), intent(inout) :: u(0:,:,:,0:)
      !< solution, in: [uᵏ⁺¹(:m-1),uᵏ(m:)], out: [uᵏ⁺¹(:m),uᵏ(m+1:)]
    real(RNP), contiguous, intent(in) :: F(0:,:,:,0:)
      !< RHS for high-order quadrature, Fᵏ
    real(RNP), contiguous, intent(in) :: F_ex(0:,:,:,0:)
      !< explicit part of low-order RHS, F_exᵏ
    real(RNP), contiguous, intent(in) :: F_im(0:,:,:,0:)
      !< implicit part of low-order RHS, F_imᵏ
    real(RNP), contiguous, intent(inout) :: F_ex_new(0:,:,:,0:)
      !< explicit of new low-order RHS, in: F_exᵏ⁺¹(0:m-1), out: F_exᵏ⁺¹(0:m)
    real(RNP), contiguous, intent(inout) :: F_im_new(0:,:,:,0:)
      !< implicit of new low-order RHS, in: F_imᵏ⁺¹(0:m-1), out: F_imᵏ⁺¹(0:m)
    real(RNP), contiguous, optional, intent(in) :: G(0:,:,:,0:)
      !< FAS defect correction term

    ! auxiliary variables ......................................................

    real(RNP), allocatable, save :: S(:,:,:)
    real(RNP), allocatable, save :: r_c(:,:,:)
    real(RNP), allocatable, save :: r_d(:,:,:)
    real(RNP), allocatable, save :: r_sd(:,:,:)
    real(RNP), allocatable, save :: u_0(:,:,:)
    real(RNP), allocatable, save :: u_1(:,:,:)
    real(RNP), allocatable, save :: u_i(:,:,:)
    real(RNP), allocatable, save :: bv(:,:)

    real(RNP), allocatable :: Me_inv(:)
    real(RNP) :: dt, dt_sub
    integer   :: e, i, k

    associate( n_sub    => this % n_sub           &
             , w_sub    => this % w_sub           &
             , nc       => cl_problem  % nc       &
             , eop      => cl_operator % eop      &
             , dx       => cl_operator % dx       &
             , po       => cl_operator % eop % po &
             , ne       => cl_operator % ne       &
             , activity => cl_operator % activity )

      ! initialization .........................................................

      allocate(S(0:po,ne,nc))
      allocate(r_c , mold = S)
      allocate(r_d , mold = S)
      allocate(r_sd, mold = S)
      allocate(u_0 , mold = S)
      allocate(u_1 , mold = S)
      allocate(u_i , mold = S)
      allocate(bv(nc,2))

      allocate(Me_inv(0:po), source = ONE / (dx/2 * eop%w))

      dt      =  t(n_sub) - t(0)    ! full interval length
      dt_sub  =  t(m)     - t(m-1)  ! subinterval length

      call SetArray(u_0, u(:,:,:,m-1), multi = .true.)
      call SetArray(u_1, u(:,:,:,m  ), multi = .true.)

      ! high-order quadrature ..................................................

      do k = 1, nc
        select case(this % point_set)
        case('RR')
          ! omit left point with Radau-right
          call SetArray(S(:,:,k), ZERO)
        case default
          ! initialize with contribution of left point (i = 0)
          do e = 1, ne
            S(:,e,k) = dt * w_sub(0,m) * F(:,e,k,0)
          end do
        end select
        ! add contribution of remaining points
        do i = 1, n_sub
        do e = 1, ne
          S(:,e,k) = S(:,e,k) + dt * w_sub(i,m) * F(:,e,k,i)
        end do
        end do
      end do

      ! IMEX ISD1 correction ...................................................

      ! intermediate solution
      do k = 1, nc
      do e = 1, ne
        if (activity(e) > 0) then
          u_i(:,e,k) = u_0(:,e,k)                       &
                     + S  (:,e,k)                       &
                     + dt_sub * ( F_ex_new (:,e,k,m-1)  &
                                - F_ex     (:,e,k,m-1)  &
                                - F_im     (:,e,k,m)    )
          if (present(G)) then
            u_i(:,e,k) = u_i(:,e,k) + Me_inv * (G(:,e,k,m) - G(:,e,k,m-1))
          end if
        else
          u_i(:,e,k) = u_0(:,e,k)
        end if
      end do
      end do

      ! implicit diffusion step
      call cl_problem % GetBoundaryValues(t(m), bv)
      call cl_problem % DiffusionSolver( cl_operator, dt_sub, dt_sub, bv  &
                                       , f      = u_i                     &
                                       , u_0    = u_0                     &
                                       , u      = u_1                     &
                                       , method = this % diffusion_method &
                                       , i_max  = this % diffusion_i_max  &
                                       , r_red  = this % diffusion_r_red  &
                                       , r_max  = this % diffusion_r_max  )

      ! update solution and corrector RHS ......................................

      call cl_problem % GetConvectionTerm(cl_operator, bv, u_1, r_c)
      call cl_problem % GetDiffusionTerm (cl_operator, bv, u_1, r_d)
      call cl_problem % GetSDTerm(cl_operator, dt_sub, bv, u_0, u_1, r_sd)

      do k = 1, nc
      do e = 1, ne
        if (activity(e) > 0) then
          u(:,e,k,m) = u_1(:,e,k)
          F_ex_new(:,e,k,m) = Me_inv * r_c(:,e,k)
          F_im_new(:,e,k,m) = Me_inv * (r_d(:,e,k) + r_sd(:,e,k))
        else
          F_ex_new(:,e,k,m) = 0
          F_im_new(:,e,k,m) = 0
        end if
      end do
      end do

      ! clean-up ...............................................................

      deallocate(S, r_c, r_d, r_sd, u_0, u_1, u_i, bv)

    end associate

  end subroutine CorrectorStep

  !=============================================================================

end module CL__SDC__Method__ISD1__1D
