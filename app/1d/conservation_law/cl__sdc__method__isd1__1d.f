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
    integer :: n_stages !< number of corrector stages
  contains
    procedure :: Init_CL_SDC_Method_ISD1_1D
    procedure :: Show => Show_CL_SDC_Method_ISD1_1D
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
    integer :: n_stages = 2 !< number of corrector stages {1,2}
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

    this % n_stages = sdc_opt % n_stages

    select case(this % n_stages)
    case(2)
      this % corrector_name = 'ISD method of order 1 with two stages'
    case default
      this % corrector_name = 'ISD method of order 1 with one stage'
    end select

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

  end subroutine Show_CL_SDC_Method_ISD1_1D

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
      allocate(bv(nc,2))

      allocate(Me_inv(0:po), source = ONE / (dx/2 * eop%w))

      call cl_problem % GetBoundaryValues(t, bv)
      call cl_problem % GetConvectionTerm(cl_operator, bv, u, r_c)
      call cl_problem % GetHybridDiffusionTerm &
                            (cl_operator, 'T', dt, bv, u_0, u, r_d)

      do k = 1, nc
      do e = 1, ne
        if (activity(e) > 0) then
          F_ex(:,e,k) = Me_inv * r_c(:,e,k)
          F_im(:,e,k) = Me_inv * r_d(:,e,k)
        else
          F_ex(:,e,k) = 0
          F_im(:,e,k) = 0
        end if
      end do
      end do

      deallocate(r_c, r_d, bv)

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
    real(RNP), allocatable, save :: u_0(:,:,:)
    real(RNP), allocatable, save :: u_i(:,:,:)
    real(RNP), allocatable, save :: u_m(:,:,:)
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
      allocate(u_0 , mold = S)
      allocate(u_i , mold = S)
      allocate(u_m , mold = S)
      allocate(bv(nc,2))

      allocate(Me_inv(0:po), source = ONE / (dx/2 * eop%w))

      dt      =  t(n_sub) - t(0)    ! full interval length
      dt_sub  =  t(m)     - t(m-1)  ! subinterval length

      call SetArray(u_0, u(:,:,:,m-1), multi = .true.)
      call SetArray(u_m, u(:,:,:,m  ), multi = .true.)

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

      ! correction stage 1: u_m = u₁ ...........................................

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
            u_i(:,e,k) = u_i(:,e,k) + Me_inv * G(:,e,k,m)
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
                                       , u      = u_m                     &
                                       , method = this % diffusion_method &
                                       , i_max  = this % diffusion_i_max  &
                                       , r_red  = this % diffusion_r_red  &
                                       , r_max  = this % diffusion_r_max  )

      if (this % n_stages == 2) then

        ! correction stage 2: u_m = u₁ → u₂ ....................................

        if (cl_problem % limiting_scope > 1) then
          if (cl_problem % limiting_method == 1) then
            call cl_problem % MomentLimiter(cl_operator, u_m)
          end if
        end if

        call cl_problem % GetConvectionTerm(cl_operator, bv, u_m, r_c)

        ! intermediate solution
        do k = 1, nc
        do e = 1, ne
          if (activity(e) > 0) then
            ! add approximate ∆F_ex(m) and remove ∆F_ex(m-1)
            u_i(:,e,k) = u_i(:,e,k)                       &
                       + dt_sub * ( Me_inv * r_c(:,e,k)   & ! + F_ex(u₁  (m  ))
                                  - F_ex     (:,e,k,m)    & ! - F_ex(uᵏ  (m  ))
                                  - F_ex_new (:,e,k,m-1)  & ! - F_ex(uᵏ⁺¹(m-1))
                                  + F_ex     (:,e,k,m-1)  ) ! + F_ex(uᵏ  (m-1))
          end if
        end do
        end do

        ! implicit diffusion step
        call cl_problem % DiffusionSolver( cl_operator, dt_sub, dt_sub, bv  &
                                         , f      = u_i                     &
                                         , u_0    = u_0                     &
                                         , u      = u_m                     &
                                         , method = this % diffusion_method &
                                         , i_max  = this % diffusion_i_max  &
                                         , r_red  = this % diffusion_r_red  &
                                         , r_max  = this % diffusion_r_max  )

      end if

      ! update solution and corrector RHS ......................................

      call cl_problem % GetConvectionTerm(cl_operator, bv, u_m, r_c)
      call cl_problem % GetHybridDiffusionTerm &
                            (cl_operator, 'T', dt_sub, bv, u_0, u_m, r_d)

      do k = 1, nc
      do e = 1, ne
        if (activity(e) > 0) then
          u(:,e,k,m) = u_m(:,e,k)
          F_ex_new(:,e,k,m) = Me_inv * r_c(:,e,k)
          F_im_new(:,e,k,m) = Me_inv * r_d(:,e,k)
        else
          F_ex_new(:,e,k,m) = 0
          F_im_new(:,e,k,m) = 0
        end if
      end do
      end do

      ! limiting ...............................................................

      if (cl_problem % limiting_scope > 0) then
        if (cl_problem % limiting_method == 1) then
          call cl_problem % MomentLimiter(cl_operator, u(:,:,:,m))
        end if
      end if

      ! clean-up ...............................................................

      deallocate(S, r_c, r_d, u_0, u_i, u_m, bv)

    end associate

  end subroutine CorrectorStep

  !=============================================================================

end module CL__SDC__Method__ISD1__1D
