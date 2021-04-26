!> summary:  SDC method based on IMEX Euler with dual splitting
!> author:   Joerg Stiller
!> date:     2020/05/10
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module CART__ISP_Flow__SDC_Method__Euler
  use, intrinsic :: ISO_Fortran_Env, only: OUTPUT_UNIT
  use Kind_Parameters, only: RNP
  use Constants,       only: ONE, ZERO
  use Array_Assignments
  use XMPI

  use ISP_Flow_Problem

  use CART__ISP_Flow__Boundary_Values
  use CART__ISP_Flow__Diffusion
  use CART__ISP_Flow__Operators
  use CART__ISP_Flow__Pressure
  use CART__ISP_Flow__Projection
  use CART__ISP_Flow__Time_Derivative
  use CART__ISP_Flow__Time_Integrator

  use CART__ISP_Flow__SDC_Method

  implicit none
  private

  public :: Flow_SDC_Method_Euler
  public :: Flow_SDC_Options_Euler

  !-----------------------------------------------------------------------------
  !> IMEX Euler SDC method for incompressible flow

  type, extends(Flow_SDC_Method) :: Flow_SDC_Method_Euler
  contains
    procedure :: Init => Init_Flow_SDC_Method_Euler
    procedure :: Show => Show_Flow_SDC_Method_Euler
    procedure, private :: CorrectorRHS
    procedure, private :: CorrectorStep
  end type Flow_SDC_Method_Euler

  ! overloading the constructor
  interface Flow_SDC_Method_Euler
    module procedure New_Flow_SDC_Method_Euler
  end interface

  !-----------------------------------------------------------------------------
  !> Type providing specific options for the Euler SDC-corrector (none, so far)

  type, extends(Flow_SDC_Options) :: Flow_SDC_Options_Euler
  end type Flow_SDC_Options_Euler

  !-----------------------------------------------------------------------------
  ! Constants

  real(RNP), parameter :: chi = -ONE  !< scaling parameter of F_d3 = χ∇(ν∇⋅v)

contains

  !=============================================================================
  ! Flow_SDC_Method_Euler: type-bound procedures

  !-----------------------------------------------------------------------------
  !> Constructor for objects of type Flow_SDC_Method_Euler

  function New_Flow_SDC_Method_Euler(problem, flow_op, pre_opt, sdc_opt) &
      result(this)
    class(FlowProblem),    target, intent(in) :: problem !< flow problem
    class(FlowOperators),  target, intent(in) :: flow_op !< flow operators
    class(TimeIntegratorOptions),  intent(in) :: pre_opt !< predictor options
    class(Flow_SDC_Options_Euler), intent(in) :: sdc_opt !< SDC options
    type(Flow_SDC_Method_Euler) :: this

    call Init_Flow_SDC_Method_Euler(this, problem, flow_op, pre_opt, sdc_opt)

  end function New_Flow_SDC_Method_Euler

  !-----------------------------------------------------------------------------
  !> Initialization of a Flow_SDC_Method_Euler object

  subroutine Init_Flow_SDC_Method_Euler(this, problem, flow_op, pre_opt, sdc_opt)
    class(Flow_SDC_Method_Euler), intent(inout) :: this
    class(FlowProblem),      target, intent(in) :: problem !< flow problem
    class(FlowOperators),    target, intent(in) :: flow_op !< flow operators
    class(TimeIntegratorOptions),    intent(in) :: pre_opt !< predictor options
    class(Flow_SDC_Options_Euler),   intent(in) :: sdc_opt !< SDC options

    ! intialize parent type
    call this % Init_Flow_SDC_Method(problem, flow_op, pre_opt, sdc_opt)
    this % corrector_name = 'IMEX Euler corrector'

  end subroutine Init_Flow_SDC_Method_Euler

  !-----------------------------------------------------------------------------
  !> Output of Flow_SDC_Method_Euler settings

  subroutine Show_Flow_SDC_Method_Euler(this, unit)
    class(Flow_SDC_Method_Euler), intent(in) :: this
    integer, optional, intent(in) :: unit !< output unit

    integer :: io

    if (present(unit)) then
      io = unit
    else
      io = OUTPUT_UNIT
    end if

    ! show parent settings
    call this % Show_Flow_SDC_Method(unit)

    write(io,*)

  end subroutine Show_Flow_SDC_Method_Euler

  !---------------------------------------------------------------------------
  !> Computes RHS at the SDC nodes
  !>
  !> Input
  !>
  !>   - `t_sdc   (  0:n_sub)`    SDC time nodes
  !>   - `nu_sdc  (…,0:n_sub)`    variable diffusivity
  !>   - `u_sdc   (…,0:n_sub)`    solution variables
  !>   - `m0                 `    time node to start with [0]
  !>
  !> Output
  !>
  !>   - `F_sdc    (…,m0:n_sub)`  RHS evaluated at SDC time nodes
  !>   - `F_sdc_ex (…,m0:n_sub)`  explicit RHS of low-order scheme
  !>   - `F_sdc_im (…,m0:n_sub)`  semi-implicit RHS of low-order scheme
  !>
  !> Notes
  !>
  !>   - in case of constant diffusivity omit `nu_sdc`
  !>   - the low-order RHS is computed only if `F_lo_ex` and `F_lo_im`
  !>     are both present

  subroutine CorrectorRHS( this, t_sdc, nu_sdc, u_sdc &
                         , F_sdc, F_sdc_ex, F_sdc_im  &
                         , m0                         )
    class(Flow_SDC_Method_Euler), intent(in)    :: this
    real(RNP),                    intent(in)    :: t_sdc    (0:)
    real(RNP), target, optional,  intent(in)    :: nu_sdc   (:,:,:,:,:,0:)
    real(RNP),                    intent(in)    :: u_sdc    (:,:,:,:,:,0:)
    real(RNP),                    intent(inout) :: F_sdc    (:,:,:,:,:,0:)
    real(RNP),         optional,  intent(inout) :: F_sdc_ex (:,:,:,:,:,0:)
    real(RNP),         optional,  intent(inout) :: F_sdc_im (:,:,:,:,:,0:)
    integer  ,         optional,  intent(in)    :: m0

    real(RNP), allocatable, save :: F_c (:,:,:,:,:)
    real(RNP), allocatable, save :: F_d3(:,:,:,:,:)
    real(RNP), pointer, contiguous, save :: nu(:,:,:,:,:)
    integer :: m, m0_

    if (present(m0)) then
      m0_ = m0
    else
      m0_ = 0
    end if

    !$omp single
    allocate(F_c , mold = u_sdc(:,:,:,:,:,0))
    allocate(F_d3, mold = u_sdc(:,:,:,:,:,0))
    nu => null()
    !$omp end single

    do m = m0_, ubound(u_sdc,6)

      ! RHS for high-order sub-integrals .......................................

      if (present(nu_sdc)) nu => nu_sdc (:,:,:,:,:,m)

      call TimeDerivative( this % problem, this % flow_op &
                         , t   = t_sdc(m)                 &
                         , nu  = nu                       &
                         , chi = chi                      &
                         , u   = u_sdc(:,:,:,:,:,m)       &
                         , F   = F_sdc(:,:,:,:,:,m)       &
                         , F_c = F_c                      )

      ! RHS for low-order scheme ...............................................

      if (present(F_sdc_ex) .and. present(F_sdc_im)) then

        ! semi-implicit diffusion using ν(t(m-1))
        if (present(nu_sdc)) nu => nu_sdc(:,:,:,:,:,max(m-1,0))

        call TimeDerivative( this % problem, this % flow_op &
                           , t    = t_sdc(m)                &
                           , nu   = nu                      &
                           , chi  = chi                     &
                           , u    = u_sdc    (:,:,:,:,:,m)  &
                           , F_d1 = F_sdc_im (:,:,:,:,:,m)  &
                           , F_d2 = F_sdc_ex (:,:,:,:,:,m)  &
                           , F_d3 = F_d3                    )

        call MergeArrays(ONE, F_sdc_ex(:,:,:,:,:,m), ONE, F_c , multi=.true.)
        call MergeArrays(ONE, F_sdc_ex(:,:,:,:,:,m), ONE, F_d3, multi=.true.)

      end if

    end do

    !$omp single
    deallocate(F_c, F_d3)
    !$omp end single

  end subroutine CorrectorRHS

  !-----------------------------------------------------------------------------
  !> Computes F_ex and F_im as defined in the corrector and F for subintegrals

  subroutine CorrectorStep( this, m, t_sdc, nu_sdc, u_sdc &
                          , F_sdc, F_sdc_ex, F_sdc_im     &
                          , F_sdc_ex_new, F_sdc_im_new    )

    class(Flow_SDC_Method_Euler), intent(inout) :: this
    integer  ,                    intent(in)    :: m
    real(RNP),                    intent(in)    :: t_sdc        (0:)
    real(RNP), target, optional,  intent(inout) :: nu_sdc       (:,:,:,:,:,0:)
    real(RNP),                    intent(inout) :: u_sdc        (:,:,:,:,:,0:)
    real(RNP),                    intent(in)    :: F_sdc        (:,:,:,:,:,0:)
    real(RNP),                    intent(in)    :: F_sdc_ex     (:,:,:,:,:,0:)
    real(RNP),                    intent(in)    :: F_sdc_im     (:,:,:,:,:,0:)
    real(RNP),                    intent(inout) :: F_sdc_ex_new (:,:,:,:,:,0:)
    real(RNP),                    intent(inout) :: F_sdc_im_new (:,:,:,:,:,0:)

    ! local variables  .........................................................

    real(RNP), allocatable, save :: u_i  (:,:,:,:,:) ! intermediate solution
    real(RNP), allocatable, save :: F_im (:,:,:,:,:) ! aproximate F_im
    real(RNP), allocatable, save :: F_d3 (:,:,:,:,:) ! aproximate F_d3
    real(RNP), allocatable, save :: w    (:,:,:,:,:) ! workspace for u
    real(RNP), allocatable, save :: dp   (:,:,:,:)   ! pressure correction

    real(RNP), contiguous, pointer, save :: nu_0(:,:,:,:,:) ! νᵏ⁺¹(t(m-1))
    real(RNP), contiguous, pointer, save :: nu  (:,:,:,:,:) ! νᵏ⁺¹(t(m))

    real(RNP) :: t, dt, tau
    integer   :: i

    associate( n_sub   => this % n_sub            &
             , problem => this % problem          &
             , flow_op => this % flow_op          &
             , mesh    => this % flow_op % mesh   &
             , u_0     => u_sdc(:,:,:,:,:,m-1)    &
             , u       => u_sdc(:,:,:,:,:,m)      &
             , p       => u_sdc(:,:,:,:,4,m)      )

      ! initialization .........................................................

      !$omp single

      allocate(u_i , mold = u)
      allocate(F_im, mold = u)
      allocate(F_d3, mold = u)
      allocate(w   , mold = u)
      allocate(dp  , mold = p)

      if (present(nu_sdc)) then
        nu_0 => nu_sdc(:,:,:,:,:,m-1)
        nu   => nu_sdc(:,:,:,:,:,m-1)
      else
        nu_0 => null()
        nu   => null()
      end if

      t  = t_sdc(m)
      dt = t_sdc(m) - t_sdc(m-1)

      !$omp end single

      ! boundary conditions ....................................................

      call GetBoundaryValues(problem, mesh, flow_op%bv_x, t, flow_op%bv_u)

      ! explicit + extrapolated diffusive parts ................................

      ! compute F_im = F_d1(ν₀ᵏ⁺¹,uᵏ) and F_d3(ν₀ᵏ⁺¹,uᵏ)
      call TimeDerivative( problem, flow_op, t, u = u, nu = nu_0, chi = chi &
                         , F_d1 = F_im, F_d3 = F_d3                         )

      ! u' = u₀ + ∆t [ F^ex(t₀)ᵏ⁺¹    - F^ex(t₀)ᵏ
      !              + F^im(ν₀ᵏ⁺¹,uᵏ) - F^im(ν₀ᵏ,uᵏ)
      !              + F_d3(ν₀ᵏ⁺¹,uᵏ) ]
      call SetArray(u_i, u_0, multi=.true.)
      do i = 1, problem % nc
        if (i == 4) cycle
        call MergeArrays(ONE, u_i(:,:,:,:,i),  dt, F_sdc_ex_new (:,:,:,:,i,m-1))
        call MergeArrays(ONE, u_i(:,:,:,:,i),  dt, F_im         (:,:,:,:,i))
        call MergeArrays(ONE, u_i(:,:,:,:,i), -dt, F_sdc_ex     (:,:,:,:,i,m-1))
        call MergeArrays(ONE, u_i(:,:,:,:,i), -dt, F_sdc_im     (:,:,:,:,i,m))
        if (i < 4) then
          call MergeArrays(ONE, u_i(:,:,:,:,i), dt, F_d3(:,:,:,:,i))
        end if
      end do

      ! u' = u' + Sᵏ
      do i = 0, n_sub
        tau = (t_sdc(n_sub) - t_sdc(0)) * this % w_sub(i,m)
        call MergeArrays(ONE, u_i, tau, F_sdc(:,:,:,:,:,i), multi=.true.)
      end do

      ! projection .............................................................

      ! solve for p = p"
      call PressureSolver(problem, flow_op, dt, u_i, p, w)

      ! v" = v' - 1/∆t ∇p" - J(v")
      call ProjectionStep(problem, flow_op, dt, p, u_i, w)

      ! diffusion ...............................................................

      ! starting values and RHS for diffusion
      do i = 1, problem % nc
        if (i == 4) cycle
        call SetArray(u(:,:,:,:,i), u_i(:,:,:,:,i))
        call MergeArrays(ONE, u_i(:,:,:,:,i), -dt, F_im(:,:,:,:,i))
        if (i < 4) then
          call MergeArrays(ONE, u_i(:,:,:,:,i), -dt, F_d3(:,:,:,:,i))
        end if
      end do

      ! solve u'''/∆t - ν(t₀)ᵏ⁺¹ u''' = u_i/∆t
      call DiffusionStep(problem, flow_op, dt, f=u_i, u=u, w=w, nu=nu_0)

      ! final projection .......................................................

      if (this % project > 0) then

        ! p = p" + ∆p  with ∆p such that  ∇²∆p = ∆t ∇⋅(v''' - v) = 0
        call SetArray(dp, ZERO)
        call PressureSolver(problem, flow_op, dt, u, dp, w)
        call MergeArrays(ONE, p, ONE, dp)

        ! v = v''' - 1/∆t ∇p - J(v)
        call ProjectionStep(problem, flow_op, dt, dp, u, w)

      end if

      ! RHS and diffusivity updates ............................................

      call TimeDerivative( problem, flow_op, t, u = u, nu = nu_0, chi = chi &
                         , F_c  = F_sdc_ex_new(:,:,:,:,:,m)                 &
                         , F_d1 = F_sdc_im_new(:,:,:,:,:,m)                 &
                         , F_d2 = w                                         &
                         , F_d3 = F_d3                                      )

      do i = 1, problem % nc
        if (i == 4) cycle
        call MergeArrays(ONE, F_sdc_ex_new(:,:,:,:,i,m), ONE, w(:,:,:,:,i))
        if (i >  4) cycle
        call MergeArrays(ONE, F_sdc_ex_new(:,:,:,:,i,m), ONE, F_d3(:,:,:,:,i))
      end do

      if (associated(nu)) then
        call problem % GetDiffusivity(flow_op%x, t, u, nu)
      end if

      ! clean-up ...............................................................

      !$omp barrier
      !$omp single
      deallocate(u_i )
      deallocate(F_im)
      deallocate(F_d3)
      deallocate(w   )
      deallocate(dp  )
      !$omp end single

    end associate

  end subroutine CorrectorStep

  !=============================================================================

end module CART__ISP_Flow__SDC_Method__Euler
