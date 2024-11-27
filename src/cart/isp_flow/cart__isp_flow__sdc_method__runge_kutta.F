!> summary:  SDC method based on IMEX Runge-Kutta with dual splitting
!> author:   Joerg Stiller
!> date:     2021/04/10
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!> @todo
!>   - integration of SGS model
!>   - revise naming and style
!>   - improve data reuse
!===============================================================================

module CART__ISP_Flow__SDC_Method__Runge_Kutta
  use, intrinsic :: ISO_Fortran_Env, only: OUTPUT_UNIT
  use Kind_Parameters, only: RNP
  use Constants,       only: ONE, ZERO
  use Lagrange_Interpolation
  use IMEX_Runge_Kutta_Method
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

  use IEEE_Arithmetic

  implicit none
  private

  public :: Flow_SDC_Method_RungeKutta
  public :: Flow_SDC_Options_RungeKutta

  !-----------------------------------------------------------------------------
  !> IMEX Runge-Kutta SDC corrector for incompressible flow

  type, extends(Flow_SDC_Method) :: Flow_SDC_Method_RungeKutta

    type(IMEX_RK_Method) :: imex_rk !< IMEX Runge-Kutta method

    integer :: diffusion !< method for computing the diffusive RHS (F_d1)

    real(RNP), allocatable :: w_rk(:,:,:) !< quadrature weights at RK nodes
    real(RNP), allocatable :: l_rk(:,:,:) !< interpolation weights at RK nodes

  contains

    procedure :: Init => Init_Flow_SDC_Method_RungeKutta
    procedure :: Show => Show_Flow_SDC_Method_RungeKutta
    procedure :: CorrectorRHS
    procedure :: CorrectorStep

  end type Flow_SDC_Method_RungeKutta

  ! overloading the constructor
  interface Flow_SDC_Method_RungeKutta
    module procedure New_Flow_SDC_Method_RungeKutta
  end interface

  !-----------------------------------------------------------------------------
  !> Type for providing Runge-Kutta SDC-corrector options (none, so far)

  type, extends(Flow_SDC_Options) :: Flow_SDC_Options_RungeKutta
    integer :: n_stage   = 3 !< number of stages
    integer :: method    = 1 !< RK method selector, if more than one exist
  end type Flow_SDC_Options_RungeKutta

  !-----------------------------------------------------------------------------
  ! Constants

  real(RNP), parameter :: chi = -ONE  !< scaling parameter of F_d3 = χ∇(ν∇⋅v)

contains

  !=============================================================================
  ! Flow_SDC_Method_RungeKutta: type-bound procedures

  !-----------------------------------------------------------------------------
  !> Constructor for objects of type Flow_SDC_Method_RungeKutta

  function New_Flow_SDC_Method_RungeKutta(problem, flow_op, pre_opt, sdc_opt) &
      result(this)
    class(FlowProblem),         target, intent(in) :: problem !< flow problem
    class(FlowOperators),       target, intent(in) :: flow_op !< flow operators
    class(TimeIntegratorOptions),       intent(in) :: pre_opt !< predictor options
    class(Flow_SDC_Options_RungeKutta), intent(in) :: sdc_opt !< SDC options
    type(Flow_SDC_Method_RungeKutta) :: this

    call Init_Flow_SDC_Method_RungeKutta(this, problem, flow_op, pre_opt, sdc_opt)

  end function New_Flow_SDC_Method_RungeKutta

  !-----------------------------------------------------------------------------
  !> Initialization of a Flow_SDC_Method_RungeKutta object

  subroutine Init_Flow_SDC_Method_RungeKutta( this, problem, flow_op &
                                            , pre_opt, sdc_opt       )

    class(Flow_SDC_Method_RungeKutta), intent(inout) :: this
    class(FlowProblem),         target, intent(in) :: problem !< flow problem
    class(FlowOperators),       target, intent(in) :: flow_op !< flow operators
    class(TimeIntegratorOptions),       intent(in) :: pre_opt !< predictor options
    class(Flow_SDC_Options_RungeKutta), intent(in) :: sdc_opt !< SDC options

    real(RNP) :: t0, t1, ti
    integer   :: i, j, m

    ! parent type ..............................................................

    call this % Init_Flow_SDC_Method(problem, flow_op, pre_opt, sdc_opt)
    this % corrector_name = 'IMEX Runge-Kutta corrector'

    ! Runge-Kutta method .......................................................

    call this % imex_rk % Init_IMEX_RK_Method( sdc_opt % n_stage &
                                             , sdc_opt % method )

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

  end subroutine Init_Flow_SDC_Method_RungeKutta

  !-----------------------------------------------------------------------------
  !> Output of Flow_SDC_Method_RungeKutta settings

  subroutine Show_Flow_SDC_Method_RungeKutta(this, unit)
    class(Flow_SDC_Method_RungeKutta), intent(in) :: this
    integer, optional, intent(in) :: unit  !< output unit

    integer :: io

    if (present(unit)) then
      io = unit
    else
      io = OUTPUT_UNIT
    end if

    ! show parent settings
    call this % Show_Flow_SDC_Method(unit)

    ! show IMEX RK settings
    call this % imex_rk % Show(unit)

    write(io,*)

  end subroutine Show_Flow_SDC_Method_RungeKutta

  !---------------------------------------------------------------------------
  !> Computes RHS at the SDC nodes
  !>
  !> Input
  !>
  !>   - `t_sdc   (  0:n_sub)`     SDC time nodes
  !>   - `nu_sdc  (…,0:n_sub)`     variable diffusivity
  !>   - `u_sdc   (…,0:n_sub)`     solution variables
  !>   - `m0                 `     time node to start with [0]
  !>
  !> Output
  !>
  !>   - `F_sdc    (…,m,0:n_sub)`  RHS evaluated at SDC time nodes
  !>   - `F_sdc_ex (…,m,0:n_sub)`  explicit RHS of low-order scheme
  !>   - `F_sdc_im (…,m,0:n_sub)`  semi-implicit RHS of low-order scheme
  !>
  !> Notes
  !>
  !>   - in case of constant diffusivity omit `nu_sdc`
  !>   - the low-order RHS is computed only if `F_lo_ex` and `F_lo_im`
  !>     are both present

  subroutine CorrectorRHS( this, t_sdc, nu_sdc, u_sdc &
                         , F_sdc, F_sdc_ex, F_sdc_im  &
                         , m0                         )
    class(Flow_SDC_Method_RungeKutta), intent(in) :: this
    real(RNP),                    intent(in)    :: t_sdc    (0:)
    real(RNP), target, optional,  intent(in)    :: nu_sdc   (:,:,:,:,:,0:)
    real(RNP),                    intent(in)    :: u_sdc    (:,:,:,:,:,0:)
    real(RNP),                    intent(inout) :: F_sdc    (:,:,:,:,:,0:)
    real(RNP),         optional,  intent(inout) :: F_sdc_ex (:,:,:,:,:,0:)
    real(RNP),         optional,  intent(inout) :: F_sdc_im (:,:,:,:,:,0:)
    integer  ,         optional,  intent(in)    :: m0

    real(RNP), allocatable, save :: F_d2(:,:,:,:,:)
    real(RNP), allocatable, save :: F_d3(:,:,:,:,:)
    real(RNP), pointer, contiguous, save :: nu(:,:,:,:,:)
    integer :: m, m0_

    if (present(m0)) then
      m0_ = m0
    else
      m0_ = 0
    end if

    !$omp single
    allocate(F_d2, mold = u_sdc(:,:,:,:,:,0))
    allocate(F_d3, mold = u_sdc(:,:,:,:,:,0))
    nu => null()
    !$omp end single

    do m = m0_, ubound(u_sdc,6)

      if (present(nu_sdc)) nu => nu_sdc (:,:,:,:,:,m)

      call TimeDerivative( this % problem, this % flow_op &
                         , t    = t_sdc(m)                &
                         , nu   = nu                      &
                         , chi  = chi                     &
                         , u    = u_sdc(:,:,:,:,:,m)      &
                         , F    = F_sdc(:,:,:,:,:,m)      &
                         , F_c  = F_sdc_ex (:,:,:,:,:,m)  &
                         , F_d1 = F_sdc_im (:,:,:,:,:,m)  &
                         , F_d2 = F_d2                    &
                         , F_d3 = F_d3                    )

      call MergeArrays(ONE, F_sdc_ex(:,:,:,:,:,m), ONE, F_d2, multi=.true.)
      call MergeArrays(ONE, F_sdc_ex(:,:,:,:,:,m), ONE, F_d3, multi=.true.)

    end do

    !$omp single
    deallocate(F_d2, F_d3)
    !$omp end single

  end subroutine CorrectorRHS

  !-----------------------------------------------------------------------------
  !> Execution of a single correction step
  !>
  !> RK-SDC corrector notes:
  !>   -  subinterval integral S
  !>      *  computed on demand in RK stage i > 1
  !>      *  realized in identical way as in Euler corrector
  !>      *  FSAL property → S for assembly is identical to S in last stage
  !>   -  G = F^{k+1} - F^k takes over the role of F
  !>   -  for stage i
  !>      *  G_rk_im(:,:,:,:,:,i) ≈ F_im^{k+1} - F_im^{k}
  !>      *  G_rk_ex(:,:,:,:,:,i) ≈ F_ex^{k+1} - F_ex^{k}
  !>   -  dimensions 1 to 5 of G_rk_im/ex identical to those of u_sdc
  !>   -  RK stage 1 coincides with SDC node m-1:
  !>      *  G_rk_im(*,1) = F_sdc_im_new(*,m-1) - F_sdc_im(*,m-1)
  !>      *  G_rk_ex(*,1) = F_sdc_ex_new(*,m-1) - F_sdc_ex(*,m-1)
  !>   -  RK stages i > 1
  !>      *  initialization of G(*,i) with interpolated F_sdc_
  !>           +  G_rk_im(*,i) = -F_rk_im ← interpolated from F_sdc_im
  !>           +  G_rk_ex(*,i) = -F_rk_ex ← interpolated from F_sdc_ex
  !>           +  no interpolation required in last stage if FSAL
  !>      *  proceed within the stage as usual
  !>      *  at end of the stage compute
  !>           +  F_rk_im = F_d1
  !>           +  F_rk_ex = F_c + F_d2 + F_d3
  !>      *  complete G(*,i)
  !>           +  G_rk_im(*,i) = G_rk_im(*,i) + F_rk_im
  !>           +  G_rk_ex(*,i) = G_rk_ex(*,i) + F_rk_ex
  !>      *  if b == a(n_stage,:), i.e. no assembly
  !>           +  F_sdc_im_new(*,m) = F_rk_im
  !>           +  F_sdc_ex_new(*,m) = F_rk_ex
  !>   -  if assembly is required
  !>      *  compute solution
  !>      *  recompute F_rk_im/ex
  !>   -  return
  !>      *  u_sdc        (*,m) = u
  !>      *  nu_sdc       (*,m) = nu           (if required)
  !>      *  F_sdc_im_new (*,m) = F_rk_im
  !>      *  F_sdc_ex_new (*,m) = F_rk_ex

  subroutine CorrectorStep( this, m, t_sdc, nu_sdc, u_sdc &
                          , F_sdc, F_sdc_ex, F_sdc_im     &
                          , F_sdc_ex_new, F_sdc_im_new    )

    class(Flow_SDC_Method_RungeKutta), intent(inout) :: this
    integer  ,                   intent(in)    :: m
    real(RNP),                   intent(in)    :: t_sdc        (0:)
    real(RNP), target, optional, intent(inout) :: nu_sdc       (:,:,:,:,:,0:)
    real(RNP),                   intent(inout) :: u_sdc        (:,:,:,:,:,0:)
    real(RNP),                   intent(in)    :: F_sdc        (:,:,:,:,:,0:)
    real(RNP),                   intent(in)    :: F_sdc_ex     (:,:,:,:,:,0:)
    real(RNP),                   intent(in)    :: F_sdc_im     (:,:,:,:,:,0:)
    real(RNP),                   intent(inout) :: F_sdc_ex_new (:,:,:,:,:,0:)
    real(RNP),                   intent(inout) :: F_sdc_im_new (:,:,:,:,:,0:)

    ! local variables  .........................................................

    real(RNP), dimension(:,:,:,:,:)  , allocatable, save :: S_i
    real(RNP), dimension(:,:,:,:,:,:), allocatable, save :: G_ex
    real(RNP), dimension(:,:,:,:,:,:), allocatable, save :: G_im
    real(RNP), dimension(:,:,:,:,:,:), allocatable, save :: F_d3

    real(RNP), dimension(:,:,:,:,:), contiguous, pointer, save :: nu_i

    logical :: fsal

    ! auxiliary
    real(RNP) :: t_i, t_m, tau, dt
 !compare to EU Corrector: dt -> dt_m, dt_sdc -> dt
 real(RNP) :: dt_sdc
    integer   :: nc, ne, np
    integer   :: i, l

    associate( problem => this % problem           &
             , flow_op => this % flow_op           &
             , mesh    => this % flow_op % mesh    &
             , x       => this % flow_op % x       &
             , u_1     => u_sdc(:,:,:,:,:,m-1)     & ! u₁ = u_{m-1}^{k+1}
             , u_i     => u_sdc(:,:,:,:,:,m)       & ! uᵢ = u_{m  }^{k+1}
             , psi     => u_sdc(:,:,:,:,4,m)       &
             , a_im    => this % imex_rk % a_im    &
             , a_ex    => this % imex_rk % a_ex    &
             , c       => this % imex_rk % c       &
             , b_im    => this % imex_rk % b_im    &
             , b_ex    => this % imex_rk % b_ex    &
             , n_stage => this % imex_rk % n_stage &
             , w_rk    => this % w_rk              &
             , l_rk    => this % l_rk              &
             , w_sdc   => this % w_nn              &
             , eop     => this % flow_op % eop_u   )

      ! initialization .........................................................

      ! time step
      dt = t_sdc(m) - t_sdc(m-1)

      ! dimensions
      np = size(u_1,1)  ! elements points per direction
      ne = size(u_1,4)  ! number of local elements
      nc = size(u_1,5)  ! number of solution components

      fsal = c(n_stage) == ONE

      !$omp single
      allocate(S_i, mold = u_1)
      allocate(G_ex(np, np, np, ne, nc, n_stage))
      allocate(G_im, mold = G_ex)
      allocate(F_d3, mold = G_ex)
      if (present(nu_sdc)) then
        nu_i => nu_sdc(:,:,:,:,:,m)
      else
        nu_i => null()
      end if
      !$omp end single

      ! stage 1 ................................................................

      t_i = t_sdc(m-1)

      G_im(:,:,:,:,:,1) = F_sdc_im_new(:,:,:,:,:,m-1) - F_sdc_im(:,:,:,:,:,m-1)
      G_ex(:,:,:,:,:,1) = F_sdc_ex_new(:,:,:,:,:,m-1) - F_sdc_ex(:,:,:,:,:,m-1)

      if (associated(nu_i)) then
        call problem % GetDiffusivity(flow_op%x, t_i, u_1, nu_i)
      end if

      call TimeDerivative( problem, flow_op, t_i      &
                         , u    = u_1                 &
                         , nu   = nu_i                &
                         , chi  = chi                 &
                         , F_d3 = F_d3 (:,:,:,:,:,1)  &
                         )

      ! stages 2 to n_stage ....................................................

      do i = 2, n_stage

    ! example of debug
    ! if (any(ieee_is_nan(G_im)))then
    ! print *, "there is nan"
    ! end if


       t_i = t_sdc(m-1) + c(i) * dt
       dt_sdc = t_sdc(this % n_sub) - t_sdc(0)
       S_i = dt_sdc * w_rk(0,i,m) * F_sdc(:,:,:,:,:,0)
        do l = 1, this % n_sub
          S_i = S_i + dt_sdc * w_rk(l,i,m) * F_sdc(:,:,:,:,:,l)
        end do

        G_ex(:,:,:,:,:,i) = -l_rk(0,i,m) * F_sdc_ex(:,:,:,:,:,0)
        do l = 1, this % n_sub
          G_ex(:,:,:,:,:,i) = G_ex(:,:,:,:,:,i) - l_rk(l,i,m) * F_sdc_ex(:,:,:,:,:,l)
        end do

        G_im(:,:,:,:,:,i) = -l_rk(0,i,m) * F_sdc_im(:,:,:,:,:,0)
        do l = 1, this % n_sub
          G_im(:,:,:,:,:,i) = G_im(:,:,:,:,:,i) - l_rk(l,i,m) * F_sdc_im(:,:,:,:,:,l)
        end do

        call CorrectorStage( this, i, t_i, dt, nu_i, u_1, u_i, S_i &
                           , G_im, G_ex, F_d3                      )

      end do

      ! assembly ...............................................................

      t_m = t_sdc(m)

      ! u = u_s + S(t_m) - S(t_s)
      if (.not. fsal) then
        u_i = u_i - S_i + w_sdc(m,0) * F_sdc(:,:,:,:,:,0)
        do l = 1, this % n_sub
          u_i = u_i + w_sdc(m,l) * F_sdc(:,:,:,:,:,l)
        end do
      end if

      ! corrections required if b ≠ a(s,:)
      if (any(b_im /= a_im(n_stage,:)) .or. any(b_ex /= a_ex(n_stage,:))) then

         do i = 1, n_stage

           ! adjust implicit contributions
           tau = (b_im(i) -  a_im(n_stage,i)) * dt
           if (abs(tau) > epsilon(ONE)) then
             call MergeArrays(ONE, u_i, tau, G_im (:,:,:,:,:,i), multi=.true.)
           end if

           ! adjust explicit contributions
           tau = (b_ex(i) - a_ex(n_stage,i)) * dt
           if (abs(tau) > epsilon(ONE)) then
             call MergeArrays(ONE, u_i, tau, G_ex (:,:,:,:,:,i), multi=.true.)
           end if

        end do
      end if

      ! enforce continuity, if required
      if (this % project > 0 .or. .not. fsal) then

        if (.not. fsal) then
          ! update boundary conditions
          call GetBoundaryValues(problem, mesh, flow_op%bv_x, t_m, flow_op%bv_u)
        end if

        !  ∇²ψ= ∇⋅ṽ/∆t
        call SetArray(psi, ZERO)
        call PressureSolver(problem, flow_op, dt, u_i, psi, w=S_i)

        ! v = ṽ - 1/∆t ∇ψ - J(v)
        call ProjectionStep(problem, flow_op, dt, psi, u_i, w=S_i)

      end if

      ! RHS and diffusivity updates ............................................

      if (associated(nu_i) .and. .not. fsal) then
        call problem % GetDiffusivity(flow_op%x, t_m, u_i, nu_i)
      end if

      ! RHS contributions, using F_d3 as workspace

      call TimeDerivative( problem, flow_op, t_m, u = u_i, nu = nu_i, chi = chi &
                         , F_c  = F_sdc_ex_new(:,:,:,:,:,m)                     &
                         , F_d1 = F_sdc_im_new(:,:,:,:,:,m)                     &
                         , F_d2 = F_d3(:,:,:,:,:,1)                             &
                         , F_d3 = F_d3(:,:,:,:,:,2)                             )

      do i = 1, problem % nc
        if (i == 4) cycle ! here: F_d2 = F_d3(*,1)
        call MergeArrays(ONE, F_sdc_ex_new(:,:,:,:,i,m), ONE, F_d3(:,:,:,:,i,1))
        if (i >  4) cycle ! here: F_d3 = F_d3(*,2)
        call MergeArrays(ONE, F_sdc_ex_new(:,:,:,:,i,m), ONE, F_d3(:,:,:,:,i,2))
      end do

    end associate

    ! clean up .................................................................

    !$omp barrier
    !$omp single
    if(allocated(S_i ))   deallocate(S_i )
    if(allocated(G_im))   deallocate(G_im)
    if(allocated(G_ex))   deallocate(G_ex)
    if(allocated(F_d3))   deallocate(F_d3)
    !$omp end single

  end subroutine CorrectorStep

  subroutine CorrectorStage(this, i, t_i, dt, nu_i, u_1, u_i, S_i, G_im, G_ex, F_d3)

    ! arguments ................................................................

    class(Flow_SDC_Method_RungeKutta), intent(inout) :: this
    integer,   intent(in)    :: i                !< stage i
    real(RNP), intent(in)    :: t_i              !< stage time tᵢ
    real(RNP), intent(in)    :: dt               !< time step width
    real(RNP), intent(inout) :: nu_i(:,:,:,:,:)  !< variable ν(x,tᵢ)
    real(RNP), intent(in)    :: u_1 (:,:,:,:,:)  !< u(x,t₁)
    real(RNP), intent(out)   :: u_i (:,:,:,:,:)  !< u(x,tᵢ)
    real(RNP), intent(in)    :: S_i (:,:,:,:,:)
    real(RNP), intent(inout) :: F_d3(:,:,:,:,:,:)  !< F_d3 = χ∇(ν∇⋅v)
    real(RNP), intent(inout) :: G_im(:,:,:,:,:,:)
    real(RNP), intent(inout) :: G_ex(:,:,:,:,:,:)

    optional :: nu_i

    ! local variables...........................................................

    real(RNP), dimension(:,:,:,:,:), allocatable, save :: F_c, F_d1, F_d2, q, w

    real(RNP) :: tau
    integer :: j, k

    ! workspace ................................................................

    !$omp single
    allocate(F_c , mold = u_1)
    allocate(F_d1, mold = u_1)
    allocate(F_d2, mold = u_1)
    allocate(q   , mold = u_1)
    allocate(w   , mold = u_1)
    !$omp end single


    ! ..........................................................................

    associate( problem => this % problem           &
             , flow_op => this % flow_op           &
             , mesh    => this % flow_op % mesh    &
             , psi     => u_i(:,:,:,:,4)           &
             , a_ex    => this % imex_rk % a_ex    &
             , a_im    => this % imex_rk % a_im    &
             , eop     => this % flow_op % eop_u   &
             )

      ! initialization .........................................................

      call GetBoundaryValues(problem, mesh, flow_op % bv_x, t_i, flow_op % bv_u)

      ! extrapolation ..........................................................

      ! u' = u₁ + S_i + ∆t ∑ⁱ⁻¹ aᵉˣ(G_ex + G_im + F_d3)
      call SetArray(u_i, u_1, multi=.true.)

      do k = 1, problem % nc
        if (k == 4) cycle ! skip pressure
        call MergeArrays(ONE, u_i(:,:,:,:,k), ONE, S_i(:,:,:,:,k))
        do j = 1, i-1
          tau = dt * a_ex(i,j)
          if (tau /= 0) then
            call MergeArrays(ONE, u_i(:,:,:,:,k), tau, G_ex(:,:,:,:,k,j))
            call MergeArrays(ONE, u_i(:,:,:,:,k), tau, G_im(:,:,:,:,k,j))
            call MergeArrays(ONE, u_i(:,:,:,:,k), tau, F_d3(:,:,:,:,k,j))
          end if
        end do
      end do

      ! projection: u_i ← u", psi ← ψ" .........................................

      call PressureSolver(problem, flow_op, dt, u_i, psi, w) ! ∇²ψ" = ∇⋅v'
      call ProjectionStep(problem, flow_op, dt, psi, u_i, w) ! v"+J(v") = v'-∇ψ"

      ! diffusion: u_i ← u''' ..................................................

      if (present(nu_i) .or. any(problem % nu_ref > 0)) then

        ! update diffusivity
        if (present(nu_i)) then
          call problem % GetDiffusivity(flow_op%x, t_i, u_i, nu_i) ! νᵢ = ν(x,tᵢ,u")
        end if

        ! RHS: q = uᵢ + ∆t[ ∑ⁱ⁻¹(aⁱᵐ-aᵉˣ)G_im - ∑ⁱ⁻¹ aᵉˣF_d3 ]
        do k = 1, problem % nc
          if (k == 4) cycle
          call SetArray(q(:,:,:,:,k), u_i(:,:,:,:,k))
          do j = 1, i
            tau = dt * (a_im(i,j) - a_ex(i,j))
            if (tau /= 0) then
              call MergeArrays(ONE, q(:,:,:,:,k), tau, G_im(:,:,:,:,k,j))
            end if
            tau = -dt * a_ex(i,j)
            if (tau /= 0 .and. k < 4) then
              call MergeArrays(ONE, q(:,:,:,:,k), tau, F_d3(:,:,:,:,k,j))
            end if
          end do
        end do

        ! uᵢ/τ - ∇⋅(ν∇uᵢ) = f/τ,  τ = ∆t aⁱᵐ(i,i)
        tau = dt * a_im(i,i)
        call DiffusionStep(problem, flow_op, tau, q, u_i, w, nu_i)

       end if

      ! update G and F_d3 ....................................................

      ! NOTE: recomputation of nu -- probably not required
      ! if (present(nu_i)) then
      !  if (problem % HasVariableProperties()) then
      !    call problem % GetDiffusivity(flow_op%x, t, u_i, nu_i) ! νᵢ = ν(x,tᵢ,u)
      !  end if
      ! end if

      call TimeDerivative( problem, flow_op, t_i     &
                         , u    = u_i                &
                         , nu   = nu_i               &
                         , chi  = chi                &
                         , F_c  = F_c  (:,:,:,:,:)   &
                         , F_d1 = F_d1 (:,:,:,:,:)   &
                         , F_d2 = F_d2 (:,:,:,:,:)   &
                         , F_d3 = F_d3 (:,:,:,:,:,i) &
                         )

       ! G_im(*,i) += F_d1
       call MergeArrays(ONE, G_im(:,:,:,:,:,i), ONE, F_d1, multi=.true.)

       ! G_ex(*,i) += F_c + F_d2 + F_d3
       call MergeArrays( ONE, G_ex(:,:,:,:,:,i), ONE, F_c , multi=.true.)
       call MergeArrays( ONE, G_im(:,:,:,:,:,i), ONE, F_d2, multi=.true.)
       call MergeArrays( ONE, G_im(:,:,:,:,:,i), ONE, F_d3(:,:,:,:,:,i) &
                       , multi=.true. )

    end associate

    ! clean up .................................................................

    !$omp barrier
    !$omp single
    if (allocated(F_c )) deallocate(F_c  )
    if (allocated(F_d1)) deallocate(F_d1 )
    if (allocated(F_d2)) deallocate(F_d2 )
    if (allocated(w   )) deallocate(w    )
    if (allocated(q   )) deallocate(q    )
    !$omp end single

  end subroutine CorrectorStage

  !=============================================================================

end module CART__ISP_Flow__SDC_Method__Runge_Kutta
