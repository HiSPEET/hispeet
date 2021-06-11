!> summary:  SDC method based on IMEX Runge-Kutta with dual splitting
!> author:   Joerg Stiller
!> date:     2021/04/10
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module CART__ISP_Flow__SDC_Method__Runge_Kutta
  use, intrinsic :: ISO_Fortran_Env, only: OUTPUT_UNIT
  use Kind_Parameters, only: RNP
  use Constants,       only: ONE, ZERO
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

  implicit none
  private

  public :: Flow_SDC_Method_RungeKutta
  public :: Flow_SDC_Options_RungeKutta

  !-----------------------------------------------------------------------------
  !> IMEX Runge-Kutta SDC corrector for incompressible flow

  type, extends(Flow_SDC_Method) :: Flow_SDC_Method_RungeKutta

    type(IMEX_RK_Method) :: imex_rk !< IMEX Runge-Kutta method

    integer :: diffusion !< method for computing the diffusive RHS (F_d1)
    integer :: assembly  !< method for constructing the final solution

    real(RNP), allocatable :: w_rk(:,:,:) !< quadrature weights at RK nodes
    real(RNP), allocatable :: l_rk(:,:,:) !< interpolation weights at RK nodes

  contains

    procedure :: Init => Init_Flow_SDC_Method_RungeKutta
    procedure :: Show => Show_Flow_SDC_Method_RungeKutta
    procedure, private :: CorrectorRHS
    procedure, private :: CorrectorStep

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
    integer :: diffusion = 3 !< method for computing the diffusive RHS (F_d1)
    integer :: assembly  = 1 !< method for constructing the final solution
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
    this % diffusion = sdc_opt % diffusion
    this % assembly  = sdc_opt % assembly

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

    ! specific settings
    write(io,'(2X,A,T15,I0)') 'diffusion:', this % diffusion
    write(io,'(2X,A,T15,I0)') 'assembly:' , this % assembly

    ! show IMEX RK settings
    call this % imex_rk % Show(unit)

    write(io,*)

  end subroutine Show_Flow_SDC_Method_RungeKutta

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
  !> Note: See parent type for argument description.

!TODO! RK-SDC corrector notes:
!TODO!   *  subinterval integral S
!TODO!      •  computed on demand in RK stage i > 1
!TODO!      •  realized in identical way as in Euler corrector
!TODO!      •  FSAL property → S for assembly is identical to S in last stage
!TODO!   *  G = F^{k+1} - F^k takes over the role of F
!TODO!   *  for stage i
!TODO!      •  G_rk_im(:,:,:,:,:,i) ≈ F_im^{k+1} - F_im^{k}
!TODO!      •  G_rk_ex(:,:,:,:,:,i) ≈ F_ex^{k+1} - F_ex^{k}
!TODO!   *  dimensions 1 to 5 of G_rk_im/ex identical to those of u_sdc
!TODO!   *  RK stage 1 coincides with SDC node m-1:
!TODO!      •  G_rk_im(*,1) = F_sdc_im_new(*,m-1) - F_sdc_im(*,m-1)
!TODO!      •  G_rk_ex(*,1) = F_sdc_ex_new(*,m-1) - F_sdc_ex(*,m-1)
!TODO!   *  RK stages i > 1
!TODO!      •  initialization of G(*,i) with interpolated F_sdc_
!TODO!           ‣  G_rk_im(*,i) = -F_rk_im ← interpolated from F_sdc_im
!TODO!           ‣  G_rk_ex(*,i) = -F_rk_ex ← interpolated from F_sdc_ex
!TODO!           ‣  no interpolation required in last stage if FSAL
!TODO!      •  proceed within the stage as usual
!TODO!      •  at end of the stage compute
!TODO!           ‣  F_rk_im = F_d1                    [+ F_p ???]
!TODO!           ‣  F_rk_ex = F_c + F_d2 + F_d3
!TODO!      •  complete G(*,i)
!TODO!           ‣  G_rk_im(*,i) = G_rk_im(*,i) + F_rk_im
!TODO!           ‣  G_rk_ex(*,i) = G_rk_ex(*,i) + F_rk_ex
!TODO!      •  if b == a(n_stage,:), i.e. no assembly
!TODO!           ‣  F_sdc_im_new(*,m) = F_rk_im
!TODO!           ‣  F_sdc_ex_new(*,m) = F_rk_ex
!TODO!   *  if assembly is required
!TODO!      •  compute solution
!TODO!      •  recompute F_rk_im/ex
!TODO!   *  return
!TODO!      •  u_sdc        (*,m) = u
!TODO!      •  nu_sdc       (*,m) = nu           (if required)
!TODO!      •  F_sdc_im_new (*,m) = F_rk_im
!TODO!      •  F_sdc_ex_new (*,m) = F_rk_ex

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

    real(RNP), dimension(:,:,:,:,:)  , allocatable, save :: u_i
    real(RNP), dimension(:,:,:,:,:)  , allocatable, save :: S_i

!TODO! which of following are really needed ?
    real(RNP), dimension(:,:,:,:,:,:), allocatable, save :: F_c
    real(RNP), dimension(:,:,:,:,:,:), allocatable, save :: F_d1, F_d2, F_d3
    real(RNP), dimension(:,:,:,:,:,:), allocatable, save :: F_p, F_s
    real(RNP), dimension(:,:,:,:,:)  , allocatable, save :: nu

    !new
    real(RNP), dimension(:,:,:,:,:,:), allocatable, save :: G_ex
    real(RNP), dimension(:,:,:,:,:,:), allocatable, save :: G_im
    !end new

    real(RNP), save :: t_0 = -huge(ONE)
    logical,   save :: fsal

    ! auxiliary
    real(RNP) :: t_i, tau, dt
    integer   :: i, l

    associate( problem => this % problem           &
             , flow_op => this % flow_op           &
             , mesh    => this % flow_op % mesh    &
             , x       => this % flow_op % x       &
             , u_0     => u_sdc(:,:,:,:,:,m-1)     &
             , p       => u_sdc(:,:,:,:,4,m)       &
             , a_im    => this % imex_rk % a_im    &
             , a_ex    => this % imex_rk % a_ex    &
             , c       => this % imex_rk % c       &
             , b_im    => this % imex_rk % b_im    &
             , b_ex    => this % imex_rk % b_ex    &
             , n_stage => this % imex_rk % n_stage &
             , w_rk    => this % w_rk              &
             , l_rk    => this % l_rk              &
             , eop     => this % flow_op % eop_u   )

      ! initialization .........................................................

      !$omp single

      ! check for reusable workspace
!TODO! needs to updated
      if (allocated(F_c)) then
        if (any(shape(F_c) /= [shape(u_0),n_stage])) then
          deallocate(u_i, S_i, F_c, F_d1, F_d2, F_d3, F_p, F_s)  ! G_ex, G_im ?
          if(allocated(nu))  deallocate(nu)
        end if
      end if

      ! check first-same-as-last condition
      fsal = c(n_stage) == ONE .and. t == t_0 .and. allocated(F_c)

      if (.not. allocated(F_c)) then
        ! allocate new workspace
        allocate(u_i, mold = u_0)
        allocate(S_i, mold = u_0)
        allocate(F_c(size(u_0,1), size(u_0,2), size(u_0,3), size(u_0,4), size(u_0,5), n_stage))
        allocate(F_d1, F_d2, F_d3, F_p, F_s, mold = F_c)
        !new
        allocate(G_ex, mold = F_c)
        allocate(G_im, mold = F_c)
        !endnew
        if (problem % HasVariableProperties()) then
          allocate(nu, mold = u_0)
        end if
      end if

      ! time
      t_0 = t_sdc(m-1)
      dt  = t_sdc(m) - t_sdc(m-1)
!?!   t   = t + dt  ! needed ?

      !$omp end single

      ! stage 1 ................................................................

      t_i = t_0

      if (fsal) then

        ! reuse last stage of previous step
!TODO! which of following are really needed ?
        call SetArray(F_c (:,:,:,:,:,1), F_c (:,:,:,:,:,n_stage), multi = .true.)
        call SetArray(F_d1(:,:,:,:,:,1), F_d1(:,:,:,:,:,n_stage), multi = .true.)
        call SetArray(F_d2(:,:,:,:,:,1), F_d2(:,:,:,:,:,n_stage), multi = .true.)
        call SetArray(F_d3(:,:,:,:,:,1), F_d3(:,:,:,:,:,n_stage), multi = .true.)
        call SetArray(F_p (:,:,:,:,:,1), F_p (:,:,:,:,:,n_stage), multi = .true.)
        call SetArray(F_s (:,:,:,:,:,1), F_s (:,:,:,:,:,n_stage), multi = .true.)

!TODO! question: can G_ex and G_im also be reused ?

      else

        if (problem % HasVariableProperties()) then
          call problem % GetDiffusivity(flow_op%x, t_i, u, nu)
        end if

        call TimeDerivative( problem, flow_op, t_i      &
                           , u    = u                   &
                           , p    = p                   &
                           , nu   = nu                  &
                           , chi  = chi                 &
                           , F_c  = F_c  (:,:,:,:,:,1)  &
                           , F_d1 = F_d1 (:,:,:,:,:,1)  &
                           , F_d2 = F_d2 (:,:,:,:,:,1)  &
                           , F_d3 = F_d3 (:,:,:,:,:,1)  &
                           , F_p  = F_p  (:,:,:,:,:,1)  &
                           , F_s  = F_s  (:,:,:,:,:,1)  &
                           )
      end if

!TODO! interpolation of  F_im/ex^k(*,i)              using this%l_rk(:,:,i)

       !new
       G_im(:,:,:,:,:,1) = F_sdc_im_new(:,:,:,:,:,m-1) - F_sdc_im(:,:,:,:,:,m-1)
       G_ex(:,:,:,:,:,1) = F_sdc_ex_new(:,:,:,:,:,m-1) - F_sdc_ex(:,:,:,:,:,m-1)
       !endnew

      ! stages 2 to n_stage ....................................................

      do i = 2, n_stage

        t_i = t_0 + c(i) * dt

!TODO! add computation   of  S_i = \int_{t_0}^{t_i} F^k  using this%w_rk(:,:,i)
!TODO! if t_i = t → S_i = S
        S_i = w_rk(0,i,m) * F_sdc(:,:,:,:,:,0)
        do l = 1, this % n_sub
          S_i = S_i + w_rk(l,i,m) * F_sdc(:,:,:,:,:,l)
        end o
!new
!nicht fertig!
        G_im = -l_rk(0,i,m) * F_sdc_im(:,:,:,:,:,0)
        do l = 1, this % n_sub
          G_im = G_im - l_rk(l,i,m) * F_sdc_im(:,:,:,:,:,l)
        end do

!TODO!  G_ex = ...
!endnew

!TODO! needs to be further adapted ?
        call CorrectorStage(this, i, t_i, dt, nu, u_0, u_i, S_i, G_im, G_ex, F_d3)

      end do

      ! assembly ...............................................................

!TODO! *** LATER *** adapt as required

      if ( all(b_im ==  this % imex_rk % a_im(n_stage,:)) .and.   &
           all(b_ex ==  this % imex_rk % a_ex(n_stage,:))       ) &
      then

        call SetArray(u, u_i, multi=.true.)

      else

        select case(this % assembly)

        case(1) ! starting from u(x,t₀)

!TODO! add corrector contributions including F_im/ex^k and S^k
          do i = 1, n_stage
            tau = b_ex(i) * dt
            call MergeArrays(ONE, u, tau, F_c (:,:,:,:,:,i), multi=.true.)
            call MergeArrays(ONE, u, tau, F_d2(:,:,:,:,:,i), multi=.true.)
            call MergeArrays(ONE, u, tau, F_d3(:,:,:,:,:,i), multi=.true.)
            tau = b_im(i) * dt
            call MergeArrays(ONE, u, tau, F_d1(:,:,:,:,:,i), multi=.true.)
            call MergeArrays(ONE, u, tau, F_p (:,:,:,:,:,i), multi=.true.)
            call MergeArrays(ONE, u, tau, F_s (:,:,:,:,:,i), multi=.true.)
          end do

        case default ! starting from final stage

          call SetArray(u, u_i, multi=.true.)

          do i = 1, n_stage

!TODO! add corrector contributions including F_im/ex^k and S^k
            ! adjust implicit contributions
            tau = (b_im(i) -  this % imex_rk % a_im(n_stage,i)) * dt
            if (abs(tau) > epsilon(ONE)) then
              call MergeArrays(ONE, u, tau, F_d1(:,:,:,:,:,i), multi=.true.)
              call MergeArrays(ONE, u, tau, F_p (:,:,:,:,:,i), multi=.true.)
              call MergeArrays(ONE, u, tau, F_s (:,:,:,:,:,i), multi=.true.)
            end if

            ! adjust explicit contributions
            tau = (b_ex(i) -  this % imex_rk % a_ex(n_stage,i)) * dt
            if (abs(tau) > epsilon(ONE)) then
              call MergeArrays(ONE, u, tau, F_c (:,:,:,:,:,i), multi=.true.)
              call MergeArrays(ONE, u, tau, F_d2(:,:,:,:,:,i), multi=.true.)
              call MergeArrays(ONE, u, tau, F_d3(:,:,:,:,:,i), multi=.true.)
            end if

          end do

        end select

      end if

!TODO! revise if necessary
      ! enforce continuity:  ∇²δp = ∇⋅ṽ/∆t, v + J(v) = ṽ - ∆t∇δp
      if (this % project > 0) then
        associate(dp => F_p(:,:,:,:,4,1), w => F_s(:,:,:,:,:,1))
          call PressureSolver(problem, flow_op, dt, u, dp, w=w)
          call ProjectionStep(problem, flow_op, dt, dp, u, w=w)
        end associate
      end if

      ! pressure ...............................................................

!TODO! revise if necessary
      select case(this % pressure)
      case(1)
        ! recompute pressure using F_p and F_s as workspace
        associate(F_v => F_p(:,:,:,:,:,1), w => F_s(:,:,:,:,:,1))
          if (c(n_stage) /= ONE .and. problem % HasVariableProperties()) then
            call problem % GetDiffusivity(flow_op%x, t, u, nu)
          end if
          call TimeDerivative(problem, flow_op, t, u, nu=nu, chi=chi, F=F_v)
          call ComputePressure(problem, flow_op, t, F_v, p, w)
        end associate
      case default
        call SetArray(p, u_i(:,:,:,:,4))
      end select

    end associate

    ! save final time for FSAL check ...........................................

    t_0 = t

    ! clean up .................................................................

    ! keep workspace in case of standby
    if (present(standby)) then
      if (standby) return
    end if

    !$omp single
    if(allocated(u_i  ))   deallocate(u_i )
    if(allocated(S_i  ))   deallocate(S_i )
    if(allocated(F_c  ))   deallocate(F_c )
    if(allocated(F_d1 ))   deallocate(F_d1)
    if(allocated(F_d2 ))   deallocate(F_d2)
    if(allocated(F_d3 ))   deallocate(F_d3)
    if(allocated(F_p  ))   deallocate(F_p )
    if(allocated(F_s  ))   deallocate(F_s )
    if(allocated(nu   ))   deallocate(nu  )
    !$omp end single

  end subroutine CorrectorStep

  subroutine CorrectorStage(this, i, t_i, dt, nu, u_0, u_i, S_i, G_im, G_ex, F_d3)
!TODO! 1) complete argument list based on algorithm and on declarations in CorrectorStep
!TODO! 2) expand step by step based on algorithm using RungeKuttaStage_d1 as a model
!TODO! 3) add local declarations as required, e.g:
    real(RNP), allocatable, save :: F_d1 (:,:,:,:,:)  !?  is local
    allocate(F_d1, mold = u_i)                        !?  etc.

  end subroutine CorrectorStage

!   subroutine RungeKuttaStage_d1( this, i, t, dt, nu, u_0, u_i    &
!                                , F_c, F_d1, F_d2, F_d3, F_p, F_s )
!
!     ! arguments ................................................................
!
!     class(TimeIntegrator_RungeKutta), intent(inout) :: this
!     integer,   intent(in)    :: i                 !< stage i
!     real(RNP), intent(in)    :: t                 !< stage time tᵢ
!     real(RNP), intent(in)    :: dt                !< time step width
!     real(RNP), intent(inout) :: nu  (:,:,:,:,:)   !< variable ν(x,t)
!     real(RNP), intent(in)    :: u_0 (:,:,:,:,:)   !< u(x,t₀)
!     real(RNP), intent(out)   :: u_i (:,:,:,:,:)   !< u(x,tᵢ)
!     real(RNP), intent(inout) :: F_c (:,:,:,:,:,:) !< F_c  = -∇⋅(vu)
!     real(RNP), intent(inout) :: F_d1(:,:,:,:,:,:) !< F_d1 =  ∇⋅(ν∇u)
!     real(RNP), intent(inout) :: F_d2(:,:,:,:,:,:) !< F_d2 =  ∇⋅(ν∇v)ᵀ
!     real(RNP), intent(inout) :: F_d3(:,:,:,:,:,:) !< F_d3 = -∇(ν∇⋅v)
!     real(RNP), intent(inout) :: F_p (:,:,:,:,:,:) !< F_p  = -∇p
!     real(RNP), intent(inout) :: F_s (:,:,:,:,:,:) !< F_s  =  f(x,t)
!
!     optional :: nu
!
!     ! local variables...........................................................
!
!     real(RNP), allocatable, save :: f (:,:,:,:,:)
!     real(RNP), allocatable, save :: w (:,:,:,:,:)
!     real(RNP), allocatable, save :: dp(:,:,:,:  )
!
!     real(RNP) :: tau
!     integer   :: j, k
!
!     ! workspace ................................................................
!
!     !$omp single
!     allocate(f , mold = u_0)
!     allocate(w , mold = u_0)
!     allocate(dp, mold = u_0(:,:,:,:,4))
!     !$omp end single
!
!     ! ..........................................................................
!
!     associate( problem => this % problem           &
!              , flow_op => this % flow_op           &
!              , mesh    => this % flow_op % mesh    &
!              , x       => this % flow_op % x       &
!              , p       => u_i(:,:,:,:,4)           &
!              , a_ex    => this % imex_rk % a_ex    &
!              , a_im    => this % imex_rk % a_im    &
!              , eop     => this % flow_op % eop_u   &
!              )
!
!       ! initialization .........................................................
!
!       call GetBoundaryValues(problem, mesh, flow_op % bv_x, t, flow_op % bv_u)
!       call problem % GetExternalSources(flow_op%x, t, F_s(:,:,:,:,:,i))
!       call SetArray(u_i, u_0, multi=.true.)
!
!       ! extrapolation: u_i ← u' ................................................
!TODO! This corresponds to line 3 in CorrectorStage algorithm
!
!       ! u' = u₀ + ∆t ∑ⁱ⁻¹ [aᵉˣ(F_c + F_d) + aⁱᵐ F_p] + ∆t ∑ⁱ aⁱᵐ F_s
!
!       do k = 1, problem % nc
!         if (k == 4) cycle ! skip pressure
!         do j = 1, i-1
!
!           tau = dt * a_ex(i,j)
!           if (tau /= 0) then
!             call MergeArrays(ONE, u_i(:,:,:,:,k), tau, F_c (:,:,:,:,k,j))
!             call MergeArrays(ONE, u_i(:,:,:,:,k), tau, F_d1(:,:,:,:,k,j))
!             if (k < 4) then
!               call MergeArrays(ONE, u_i(:,:,:,:,k),   tau, F_d2(:,:,:,:,k,j))
!               call MergeArrays(ONE, u_i(:,:,:,:,k), 2*tau, F_d3(:,:,:,:,k,j))
!             end if
!           end if
!
!           tau = dt * a_im(i,j)
!           if (tau /= 0) then
!             call MergeArrays(ONE, u_i(:,:,:,:,k), tau, F_p (:,:,:,:,k,j))
!             call MergeArrays(ONE, u_i(:,:,:,:,k), tau, F_s (:,:,:,:,k,j))
!           end if
!
!         end do
!       end do
!
!       tau = dt * a_im(i,i)
!       do k = 1, problem % nc
!         if (k == 4) cycle ! skip pressure
!         call MergeArrays(ONE, u_i(:,:,:,:,k), tau, F_s (:,:,:,:,k,i))
!       end do
!
!       ! projection: u_i ← u", p ← p" ...........................................
!TODO! This corresponds to lines 4-5 in CorrectorStage algorithm
!
!       tau = dt * a_im(i,i)
!       call PressureSolver(problem, flow_op, tau, u_i, p, w) ! ∇²p" = ∇⋅v'/τ
!       call ProjectionStep(problem, flow_op, tau, p, u_i, w) ! v"+J(v") = v'-τ∇p"
!
!       ! diffusion: u_i ← u''' ..................................................
!TODO! This corresponds to line 6 in CorrectorStage algorithm
!
!       if (present(nu) .or. any(problem % nu_ref > 0)) then
!
!         ! update diffusivity
!         if (present(nu)) then
!           call problem % GetDiffusivity(flow_op%x, t, u_i, nu) ! ν = ν(x,t,u")
!         end if
!
!TODO! This corresponds to line 7 in CorrectorStage algorithm
!         ! RHS: f = uᵢ + τ∆F = uᵢ + ∆t[ ∑ⁱ⁻¹(aⁱᵐ-aᵉˣ)F_d1 - ∑ⁱ⁻¹ aᵉˣF_d3 ]
!         do k = 1, problem % nc
!           if (k == 4) cycle
!           call SetArray(f(:,:,:,:,k), u_i(:,:,:,:,k))
!           do j = 1, i-1
!             tau = dt * (a_im(i,j) - a_ex(i,j))
!             if (tau /= 0) then
!               call MergeArrays(ONE, f(:,:,:,:,k), tau, F_d1(:,:,:,:,k,j))
!             end if
!             tau = -dt * a_ex(i,j)
!             if (tau /= 0 .and. k < 4) then
!               call MergeArrays(ONE, f(:,:,:,:,k), tau, F_d3(:,:,:,:,k,j))
!             end if
!           end do
!         end do
!
!TODO! This corresponds to line 8 in CorrectorStage algorithm
!         ! uᵢ/τ - ∇⋅(ν∇uᵢ) = f/τ
!         tau = dt * a_im(i,i)
!         call DiffusionStep(problem, flow_op, tau, f, u_i, w, nu)
!
!       end if
!
!       ! final projection: u_i ← u ..............................................
!
!       if (this % project > 1) then
!
!         tau = dt * a_im(i,i)
!
!         ! ∇²δp = ∇⋅v'''/τ
!         call SetArray(dp, ZERO)
!         call PressureSolver(problem, flow_op, tau, u_i, dp, w)
!         call MergeArrays(ONE, p, ONE, dp)
!
!         ! v + J(v) = v''' - 1/τ ∇p
!         call ProjectionStep(problem, flow_op, tau, dp, u_i(:,:,:,:,1:3), w)
!
!       end if
!
!       ! RHS contributions ......................................................
!TODO! This corresponds to lines 9-11 in CorrectorStage algorithm
!
!       if (present(nu)) then
!         if (problem % HasVariableProperties()) then
!           call problem % GetDiffusivity(flow_op%x, t, u_i, nu) ! ν = ν(x,t,u)
!         end if
!       end if
!
!       ! contributions to time derivative
!       call TimeDerivative( problem, flow_op, t        &
!                          , u    = u_i                 &
!                          , p    = p                   &
!                          , nu   = nu                  &
!                          , chi  = chi                 &
!                          , F_c  = F_c  (:,:,:,:,:,i)  &
!                          , F_d1 = F_d1 (:,:,:,:,:,i)  &
!                          , F_d2 = F_d2 (:,:,:,:,:,i)  &
!                          , F_d3 = F_d3 (:,:,:,:,:,i)  &
!                          , F_p  = F_p  (:,:,:,:,:,i)  &
!                          )
!
!     end associate
!
!     ! clean up .................................................................
!
!     !$omp barrier
!     !$omp single
!     if (allocated(f )) deallocate(f )
!     if (allocated(w )) deallocate(w )
!     if (allocated(dp)) deallocate(dp)
!     !$omp end single
!
!   end subroutine RungeKuttaStage_d1


  !=============================================================================

end module CART__ISP_Flow__SDC_Method__Runge_Kutta
