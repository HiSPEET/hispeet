!> summary:  Runge-Kutta method for incompressible flows with dual splitting
!> author:   Joerg Stiller, Montadhar Guesmi
!> date:     2020/03/05
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> @note
!>    *  This implementation assumes c(1) = 0
!>    *  The pressure is returned as computed in last stage and, hence,
!>       of order 1 or 2 only
!>    *  Schemes with c(s) ≠ 1 may require a recomputation of p, which
!>       is not included in the present scheme
!>    *  Schemes with c(s) = 1 allow for reusing the RHS, which is not
!>       exploited, so far
!> @endnote
!>
!===============================================================================

module CART__ISP_Flow__Time_Integrator__Runge_Kutta
  use, intrinsic :: ISO_Fortran_Env, only: OUTPUT_UNIT
  use Kind_Parameters,   only: RNP
  use Constants,         only: ZERO, ONE
  use Execution_Control, only: Error

  use XMPI
  use Array_Assignments
  use IMEX_Runge_Kutta_Method

  use ISP_Flow_Problem

  use CART__ISP_Flow__Boundary_Values
  use CART__ISP_Flow__Diffusion
  use CART__ISP_Flow__Operators
  use CART__ISP_Flow__Pressure
  use CART__ISP_Flow__Projection
  use CART__ISP_Flow__SGS_Diffusivity
  use CART__ISP_Flow__Time_Derivative
  use CART__ISP_Flow__Time_Integrator
  use CART__DG_Weak_Gradient

  implicit none
  private

  public :: TimeIntegrator_RungeKutta
  public :: TimeIntegrator_RungeKutta_Options

  !-----------------------------------------------------------------------------
  !> IMEX Runge-Kutta method for incompressible flows with dual splitting

  type, extends(TimeIntegrator) :: TimeIntegrator_RungeKutta
    type(IMEX_RK_Method) :: imex_rk  !< IMEX Runge-Kutta method
    integer :: diffusion   !< method for computing diffusive RHS (F_d1)
    integer :: assembly    !< method for constructing the final solution
  contains
    procedure :: Init_TimeIntegrator_RungeKutta
    procedure :: Show => Show_TimeIntegrator_RungeKutta
    procedure :: TimeStep
  end type TimeIntegrator_RungeKutta

  ! overloading the constructor
  interface TimeIntegrator_RungeKutta
    module procedure New_TimeIntegrator_RungeKutta
  end interface

  !-----------------------------------------------------------------------------
  !> Type for providing Runge-Kutta time-integrator options
  !>
  !> Options inherited from base class
  !>
  !>   * `project`                                                            \n
  !>        0: no additional projection step                                  \n
  !>        1: additional projection at the end of the time step              \n
  !>        2: additional projection at the end of each stage
  !>
  !>   * `pressure`
  !>        0: return pressure as computed                                    \n
  !>        1: recompute pressure at the end of the time step
  !>
  !> Further options
  !>
  !>   * `diffusion`                                                          \n
  !>        1: compute diffusive RHS from stage solution: F_d1(ν,u)           \n
  !>        2: compute diffusive RHS as F_d1 = [u - (u' + τ∆F)]/τ - F_p,      \n
  !>           where                                                          \n
  !>             u' is the extrapolated solution,                             \n
  !>             u  is the stage solution,                                    \n
  !>             τ  = ∆t aⁱᵐ(i,i), and                                        \n
  !>             ∆F = ∆t/τ ∑ⁱ⁻¹ [(aⁱᵐ - aᵉˣ)F_d1 - aᵉˣ F_d3)]                 \n
  !>           with F_d3 = -∇ν(∇⋅v)
  !>        3: employs
  !>             diffusion 1 for all, but the last stage, and
  !>             diffusion 2 for the latter
  !>
  !>   * `assembly`                                                           \n
  !>        1: starting from u(x,t₀)                                          \n
  !>        2: starting from last stage
  !>
  !> Options for initializing the IMEX_RK_Method
  !>
  !>   * `n_stage` -- defines the number of stages of the Runge-Kutta method
  !>
  !>   * `method`  -- selector among RK methods with `n_stage` stages

  type, extends(TimeIntegratorOptions) :: TimeIntegrator_RungeKutta_Options
    integer :: n_stage   = 3 !< number of stages
    integer :: method    = 1 !< RK method selector, if more than one exist
    integer :: diffusion = 3 !< method for computing the diffusive RHS (F_d1)
    integer :: assembly  = 1 !< method for constructing the final solution
  contains
    procedure :: Bcast => Bcast_TimeIntegrator_RungeKutta_Options
  end type TimeIntegrator_RungeKutta_Options

  !-----------------------------------------------------------------------------
  ! Constants

  real(RNP), parameter :: chi = -ONE  !< scaling parameter of F_d3 = χ∇(ν∇⋅v)

contains

  !=============================================================================
  ! TimeIntegrator_RungeKutta: type-bound procedures

  !-----------------------------------------------------------------------------
  !> Constructor for objects of type TimeIntegrator_RungeKutta

  function New_TimeIntegrator_RungeKutta(problem, flow_op, opt) result(this)
    class(FlowProblem),                       intent(in) :: problem
    class(FlowOperators),                     intent(in) :: flow_op
    class(TimeIntegrator_RungeKutta_Options), intent(in) :: opt
    type(TimeIntegrator_RungeKutta) :: this

    call Init_TimeIntegrator_RungeKutta(this, problem, flow_op, opt)

  end function New_TimeIntegrator_RungeKutta

  !-----------------------------------------------------------------------------
  !> Initialization of Init_TimeIntegrator_RungeKutta object

  subroutine Init_TimeIntegrator_RungeKutta(this, problem, flow_op, opt)
    class(TimeIntegrator_RungeKutta),         intent(inout) :: this
    class(FlowProblem),                       intent(in)    :: problem
    class(FlowOperators),                     intent(in)    :: flow_op
    class(TimeIntegrator_RungeKutta_Options), intent(in)    :: opt

    ! intialize parent type
    call this % Init_TimeIntegrator(problem, flow_op, opt)
    this % name = 'IMEX Runge-Kutta method'

    ! initialize RK method
    call this % imex_rk % Init_IMEX_RK_Method(opt % n_stage, opt % method)

    ! specific settings
    this % diffusion = opt % diffusion
    this % assembly  = opt % assembly

  end subroutine Init_TimeIntegrator_RungeKutta

  !-----------------------------------------------------------------------------
  !> Output of TimeIntegrator_RungeKutta settings

  subroutine Show_TimeIntegrator_RungeKutta(this, unit)
    class(TimeIntegrator_RungeKutta), intent(in) :: this
    integer,                optional, intent(in) :: unit  !< output unit

    integer :: io

    if (present(unit)) then
      io = unit
    else
      io = OUTPUT_UNIT
    end if

    ! show parent settings
    call this % Show_TimeIntegrator(unit)

    write(io,'(A)')           'TimeIntegrator_RungeKutta settings'
    write(io,'(A,/)')         repeat('-',80)
    write(io,'(2X,A,T15,I0)') 'diffusion:', this % diffusion
    write(io,'(2X,A,T15,I0)') 'assembly:',  this % assembly

    ! show IMEX RK settings
    call this % imex_rk % Show(unit)

  end subroutine Show_TimeIntegrator_RungeKutta

  !-----------------------------------------------------------------------------
  !> Performs a single IMEX Runge-Kutta step

  subroutine TimeStep(this, t, dt, u, standby)
    class(TimeIntegrator_RungeKutta), intent(inout) :: this
    real(RNP),         intent(inout) :: t             !< time t₀ → t
    real(RNP),         intent(in)    :: dt            !< step size ∆t = t-t₀
    real(RNP),         intent(inout) :: u (:,:,:,:,:) !< u(x,t₀) → u(x,t)
    logical, optional, intent(in)    :: standby       !< reuse workspace [F]

    ! local variables  .........................................................

    real(RNP), dimension(:,:,:,:,:)  , allocatable, save :: u_i
    real(RNP), dimension(:,:,:,:,:,:), allocatable, save :: F_c
    real(RNP), dimension(:,:,:,:,:,:), allocatable, save :: F_d1, F_d2, F_d3
    real(RNP), dimension(:,:,:,:,:,:), allocatable, save :: F_p, F_s
    real(RNP), dimension(:,:,:,:,:)  , allocatable, save :: nu

    ! control
    real(RNP), save :: t_0 = -huge(ONE)
    logical,   save :: fsal
    logical,   save :: variable_nu

    ! auxiliary
    real(RNP) :: t_i, tau
    integer   :: i

    associate( problem => this % problem           &
             , flow_op => this % flow_op           &
             , mesh    => this % flow_op % mesh    &
             , x       => this % flow_op % x       &
             , p       => u(:,:,:,:,4)             &
             , a_im    => this % imex_rk % a_im    &
             , a_ex    => this % imex_rk % a_ex    &
             , c       => this % imex_rk % c       &
             , b_im    => this % imex_rk % b_im    &
             , b_ex    => this % imex_rk % b_ex    &
             , n_stage => this % imex_rk % n_stage &
             , eop     => this % flow_op % eop_u   )

      ! initialization .........................................................

      !$omp single

      ! check for reusable workspace
      if (allocated(F_c)) then
        if (any(shape(F_c) /= [shape(u),n_stage])) then
          deallocate(u_i, F_c, F_d1, F_d2, F_d3, F_p, F_s)
          if(allocated(nu))  deallocate(nu)
        end if
      end if

      ! check first-same-as-last condition
      fsal = c(n_stage) == ONE .and. t == t_0 .and. allocated(F_c)
      ! deactivate FSAL if the SGS model is called only at stage 1
      fsal = fsal .and. flow_op % control % sgs /= 1

      ! check for variable diffusivity
      variable_nu = problem % HasVariableProperties() .or. &
                    flow_op % control % sgs > 0

      if (.not. allocated(F_c)) then
        ! allocate new workspace
        allocate(u_i, mold = u)
        allocate(F_c(size(u,1), size(u,2), size(u,3), size(u,4), size(u,5), n_stage))
        allocate(F_d1, F_d2, F_d3, F_p, F_s, mold = F_c)
        if (variable_nu) then
          allocate(nu, mold = u)
        end if
      end if

      ! time
      t_0 = t
      t   = t + dt

      !$omp end single

      ! stage 1 ................................................................

      t_i = t_0

      if (fsal) then

        ! reuse last stage of previous step
        call SetArray(F_c (:,:,:,:,:,1), F_c (:,:,:,:,:,n_stage), multi = .true.)
        call SetArray(F_d1(:,:,:,:,:,1), F_d1(:,:,:,:,:,n_stage), multi = .true.)
        call SetArray(F_d2(:,:,:,:,:,1), F_d2(:,:,:,:,:,n_stage), multi = .true.)
        call SetArray(F_d3(:,:,:,:,:,1), F_d3(:,:,:,:,:,n_stage), multi = .true.)
        call SetArray(F_p (:,:,:,:,:,1), F_p (:,:,:,:,:,n_stage), multi = .true.)
        call SetArray(F_s (:,:,:,:,:,1), F_s (:,:,:,:,:,n_stage), multi = .true.)

      else

        if (variable_nu) then
          call problem % GetDiffusivity(flow_op%x, t_i, u, nu)  ! ν_phy
          if (flow_op % control % sgs > 0) then                 !
            call Add_SGS_Diffusivity(flow_op, u, nu)            ! ν_sgs
          end if
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

      ! stages 2 to n_stage ....................................................

      do i = 2, n_stage
        t_i = t_0 + c(i) * dt
        if (this%diffusion == 1 .or. this%diffusion == 3 .and. i < n_stage) then
          call RungeKuttaStage_d1( this, i, t_i, dt, nu, u, u_i    &
                                 , F_c, F_d1, F_d2, F_d3, F_p, F_s )
        else
          call RungeKuttaStage_d2( this, i, t_i, dt, nu, u, u_i    &
                                 , F_c, F_d1, F_d2, F_d3, F_p, F_s )
        end if
      end do

      ! assembly ...............................................................

      if ( all(b_im ==  this % imex_rk % a_im(n_stage,:)) .and.   &
           all(b_ex ==  this % imex_rk % a_ex(n_stage,:))       ) &
      then

        call SetArray(u, u_i, multi=.true.)

      else

        select case(this % assembly)

        case(1) ! starting from u(x,t₀)

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

      ! enforce continuity:  ∇²δp = ∇⋅ṽ/∆t, v + J(v) = ṽ - ∆t∇δp
      if (this % project > 0) then
        associate(dp => F_p(:,:,:,:,4,1), w => F_s(:,:,:,:,:,1))
          call PressureSolver(problem, flow_op, dt, u, dp, w=w)
          call ProjectionStep(problem, flow_op, dt, dp, u, w=w)
        end associate
      end if

      ! pressure ...............................................................

      select case(this % pressure)
      case(1)
        ! recompute pressure using F_p and F_s as workspace
        associate(F_v => F_p(:,:,:,:,:,1), w => F_s(:,:,:,:,:,1))
          if (c(n_stage) /= ONE .and. variable_nu) then
            call problem % GetDiffusivity(flow_op%x, t, u, nu)
            if (flow_op % control % sgs > 0) then
              call Add_SGS_Diffusivity(flow_op, u, nu)
            end if
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
    if(allocated(F_c  ))   deallocate(F_c )
    if(allocated(F_d1 ))   deallocate(F_d1)
    if(allocated(F_d2 ))   deallocate(F_d2)
    if(allocated(F_d3 ))   deallocate(F_d3)
    if(allocated(F_p  ))   deallocate(F_p )
    if(allocated(F_s  ))   deallocate(F_s )
    if(allocated(nu   ))   deallocate(nu  )
    !$omp end single

  end subroutine TimeStep

  !-----------------------------------------------------------------------------
  !> Executes one IMEX Runge-Kutta stage, RHS computed with TimeDerivative
  !>
  !> @note
  !>   * F_d1(:,:,:,:,:,i) could be used as workspace, e.g. instead of f
  !> @endnote

  subroutine RungeKuttaStage_d1( this, i, t, dt, nu, u_0, u_i    &
                               , F_c, F_d1, F_d2, F_d3, F_p, F_s )

    ! arguments ................................................................

    class(TimeIntegrator_RungeKutta), intent(inout) :: this
    integer,   intent(in)    :: i                 !< stage i
    real(RNP), intent(in)    :: t                 !< stage time tᵢ
    real(RNP), intent(in)    :: dt                !< time step width
    real(RNP), intent(inout) :: nu  (:,:,:,:,:)   !< variable ν(x,t)
    real(RNP), intent(in)    :: u_0 (:,:,:,:,:)   !< u(x,t₀)
    real(RNP), intent(out)   :: u_i (:,:,:,:,:)   !< u(x,tᵢ)
    real(RNP), intent(inout) :: F_c (:,:,:,:,:,:) !< F_c  = -∇⋅(vu)
    real(RNP), intent(inout) :: F_d1(:,:,:,:,:,:) !< F_d1 =  ∇⋅(ν∇u)
    real(RNP), intent(inout) :: F_d2(:,:,:,:,:,:) !< F_d2 =  ∇⋅(ν∇v)ᵀ
    real(RNP), intent(inout) :: F_d3(:,:,:,:,:,:) !< F_d3 = -∇(ν∇⋅v)
    real(RNP), intent(inout) :: F_p (:,:,:,:,:,:) !< F_p  = -∇p
    real(RNP), intent(inout) :: F_s (:,:,:,:,:,:) !< F_s  =  f(x,t)

    optional :: nu

    ! local variables...........................................................

    real(RNP), allocatable, save :: f (:,:,:,:,:)
    real(RNP), allocatable, save :: w (:,:,:,:,:)
    real(RNP), allocatable, save :: dp(:,:,:,:  )

    real(RNP) :: tau
    logical   :: update_nu
    integer   :: j, k

    ! workspace ................................................................

    !$omp single
    allocate(f , mold = u_0)
    allocate(w , mold = u_0)
    allocate(dp, mold = u_0(:,:,:,:,4))
    !$omp end single

    ! ..........................................................................

    associate( problem => this % problem           &
             , flow_op => this % flow_op           &
             , mesh    => this % flow_op % mesh    &
             , x       => this % flow_op % x       &
             , p       => u_i(:,:,:,:,4)           &
             , a_ex    => this % imex_rk % a_ex    &
             , a_im    => this % imex_rk % a_im    &
             , eop     => this % flow_op % eop_u   &
             )

      ! initialization .........................................................

      call GetBoundaryValues(problem, mesh, flow_op % bv_x, t, flow_op % bv_u)
      call problem % GetExternalSources(flow_op%x, t, F_s(:,:,:,:,:,i))
      call SetArray(u_i, u_0, multi=.true.)

      ! check if diffusivity needs to be updated
      update_nu = problem % HasVariableProperties() .or. &
                  flow_op % control % sgs > 1

      ! extrapolation: u_i ← u' ................................................

      ! u' = u₀ + ∆t ∑ⁱ⁻¹ [aᵉˣ(F_c + F_d) + aⁱᵐ F_p] + ∆t ∑ⁱ aⁱᵐ F_s

      do k = 1, problem % nc
        if (k == 4) cycle ! skip pressure
        do j = 1, i-1

          tau = dt * a_ex(i,j)
          if (tau /= 0) then
            call MergeArrays(ONE, u_i(:,:,:,:,k), tau, F_c (:,:,:,:,k,j))
            call MergeArrays(ONE, u_i(:,:,:,:,k), tau, F_d1(:,:,:,:,k,j))
            if (k < 4) then
              call MergeArrays(ONE, u_i(:,:,:,:,k),   tau, F_d2(:,:,:,:,k,j))
              call MergeArrays(ONE, u_i(:,:,:,:,k), 2*tau, F_d3(:,:,:,:,k,j))
            end if
          end if

          tau = dt * a_im(i,j)
          if (tau /= 0) then
            call MergeArrays(ONE, u_i(:,:,:,:,k), tau, F_p (:,:,:,:,k,j))
            call MergeArrays(ONE, u_i(:,:,:,:,k), tau, F_s (:,:,:,:,k,j))
          end if

        end do
      end do

      tau = dt * a_im(i,i)
      do k = 1, problem % nc
        if (k == 4) cycle ! skip pressure
        call MergeArrays(ONE, u_i(:,:,:,:,k), tau, F_s (:,:,:,:,k,i))
      end do

      ! projection: u_i ← u", p ← p" ...........................................

      tau = dt * a_im(i,i)
      call PressureSolver(problem, flow_op, tau, u_i, p, w) ! ∇²p" = ∇⋅v'/τ
      call ProjectionStep(problem, flow_op, tau, p, u_i, w) ! v"+J(v") = v'-τ∇p"

      ! diffusion: u_i ← u''' ..................................................

      if (present(nu) .or. any(problem % nu_ref > 0)) then

         if (update_nu) then
          call problem % GetDiffusivity(flow_op%x, t, u_i, nu)
          if (flow_op % control % sgs > 0) then
            call Add_SGS_Diffusivity(flow_op, u_i, nu)
          end if
        end if

        ! RHS: f = uᵢ + τ∆F = uᵢ + ∆t[ ∑ⁱ⁻¹(aⁱᵐ-aᵉˣ)F_d1 - ∑ⁱ⁻¹ aᵉˣF_d3 ]
        do k = 1, problem % nc
          if (k == 4) cycle
          call SetArray(f(:,:,:,:,k), u_i(:,:,:,:,k))
          do j = 1, i-1
            tau = dt * (a_im(i,j) - a_ex(i,j))
            if (tau /= 0) then
              call MergeArrays(ONE, f(:,:,:,:,k), tau, F_d1(:,:,:,:,k,j))
            end if
            tau = -dt * a_ex(i,j)
            if (tau /= 0 .and. k < 4) then
              call MergeArrays(ONE, f(:,:,:,:,k), tau, F_d3(:,:,:,:,k,j))
            end if
          end do
        end do

        ! uᵢ/τ - ∇⋅(ν∇uᵢ) = f/τ
        tau = dt * a_im(i,i)
        call DiffusionStep(problem, flow_op, tau, f, u_i, w, nu)

      end if

      ! final projection: u_i ← u ..............................................

      if (this % project > 1) then

        tau = dt * a_im(i,i)

        ! ∇²δp = ∇⋅v'''/τ
        call SetArray(dp, ZERO)
        call PressureSolver(problem, flow_op, tau, u_i, dp, w)
        call MergeArrays(ONE, p, ONE, dp)

        ! v + J(v) = v''' - 1/τ ∇p
        call ProjectionStep(problem, flow_op, tau, dp, u_i(:,:,:,:,1:3), w)

      end if

      ! RHS contributions ......................................................

      if (update_nu) then
        call problem % GetDiffusivity(flow_op%x, t, u_i, nu)
        if (flow_op % control % sgs > 0) then
          call Add_SGS_Diffusivity(flow_op, u_i, nu)
        end if
      end if

      ! contributions to time derivative
      call TimeDerivative( problem, flow_op, t        &
                         , u    = u_i                 &
                         , p    = p                   &
                         , nu   = nu                  &
                         , chi  = chi                 &
                         , F_c  = F_c  (:,:,:,:,:,i)  &
                         , F_d1 = F_d1 (:,:,:,:,:,i)  &
                         , F_d2 = F_d2 (:,:,:,:,:,i)  &
                         , F_d3 = F_d3 (:,:,:,:,:,i)  &
                         , F_p  = F_p  (:,:,:,:,:,i)  &
                         )

    end associate

    ! clean up .................................................................

    !$omp barrier
    !$omp single
    if (allocated(f )) deallocate(f )
    if (allocated(w )) deallocate(w )
    if (allocated(dp)) deallocate(dp)
    !$omp end single

  end subroutine RungeKuttaStage_d1

  !-----------------------------------------------------------------------------
  !> Executes one IMEX Runge-Kutta stage

  subroutine RungeKuttaStage_d2( this, i, t, dt, nu, u_0, u_i    &
                               , F_c, F_d1, F_d2, F_d3, F_p, F_s )

    ! arguments ................................................................

    class(TimeIntegrator_RungeKutta), intent(inout) :: this
    integer,   intent(in)    :: i                 !< stage i
    real(RNP), intent(in)    :: t                 !< stage time tᵢ
    real(RNP), intent(in)    :: dt                !< time step width
    real(RNP), intent(inout) :: nu  (:,:,:,:,:)   !< variable ν(x,t)
    real(RNP), intent(in)    :: u_0 (:,:,:,:,:)   !< u(x,t₀)
    real(RNP), intent(out)   :: u_i (:,:,:,:,:)   !< u(x,tᵢ)
    real(RNP), intent(inout) :: F_c (:,:,:,:,:,:) !< F_c  = -∇⋅(vu)
    real(RNP), intent(inout) :: F_d1(:,:,:,:,:,:) !< F_d1 =  ∇⋅(ν∇u)
    real(RNP), intent(inout) :: F_d2(:,:,:,:,:,:) !< F_d2 =  ∇⋅(ν∇v)ᵀ
    real(RNP), intent(inout) :: F_d3(:,:,:,:,:,:) !< F_d3 = -∇(ν∇⋅v)
    real(RNP), intent(inout) :: F_p (:,:,:,:,:,:) !< F_p  = -∇p
    real(RNP), intent(inout) :: F_s (:,:,:,:,:,:) !< F_s  =  f(x,t)

    optional :: nu

    ! local variables...........................................................

    real(RNP), allocatable, save :: f (:,:,:,:,:)
    real(RNP), allocatable, save :: w (:,:,:,:,:)
    real(RNP), allocatable, save :: dp(:,:,:,:  )

    real(RNP) :: tau
    logical   :: update_nu
    integer   :: j, k

    ! workspace ................................................................

    !$omp single
    allocate(f , mold = u_0)
    allocate(w , mold = u_0)
    allocate(dp, mold = u_0(:,:,:,:,4))
    !$omp end single

    ! ..........................................................................

    associate( problem => this % problem           &
             , flow_op => this % flow_op           &
             , mesh    => this % flow_op % mesh    &
             , x       => this % flow_op % x       &
             , p       => u_i(:,:,:,:,4)           &
             , a_ex    => this % imex_rk % a_ex    &
             , a_im    => this % imex_rk % a_im    &
             , eop     => this % flow_op % eop_u   &
             )

      ! initialization .........................................................

      call GetBoundaryValues(problem, mesh, flow_op % bv_x, t, flow_op % bv_u)
      call problem % GetExternalSources(flow_op%x, t, F_s(:,:,:,:,:,i))
      call SetArray(u_i, u_0, multi=.true.)

      ! check if diffusivity needs to be updated
      update_nu = problem % HasVariableProperties() .or. &
                  flow_op % control % sgs > 1

      ! extrapolation: u_i ← u' ................................................

      ! u' = u₀ + ∆t ∑ⁱ⁻¹ [aᵉˣ(F_c + F_d) + aⁱᵐ F_p] + ∆t ∑ⁱ aⁱᵐ F_s

      do k = 1, problem % nc
        if (k == 4) cycle ! skip pressure
        do j = 1, i-1

          tau = dt * a_ex(i,j)
          if (tau /= 0) then
            call MergeArrays(ONE, u_i(:,:,:,:,k), tau, F_c (:,:,:,:,k,j))
            call MergeArrays(ONE, u_i(:,:,:,:,k), tau, F_d1(:,:,:,:,k,j))
            if (k < 4) then
              call MergeArrays(ONE, u_i(:,:,:,:,k),   tau, F_d2(:,:,:,:,k,j))
              call MergeArrays(ONE, u_i(:,:,:,:,k), 2*tau, F_d3(:,:,:,:,k,j))
            end if
          end if

          tau = dt * a_im(i,j)
          if (tau /= 0) then
            call MergeArrays(ONE, u_i(:,:,:,:,k), tau, F_p (:,:,:,:,k,j))
            call MergeArrays(ONE, u_i(:,:,:,:,k), tau, F_s (:,:,:,:,k,j))
          end if

        end do
      end do

      tau = dt * a_im(i,i)
      do k = 1, problem % nc
        if (k == 4) cycle ! skip pressure
        call MergeArrays(ONE, u_i(:,:,:,:,k), tau, F_s (:,:,:,:,k,i))
      end do

      ! initialize diffusive RHS ...............................................

        ! RHS: f = uᵢ + τ∆F = uᵢ + ∆t[ ∑ⁱ⁻¹(aⁱᵐ-aᵉˣ)F_d1 - ∑ⁱ⁻¹ aᵉˣF_d3 ]
        do k = 1, problem % nc
          if (k == 4) cycle
          call SetArray(f(:,:,:,:,k), u_i(:,:,:,:,k))
          do j = 1, i-1
            tau = dt * (a_im(i,j) - a_ex(i,j))
            if (tau /= 0) then
              call MergeArrays(ONE, f(:,:,:,:,k), tau, F_d1(:,:,:,:,k,j))
            end if
            tau = -dt * a_ex(i,j)
            if (tau /= 0 .and. k < 4) then
              call MergeArrays(ONE, f(:,:,:,:,k), tau, F_d3(:,:,:,:,k,j))
            end if
          end do
        end do

      do k = 1, problem % nc
        if (k == 4) cycle
        call SetArray(f(:,:,:,:,k), ZERO)

        do j = 1, i-1

          ! f = ∆t ∑ⁱ⁻¹ (aⁱᵐ - aᵉˣ) F_d1
          tau = dt * (a_im(i,j) - a_ex(i,j))
          if (tau /= 0) then
            call MergeArrays(ONE, f(:,:,:,:,k), tau, F_d1(:,:,:,:,k,j))
          end if

          ! f = f - ∆t ∑ⁱ⁻¹ aᵉˣ F_d3
          tau = -dt * a_ex(i,j)
          if (tau /= 0 .and. k < 4) then
            call MergeArrays(ONE, f(:,:,:,:,k), tau, F_d3(:,:,:,:,k,j))
          end if

        end do

        ! F_d1(*,i) = -(u' + τ∆F)/τ
        tau = dt * a_im(i,i)
        call MergeArrays(ZERO, F_d1(:,:,:,:,k,i), -1/tau, u_i(:,:,:,:,k)) ! -u'/τ
        call MergeArrays(ONE , F_d1(:,:,:,:,k,i), -1/tau, f  (:,:,:,:,k)) ! -∆F

      end do

      ! projection: u_i ← u", p ← p" ...........................................

      tau = dt * a_im(i,i)
      call PressureSolver(problem, flow_op, tau, u_i, p, w) ! ∇²p" = ∇⋅v'/τ
      call ProjectionStep(problem, flow_op, tau, p, u_i, w) ! v"+J(v") = v'-τ∇p"

      ! diffusion: u_i ← u''' ..................................................

      if (present(nu) .or. any(problem % nu_ref > 0)) then

        if (update_nu) then
          call problem % GetDiffusivity(flow_op%x, t, u_i, nu)
          if (flow_op % control % sgs > 0) then
            call Add_SGS_Diffusivity(flow_op, u_i, nu)
          end if
        end if

        ! RHS: f =  u_i + f  = uᵢ + τ∆F
        do k = 1, problem % nc
          if (k == 4) cycle
          call MergeArrays(ONE, f(:,:,:,:,k), ONE, u_i(:,:,:,:,k))
        end do

        ! u'''/τ - ∇⋅(ν∇u''') = f/τ,  τ = ∆t aⁱᵐ(i,i)
        tau = dt * a_im(i,i)
        call DiffusionStep(problem, flow_op, tau, f, u_i, w, nu)

      end if

      ! final projection: u_i ← u ..............................................

      if (this % project > 1) then

        tau = dt * a_im(i,i)

        ! ∇²δp = ∇⋅v'''/τ
        call SetArray(dp, ZERO)
        call PressureSolver(problem, flow_op, tau, u_i, dp, w)
        call MergeArrays(ONE, p, ONE, dp)

        ! v + J(v) = v''' - 1/τ ∇p
        call ProjectionStep(problem, flow_op, tau, dp, u_i(:,:,:,:,1:3), w)

      end if

      ! RHS contributions ......................................................

      ! update ν if required
      if (update_nu) then
        call problem % GetDiffusivity(flow_op%x, t, u_i, nu)
        if (flow_op % control % sgs > 0) then
          call Add_SGS_Diffusivity(flow_op, u_i, nu)
        end if
      end if

      ! contributions to time derivative
      call TimeDerivative( problem, flow_op, t        &
                         , u    = u_i                 &
                         , p    = p                   &
                         , nu   = nu                  &
                         , chi  = chi                 &
                         , F_c  = F_c  (:,:,:,:,:,i)  &
                         , F_d2 = F_d2 (:,:,:,:,:,i)  &
                         , F_d3 = F_d3 (:,:,:,:,:,i)  &
                         , F_p  = F_p  (:,:,:,:,:,i)  &
                         )

      ! F_d1 += u/τ - F_p, with F_d1 as computed before first projection
      tau = dt * a_im(i,i)
      if (tau /= 0) then
        do k = 1, problem % nc
          if (k == 4) cycle
          call MergeArrays(ONE, F_d1(:,:,:,:,k,i), 1/tau, u_i(:,:,:,:,k))
          call MergeArrays(ONE, F_d1(:,:,:,:,k,i), -ONE , F_p(:,:,:,:,k,i))
        end do
      end if

    end associate

    ! clean up .................................................................

    !$omp barrier
    !$omp single
    if (allocated(f )) deallocate(f )
    if (allocated(w )) deallocate(w )
    if (allocated(dp)) deallocate(dp)
    !$omp end single

  end subroutine RungeKuttaStage_d2

  !=============================================================================
  ! TimeIntegrator_RungeKutta_Options: type-bound procedures

  !-----------------------------------------------------------------------------
  !> MPI broadcasting of time-integrator options

  subroutine Bcast_TimeIntegrator_RungeKutta_Options(this, root, comm)
    class(TimeIntegrator_RungeKutta_Options), intent(inout) :: this
    integer,        intent(in) :: root !< rank of broadcast root
    type(MPI_Comm), intent(in) :: comm !< MPI communicator

    type(MPI_Request) :: request(4)
    type(MPI_Status)  :: stat(size(request))
    integer :: n

    ! broadcast options of parent class
    call this % TimeIntegratorOptions % Bcast(root, comm)

    n = 1
    call XMPI_Ibcast( this % n_stage  , root, comm, request(n) );  n = n + 1
    call XMPI_Ibcast( this % method   , root, comm, request(n) );  n = n + 1
    call XMPI_Ibcast( this % diffusion, root, comm, request(n) );  n = n + 1
    call XMPI_Ibcast( this % assembly , root, comm, request(n) )

    call MPI_Waitall(n, request, stat)

  end subroutine Bcast_TimeIntegrator_RungeKutta_Options

  !===========================================================================

end module CART__ISP_Flow__Time_Integrator__Runge_Kutta
