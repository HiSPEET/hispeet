!> summary:  SDC corrector based on IMEX Runge-Kutta with dual splitting
!> author:   Joerg Stiller
!> date:     2021/04/10
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module CART__ISP_Flow__SDC_Corrector__Runge_Kutta
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

  use CART__ISP_Flow__SDC_Corrector

  implicit none
  private

  public :: SDC_Corrector_RungeKutta
  public :: SDC_CorrectorOptions_RungeKutta

  !-----------------------------------------------------------------------------
  !> IMEX Runge-Kutta SDC corrector for incompressible flow

  type, extends(SDC_Corrector) :: SDC_Corrector_RungeKutta

    type(IMEX_RK_Method) :: imex_rk !< IMEX Runge-Kutta method

    integer :: diffusion !< method for computing the diffusive RHS (F_d1)
    integer :: assembly  !< method for constructing the final solution

    real(RNP), allocatable :: w_rk(:,:,:) !< quadrature weights at RK nodes
    real(RNP), allocatable :: l_rk(:,:,:) !< interpolation weights at RK nodes

  contains

    procedure :: Init_SDC_Corrector_RungeKutta
    procedure :: Show => Show_SDC_Corrector_RungeKutta
    procedure :: GetCorrectorRHS
    procedure :: CorrectionStep

  end type SDC_Corrector_RungeKutta

  ! overloading the constructor
  interface SDC_Corrector_RungeKutta
    module procedure New_SDC_Corrector_RungeKutta
  end interface

  !-----------------------------------------------------------------------------
  !> Type for providing Runge-Kutta SDC-corrector options (none, so far)

  type, extends(SDC_CorrectorOptions) :: SDC_CorrectorOptions_RungeKutta
    integer :: n_stage   = 3 !< number of stages
    integer :: method    = 1 !< RK method selector, if more than one exist
    integer :: diffusion = 3 !< method for computing the diffusive RHS (F_d1)
    integer :: assembly  = 1 !< method for constructing the final solution
  end type SDC_CorrectorOptions_RungeKutta

  !-----------------------------------------------------------------------------
  ! Constants

  real(RNP), parameter :: chi = -ONE  !< scaling parameter of F_d3 = χ∇(ν∇⋅v)

contains

  !=============================================================================
  ! SDC_Corrector_RungeKutta: type-bound procedures

  !-----------------------------------------------------------------------------
  !> Constructor for objects of type SDC_Corrector_RungeKutta

  function New_SDC_Corrector_RungeKutta(problem, flow_op, opt) result(this)
    class(FlowProblem),   target, intent(in) :: problem !< flow problem
    class(FlowOperators), target, intent(in) :: flow_op !< flow operators
    class(SDC_CorrectorOptions_RungeKutta), intent(in) :: opt !< SDC options
    type(SDC_Corrector_RungeKutta) :: this

    call Init_SDC_Corrector_RungeKutta(this, problem, flow_op, opt)

  end function New_SDC_Corrector_RungeKutta

  !-----------------------------------------------------------------------------
  !> Initialization of a SDC_Corrector_RungeKutta object

  subroutine Init_SDC_Corrector_RungeKutta(this, problem, flow_op, opt)
    class(SDC_Corrector_RungeKutta), intent(inout) :: this
    class(FlowProblem),      target, intent(in)    :: problem !< flow problem
    class(FlowOperators),    target, intent(in)    :: flow_op !< flow operators
    class(SDC_CorrectorOptions_RungeKutta), intent(in) :: opt !< SDC options

    real(RNP) :: t0, t1, ti
    integer   :: i, j, m

    ! parent type ..............................................................

    call this % Init_SDC_Corrector(problem, flow_op, opt)
    this % name = 'IMEX Runge-Kutta corrector'

    ! Runge-Kutta method .......................................................

    call this % imex_rk % Init_IMEX_RK_Method(opt % n_stage, opt % method)

    this % diffusion = opt % diffusion
    this % assembly  = opt % assembly

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

  end subroutine Init_SDC_Corrector_RungeKutta

  !-----------------------------------------------------------------------------
  !> Output of SDC_Corrector_RungeKutta settings

  subroutine Show_SDC_Corrector_RungeKutta(this, unit)
    class(SDC_Corrector_RungeKutta), intent(in) :: this
    integer, optional, intent(in) :: unit  !< output unit

    integer :: io

    if (present(unit)) then
      io = unit
    else
      io = OUTPUT_UNIT
    end if

    ! show parent settings
    call this % Show_SDC_Corrector(unit)

    ! specific settings
    write(io,'(A)')           'SDC_Corrector_RungeKutta settings'
    write(io,'(A,/)')         repeat('-',80)
    write(io,'(2X,A,T15,I0)') 'diffusion:', this % diffusion
    write(io,'(2X,A,T15,I0)') 'assembly:' , this % assembly

    ! show IMEX RK settings
    call this % imex_rk % Show(unit)

    write(io,*)

  end subroutine Show_SDC_Corrector_RungeKutta

  !-----------------------------------------------------------------------------
  !> Computes F_ex and F_im as defined in the corrector and F for subintegrals

  subroutine GetCorrectorRHS(this, t, nu_c, nu_s, u, F_ex, F_im, F)

    class(SDC_Corrector_RungeKutta), intent(in) :: this
    real(RNP), intent(in)  :: t               !< time
    real(RNP), intent(in)  :: nu_c(:,:,:,:,:) !< variable ν for corrector RHS
    real(RNP), intent(in)  :: nu_s(:,:,:,:,:) !< variable ν for subintegral RHS
    real(RNP), intent(in)  :: u   (:,:,:,:,:) !< u
    real(RNP), intent(out) :: F_ex(:,:,:,:,:) !< explicit RHS for corrector
    real(RNP), intent(out) :: F_im(:,:,:,:,:) !< implicit RHS for corrector
    real(RNP), intent(out) :: F   (:,:,:,:,:) !< RHS for subintegrals

    optional :: nu_c, nu_s

    real(RNP), allocatable, save :: F_d2(:,:,:,:,:)
    real(RNP), allocatable, save :: F_d3(:,:,:,:,:)
    integer :: i

    !$omp single
    allocate(F_d2, mold = u)
    allocate(F_d3, mold = u)
    !$omp end single

    ! RHS for high-order subintegrals ..........................................

    call TimeDerivative( this % problem     &
                       , this % flow_op     &
                       , t                  &
                       , u    = u           &
                       , chi  = chi         &
                       , nu   = nu_s        &
                       , F    = F           &
                       , F_c  = F_ex        )

    ! RHS for low-order part of corrector ......................................

    call TimeDerivative( this % problem     &
                       , this % flow_op     &
                       , t                  &
                       , u    = u           &
                       , chi  = chi         &
                       , nu   = nu_c        &
                       , F_d1 = F_im        &
                       , F_d2 = F_d2        &
                       , F_d3 = F_d3        )

    do i = 1, 3
      call MergeArrays(ONE, F_ex(:,:,:,:,i), ONE, F_d2(:,:,:,:,i))
      call MergeArrays(ONE, F_ex(:,:,:,:,i), ONE, F_d3(:,:,:,:,i))
    end do

    !$omp single
    deallocate(F_d2, F_d3)
    !$omp end single

  end subroutine GetCorrectorRHS

  !-----------------------------------------------------------------------------
  !> Execution of a single correction step

  subroutine CorrectionStep( this, t, dt                   &
                           , F_ex_0_old, F_ex_0, F_ex_old  &
                           , F_im_0_old, F_im_0, F_im_old  &
                           , S, u_0, u, nu                 )

    class(SDC_Corrector_RungeKutta), intent(inout) :: this
    real(RNP), intent(inout) :: t                      !< time t₀ → t
    real(RNP), intent(in)    :: dt                     !< step size ∆t = t-t₀
    real(RNP), intent(in)    :: F_ex_0_old (:,:,:,:,:) !< F^ex (t₀)ᵏ
    real(RNP), intent(in)    :: F_ex_0     (:,:,:,:,:) !< F^ex (t₀)ᵏ⁺¹
    real(RNP), intent(in)    :: F_ex_old   (:,:,:,:,:) !< F^ex (t )ᵏ    ! unused
    real(RNP), intent(in)    :: F_im_0_old (:,:,:,:,:) !< F^im (t₀)ᵏ    ! unused
    real(RNP), intent(in)    :: F_im_0     (:,:,:,:,:) !< F^im (t₀)ᵏ⁺¹  ! unused
    real(RNP), intent(in)    :: F_im_old   (:,:,:,:,:) !< F^im (t )ᵏ
    real(RNP), intent(in)    :: S          (:,:,:,:,:) !< S    (t₀)ᵏ
    real(RNP), intent(in)    :: u_0        (:,:,:,:,:) !< u    (t₀)ᵏ⁺¹
    real(RNP), intent(inout) :: u          (:,:,:,:,:) !< u    (t )ᵏ → (t)ᵏ⁺¹
    real(RNP), intent(inout) :: nu         (:,:,:,:,:) !< v    (t₀)ᵏ⁺¹  ! unused

    optional :: nu

    ! local variables  .........................................................

    real(RNP), dimension(:,:,:,:,:)  , allocatable, save :: u_i
    real(RNP), dimension(:,:,:,:,:,:), allocatable, save :: F_c
    real(RNP), dimension(:,:,:,:,:,:), allocatable, save :: F_d1, F_d2, F_d3
    real(RNP), dimension(:,:,:,:,:,:), allocatable, save :: F_p, F_s
    real(RNP), dimension(:,:,:,:,:)  , allocatable, save :: nu
    real(RNP), dimension(:)          , allocatable, save :: ts

    real(RNP), save :: t_0 = -huge(ONE)
    logical,   save :: fsal

    ! auxiliary
    real(RNP) :: tau, d3
    integer   :: i !?!, k

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
          deallocate(ts, u_i, F_c, F_d1, F_d2, F_d3, F_p, F_s)
          if(allocated(nu))   deallocate(nu  )
        end if
      end if

      ! check first-same-as-last condition
      fsal = c(n_stage) == ONE .and. t == t_0 .and. allocated(F_c) !.and. .false.

      if (.not. allocated(F_c)) then
        ! allocate new workspace
        allocate(ts(1:n_stage))
        allocate(u_i, mold = u)
        allocate(F_c(size(u,1), size(u,2), size(u,3), size(u,4), size(u,5), n_stage))
        allocate(F_d1, F_d2, F_d3, F_p, F_s, mold = F_c)
        if (problem % HasVariableProperties()) then
          allocate(nu, mold = u)
        end if
      end if

      ! node times
      do k = 1, n_stage
        ts(k) = t + c(k) * dt
      end do

      ! time
      t_0 = t
      t   = t + dt

      !$omp end single

      ! stage 1 ................................................................

      if (fsal) then

        ! reuse last stage of previous step
        call SetArray(F_c (:,:,:,:,:,1), F_c (:,:,:,:,:,n_stage), multi = .true.)
        call SetArray(F_d1(:,:,:,:,:,1), F_d1(:,:,:,:,:,n_stage), multi = .true.)
        call SetArray(F_d2(:,:,:,:,:,1), F_d2(:,:,:,:,:,n_stage), multi = .true.)
        call SetArray(F_d3(:,:,:,:,:,1), F_d3(:,:,:,:,:,n_stage), multi = .true.)
        call SetArray(F_p (:,:,:,:,:,1), F_p (:,:,:,:,:,n_stage), multi = .true.)
        call SetArray(F_s (:,:,:,:,:,1), F_s (:,:,:,:,:,n_stage), multi = .true.)

      else

        if (problem % HasVariableProperties()) then
          call problem % GetDiffusivity(flow_op%x, ts(1), u, nu)
        end if

        call TimeDerivative( problem, flow_op, ts(1)    &
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

!EU!      ! boundary conditions ....................................................
!EU!
!EU!      call GetBoundaryValues(problem, mesh, flow_op%bv_x, t, flow_op%bv_u)
!EU!
!EU!      ! explicit + extrapolated diffusive parts ................................
!EU!
!EU!      ! compute F_im = F_d1(v₀ᵏ⁺¹,uᵏ) and F_d3(v₀ᵏ⁺¹,uᵏ)
!EU!      call TimeDerivative( this % problem     &
!EU!                         , this % flow_op     &
!EU!                         , t_0                &
!EU!                         , u    = u           &
!EU!                         , chi  = this % chi  &
!EU!                         , nu   = nu          &
!EU!                         , F_d1 = F_im        &
!EU!                         , F_d3 = F_d3        )
!EU!
!EU!      ! u' = u₀ + Sᵏ + ∆t [ F^ex(t₀)ᵏ⁺¹    - F^ex(t₀)ᵏ
!EU!      !                   + F^im(v₀ᵏ⁺¹,uᵏ) - F^im(v₀ᵏ,uᵏ)
!EU!      !                   + F_d3(v₀ᵏ⁺¹,uᵏ) ]
!EU!      call SetArray(u_i, u_0, multi=.true.)
!EU!      do i = 1, problem % nc
!EU!        if (i == 4) cycle
!EU!        call MergeArrays(ONE, u_i(:,:,:,:,i), ONE, S         (:,:,:,:,i))
!EU!        call MergeArrays(ONE, u_i(:,:,:,:,i),  dt, F_ex_0    (:,:,:,:,i))
!EU!        call MergeArrays(ONE, u_i(:,:,:,:,i), -dt, F_ex_0_old(:,:,:,:,i))
!EU!        call MergeArrays(ONE, u_i(:,:,:,:,i),  dt, F_im      (:,:,:,:,i))
!EU!        call MergeArrays(ONE, u_i(:,:,:,:,i), -dt, F_im_old  (:,:,:,:,i))
!EU!        if (i < 4 .and. this % chi /= 0) then
!EU!          call MergeArrays(ONE, u_i(:,:,:,:,i), dt*this%cd3, F_d3(:,:,:,:,i))
!EU!        end if
!EU!      end do
!EU!
!EU!      ! projection .............................................................
!EU!
!EU!      ! solve for p = p"
!EU!      call PressureSolver(problem, flow_op, dt, u_i, p, w)
!EU!
!EU!      ! v" = v' - 1/∆t ∇p" - J(v")
!EU!      call ProjectionStep(problem, flow_op, dt, p, u_i, w)
!EU!
!EU!      ! diffusion ...............................................................
!EU!
!EU!      ! starting values and RHS for diffusion
!EU!      do i = 1, problem % nc
!EU!        if (i == 4) cycle
!EU!        call SetArray(u(:,:,:,:,i), u_i(:,:,:,:,i))
!EU!        call MergeArrays(ONE, u_i(:,:,:,:,i), -dt, F_im(:,:,:,:,i))
!EU!        if (i < 4 .and. this % cd3 /= 0) then
!EU!          call MergeArrays(ONE, u_i(:,:,:,:,i), -dt*cd3, F_d3(:,:,:,:,i))
!EU!        end if
!EU!      end do
!EU!
!EU!      ! solve u'''/∆t - ν(t₀)ᵏ⁺¹ u''' = u_i/∆t
!EU!      call DiffusionStep(problem, flow_op, dt, f=u_i, u=u, w=w, nu=nu)
!EU!
!EU!      ! final projection .......................................................
!EU!
!EU!      if (this % project > 0) then
!EU!
!EU!        ! p = p" + ∆p  with ∆p such that  ∇²∆p = ∆t ∇⋅(v''' - v) = 0
!EU!        call SetArray(dp, ZERO)
!EU!        call PressureSolver(problem, flow_op, dt, u, dp, w)
!EU!        call MergeArrays(ONE, p, ONE, dp)
!EU!
!EU!        ! v = v''' - 1/∆t ∇p - J(v)
!EU!        call ProjectionStep(problem, flow_op, dt, dp, u, w)
!EU!
!EU!      end if

      ! clean-up ...............................................................

      !$omp barrier
      !$omp single
!EU!      deallocate(u_i )
!EU!      deallocate(F_im)
!EU!      deallocate(F_d3)
!EU!      deallocate(w   )
!EU!      deallocate(dp  )
      !$omp end single

    end associate

    ! touch possibly unused arguments to suppress compiler warnings ;)
!EU!    if (size(F_ex_0_old) > 0) return
!EU!    if (size(F_ex_0    ) > 0) return
!EU!    if (size(F_ex_old  ) > 0) return
!EU!    if (size(F_im_0_old) > 0) return
!EU!    if (size(F_im_0    ) > 0) return

  end subroutine CorrectionStep

  !=============================================================================

end module CART__ISP_Flow__SDC_Corrector__Runge_Kutta
