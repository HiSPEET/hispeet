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
  use Kind_Parameters,   only: RNP
  use Constants,         only: ONE, ZERO
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
    integer :: project  = 0 !< switch for additional projections
    integer :: pressure = 0 !< switch for pressure recomputation
    integer :: variant  = 1 !< RK stage variant
  contains
    procedure :: Init_TimeIntegrator_RungeKutta
    procedure :: TimeStep
  end type TimeIntegrator_RungeKutta

  ! overloading the constructor
  interface TimeIntegrator_RungeKutta
    module procedure New_TimeIntegrator_RungeKutta
  end interface

  !-----------------------------------------------------------------------------
  !> Type for providing Runge-Kutta- time-integrator options

  type, extends(TimeIntegratorOptions) :: TimeIntegrator_RungeKutta_Options
    integer :: n_stage  = 3 !< number of stages
    integer :: method   = 1 !< RK method selector, if more than one exist
    integer :: project  = 0 !< switch for additional projections
    integer :: pressure = 0 !< switch for pressure recomputation
    integer :: variant  = 2 !< RK stage variant
  contains
    procedure :: Bcast => TimeIntegrator_RungeKutta_Options_Bcast
  end type TimeIntegrator_RungeKutta_Options

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
    call this % Init_TimeIntegrator(problem, flow_op)

    ! initialize RK method
    call this % imex_rk % Init_IMEX_RK_Method(opt % n_stage, opt % method)

    ! TK stage variant
    this % variant = opt % variant

  end subroutine Init_TimeIntegrator_RungeKutta

  !-----------------------------------------------------------------------------
  !> Performs a single IMEX Runge-Kutta step

  subroutine TimeStep(this, t, dt, u)
    class(TimeIntegrator_RungeKutta), intent(inout) :: this
    real(RNP), intent(inout) :: t             !< time t₀ → t
    real(RNP), intent(in)    :: dt            !< step size ∆t = t-t₀
    real(RNP), intent(inout) :: u (:,:,:,:,:) !< u(x,t₀) → u(x,t)

    ! local variables  .........................................................

    real(RNP), dimension(:,:,:,:,:)  , allocatable, save :: u_i
    real(RNP), dimension(:,:,:,:,:,:), allocatable, save :: F_c
    real(RNP), dimension(:,:,:,:,:,:), allocatable, save :: F_d1, F_d2, F_d3
    real(RNP), dimension(:,:,:,:,:,:), allocatable, save :: F_p, F_s
    real(RNP), dimension(:,:,:,:,:)  , allocatable, save :: nu
    real(RNP), dimension(:)          , allocatable, save :: ts

    ! auxiliary
    integer   :: i, k
    real(RNP) :: tau

    associate( problem => this % problem           &
             , flow_op => this % flow_op           &
             , mesh    => this % flow_op % mesh    &
             , x       => this % flow_op % x       &
             , p       => u(:,:,:,:,4)             &
             , a_im    => this % imex_rk % a_im    &
             , a_ex    => this % imex_rk % a_ex    &
             , c       => this % imex_rk % c       &
             , b       => this % imex_rk % b       &
             , ns      => this % imex_rk % n_stage &
             , eop     => this % flow_op % eop_u   )
p = 0

      ! workspace ................................................................

      !$omp single
      allocate(ts(1:ns))
      allocate(u_i, mold = u)
      allocate(F_c(size(u,1), size(u,2), size(u,3), size(u,4), size(u,5), ns))
      allocate(F_d1, F_d2, F_d3, F_p, F_s, mold = F_c)
      if (problem % HasVariableProperties()) then
        allocate(nu, mold = u)
      end if
      !$omp end single

      ! node times .............................................................

      do k = 1, ns
        ts(k) = t + c(k) * dt
      end do

      ! stages 1 ...............................................................

      if (problem % HasVariableProperties()) then
        call problem % GetDiffusivity(flow_op%x, ts(1), u, nu)
      end if

      call TimeDerivative( problem, flow_op, ts(1)    &
                         , u_c  = u                   &
                         , u_d  = u                   &
                         , p    = p                   &
                         , nu   = nu                  &
                         , F_c  = F_c  (:,:,:,:,:,1)  &
                         , F_d1 = F_d1 (:,:,:,:,:,1)  &
                         , F_d2 = F_d2 (:,:,:,:,:,1)  &
                         , F_d3 = F_d3 (:,:,:,:,:,1)  &
                         , F_p  = F_p  (:,:,:,:,:,1)  &
                         , F_s  = F_s  (:,:,:,:,:,1)  &
                         )

      ! stages 2 to ns .........................................................

      do i = 2, ns
        select case(this % variant)
        case(2)
          call RungeKuttaStage_v2( this, i, ts(i), dt, nu, u, u_i  &
                                 , F_c, F_d1, F_d2, F_d3, F_p, F_s )
        end select
      end do

      ! assembly ...............................................................

      do i =1, ns
        tau = b(i) * dt
        call MergeArrays(ONE, u, tau, F_c (:,:,:,:,:,i), multi=.true.)
        call MergeArrays(ONE, u, tau, F_d1(:,:,:,:,:,i), multi=.true.)
        call MergeArrays(ONE, u, tau, F_d2(:,:,:,:,:,i), multi=.true.)
        call MergeArrays(ONE, u, tau, F_p (:,:,:,:,:,i), multi=.true.)
        call MergeArrays(ONE, u, tau, F_s (:,:,:,:,:,i), multi=.true.)
      end do

      ! enforce continuity:  ∇²δp = ∇⋅ṽ/∆t, v + J(v) = ṽ - ∆t∇δp
      if (this % project > 0) then
        associate(dp => F_p(:,:,:,:,4,1), w => F_s(:,:,:,:,:,1))
          call PressureSolver(problem, flow_op, dt, u, dp, w=w)
          call ProjectionStep(problem, flow_op, dt, dp, u, w=w)
        end associate
      end if

      call SetArray(p, u_i(:,:,:,:,4), multi=.true.)

      ! prepare next timestep ..................................................

      t = t + dt

      ! clean up ...............................................................

      !$omp barrier
      !$omp master
      if(allocated(ts   ))   deallocate(ts  )
      if(allocated(u_i  ))   deallocate(u_i )
      if(allocated(F_c  ))   deallocate(F_c )
      if(allocated(F_d1 ))   deallocate(F_d1)
      if(allocated(F_d2 ))   deallocate(F_d2)
      if(allocated(F_d3 ))   deallocate(F_d3)
      if(allocated(F_p  ))   deallocate(F_p )
      if(allocated(F_s  ))   deallocate(F_s )
      if(allocated(nu   ))   deallocate(nu  )
      !$omp end master

    end associate

  end subroutine TimeStep

  !-----------------------------------------------------------------------------
  !> Executes one IMEX Runge-Kutta stage

  subroutine RungeKuttaStage_v2( this, i, t, dt, nu, u_0, u_i    &
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
    real(RNP), intent(inout) :: F_d3(:,:,:,:,:,:) !< F_d3 = -χ∇(ν∇⋅v)
    real(RNP), intent(inout) :: F_p (:,:,:,:,:,:) !< F_p  = -∇p
    real(RNP), intent(inout) :: F_s (:,:,:,:,:,:) !< F_s  =  f(x,t)

    optional :: nu

    ! local variables...........................................................

    real(RNP), allocatable, save :: w (:,:,:,:,:)  ! workspace for 2nd and 3rd step
    real(RNP), allocatable, save :: dp(:,:,:,:  )  ! workspace for final step

    real(RNP) :: tau
    integer   :: j

    ! workspace ................................................................

    !$omp single
    allocate(w , mold = u_0)
    allocate(dp, mold = u_0(:,:,:,:,4))
    !$omp end single

    associate( problem => this % problem           &
             , flow_op => this % flow_op           &
             , mesh    => this % flow_op % mesh    &
             , x       => this % flow_op % x       &
             , p       => u_i(:,:,:,:,4)           &
             , a_ex    => this % imex_rk % a_ex    &
             , a_im    => this % imex_rk % a_im    &
             , ns      => this % imex_rk % n_stage &
             , eop     => this % flow_op % eop_u   &
             )

      ! initialization .........................................................

      call GetBoundaryValues(problem, mesh, flow_op % bv_x, t, flow_op % bv_u)

      call SetArray(u_i, u_0, multi=.true.)
      call SetArray(F_d1(:,:,:,:,:,i), ZERO, multi=.true.)

      ! extrapolation: u_i ← u .................................................

      do j = 1, i-1

        tau = dt * a_ex(i,j)
        call MergeArrays(ONE, u_i, tau, F_c (:,:,:,:,:,j), multi =.true.)
        call MergeArrays(ONE, u_i, tau, F_d2(:,:,:,:,:,j), multi =.true.)

        tau = dt * a_im(i,j)
        call MergeArrays(ONE, u_i, tau, F_p (:,:,:,:,:,j), multi = .true.)
        call MergeArrays(ONE, u_i, tau, F_s (:,:,:,:,:,j), multi = .true.)
        call MergeArrays(ONE, u_i, tau, F_d1(:,:,:,:,:,j), multi =.true.)

        tau = dt * (a_ex(i,j) - a_im(i,j)) / a_im(i,i)
        call MergeArrays( ONE, F_d1(:,:,:,:,:,i) &
                        , tau, F_d1(:,:,:,:,:,j) &
                        , multi =.true.)
      end do

      ! projection: u_i ← u", p ← p" ...........................................

      tau = dt * a_im(i,i)

      call PressureSolver(problem, flow_op, tau, u_i, p, w) ! ∇²p" = ∇⋅v'/τ
      call ProjectionStep(problem, flow_op, tau, p, u_i, w) ! v"+J(v") = v'-τ∇p"

      ! diffusion: u_i ← u''' ..................................................

      if (present(nu)) then
        call problem % GetDiffusivity(flow_op%x, t, u_i, nu) ! ν = ν(x,t,u")
      end if

      ! RHS: f = uᵢ - τ F_d1(*,i) → F_d1(*,i)
      call MergeArrays(tau, F_d1(:,:,:,:,:,i), ONE, u_i, multi=.true.)

      ! uᵢ/τ - ∇⋅(ν∇uᵢ) = f/τ
      call DiffusionStep(problem, flow_op, tau, F_d1(:,:,:,:,:,i), u_i, w, nu)

      ! final projection: u_i ← u ..............................................

      if (flow_op % control % div_final) then

        ! ∇²δp = ∇⋅v'''/τ
        call SetArray(dp, ZERO)
        call PressureSolver(problem, flow_op, tau, u_i, dp, w)
        call MergeArrays(ONE, p, ONE, dp)

        ! v + J(v) = v''' - 1/τ ∇p
        call ProjectionStep(problem, flow_op, tau, dp, u_i(:,:,:,:,1:3), w)

      end if

      ! RHS contributions ......................................................

      if (present(nu)) then
        if (problem % HasVariableProperties()) then
          call problem % GetDiffusivity(flow_op%x, t, u_i, nu) ! ν = ν(x,t,u")
        end if
      end if

      ! contributions to time derivative
      call TimeDerivative( problem, flow_op, t        &
                         , u_c  = u_i                 &
                         , u_d  = u_i                 &
                         , p    = p                   &
                         , nu   = nu                  &
                         , F_c  = F_c  (:,:,:,:,:,i)  &
                         , F_d1 = F_d1 (:,:,:,:,:,i)  &
                         , F_d2 = F_d2 (:,:,:,:,:,i)  &
                         , F_d3 = F_d3 (:,:,:,:,:,i)  &
                         , F_p  = F_p  (:,:,:,:,:,i)  &
                         , F_s  = F_s  (:,:,:,:,:,i)  &
                         )

    end associate

    ! clean up .................................................................

    !$omp barrier
    !$omp master
    if (allocated(w )) deallocate(w )
    if (allocated(dp)) deallocate(dp)
    !$omp end master

  end subroutine RungeKuttaStage_v2

  !=============================================================================
  ! TimeIntegrator_RungeKutta_Options: type-bound procedures

  !-----------------------------------------------------------------------------
  !> MPI broadcasting of time-integrator options

  subroutine TimeIntegrator_RungeKutta_Options_Bcast(this, root, comm)
    class(TimeIntegrator_RungeKutta_Options), intent(inout) :: this
    integer,        intent(in) :: root !< rank of broadcast root
    type(MPI_Comm), intent(in) :: comm !< MPI communicator

    type(MPI_Request)  :: request(5)
    type(MPI_Status)   :: stat(size(request))
    integer :: n

    n = 1
    call XMPI_Ibcast( this % n_stage , root, comm, request(n) );  n = n + 1
    call XMPI_Ibcast( this % method  , root, comm, request(n) );  n = n + 1
    call XMPI_Ibcast( this % pressure, root, comm, request(n) );  n = n + 1
    call XMPI_Ibcast( this % project , root, comm, request(n) );  n = n + 1
    call XMPI_Ibcast( this % variant , root, comm, request(n) )

    call MPI_Waitall(n, request, stat)

  end subroutine TimeIntegrator_RungeKutta_Options_Bcast

  !===========================================================================

end module CART__ISP_Flow__Time_Integrator__Runge_Kutta
