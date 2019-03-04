!> summary:  DGM with implicit-explicit Euler velocity correction scheme
!> author:   Joerg Stiller
!> date:     2018/04/10
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### DGM with implicit-explicit Euler velocity correction scheme
!>
!> @todo
!>   *  clean-up
!>   *  eliminate repeated computation
!> @endtodo
!===============================================================================

module CART__ISP_Flow__Euler

  use Kind_Parameters, only: RNP
  use Constants,       only: ONE
  use Array_Assignments
  use ISP_Flow_Problem
  use CART__Mesh_Partition
  use CART__Boundary_Variable
  use CART__Weak_Divergence
  use CART__Weak_Gradient
  use CART__ISP_Flow__Operators
  use CART__ISP_Flow__Boundary_Values
  use CART__ISP_Flow__Convection
  use CART__ISP_Flow__Diffusion
  use CART__ISP_Flow__Pressure
  use CART__ISP_Flow__Time_Derivative

  implicit none
  private

  public :: EulerVC
  public :: EulerCS
  public :: EulerPC

contains

!===============================================================================
! Velocity correction

!-------------------------------------------------------------------------------
!>

subroutine EulerVC(problem, flow_op, t, dt, u_0, u, F, G, S, n_cpi)

  ! arguments ..................................................................

  class(FlowProblem),   intent(in)    :: problem !< flow problem
  class(FlowOperators), intent(inout) :: flow_op !< flow operators
  real(RNP),            intent(inout) :: t       !< time
  real(RNP),            intent(in)    :: dt      !< time step size
  real(RNP),            intent(in)    :: u_0     !< solution u(t)
  real(RNP),            intent(inout) :: u       !< solution u(t+dt)
  real(RNP),  optional, intent(out)   :: F       !< ∂u/∂t(t+dt)
  real(RNP),  optional, intent(inout) :: G       !< δu/δt(t) → δu/δt(t+dt)
  real(RNP),  optional, intent(in)    :: S       !< ∫∂u/dt over (t,t+dt)
  integer,    optional, intent(in)    :: n_cpi   !< = num consist p-iterations

  dimension :: u_0 (:,:,:,:,:)
  dimension :: u   (:,:,:,:,:)
  dimension :: F   (:,:,:,:,:)
  dimension :: G   (:,:,:,:,:)
  dimension :: S   (:,:,:,:,:)

  ! local variables  ...........................................................

  real(RNP), dimension(:,:,:,:,:), allocatable, save :: F_d, u_i, w

  real(RNP) :: t_0
  integer   :: np, ne, nc

  associate( mesh => flow_op % mesh  &
           , x    => flow_op % x     &
           , eop  => flow_op % eop_u &
           , p    => u(:,:,:,:,4)    )

    ! initialization ...........................................................

    t_0 = t
    t   = t + dt

    ! workspace
    !$omp single
    allocate(F_d, mold = u)
    allocate(u_i, mold = u)
    allocate(w  , mold = u)
    !$omp end single

    np = size(u,1)
    ne = size(u,4)
    nc = size(u,5)

    ! boundary conditions ......................................................

    call GetBoundaryValues(problem, mesh, flow_op%bv_x, t, flow_op%bv_u)

    ! predictor ................................................................

    ! u_i = u_0
    call AssignArray(u_i, u_0, multi=.true.)

    ! SDC: u_i += -dt δu/δt(t) + ∫∂u/∂t dt
    if (present(G) .and. present(S)) then
      call MergeArrays(ONE, u_i, -dt, G, multi=.true.)
      call MergeArrays(ONE, u_i, ONE, S, multi=.true.)
    end if

    ! w= F_c + F_d + F_s
    call TimeDerivative( problem, flow_op, t  &
                       , u_c  = u_0           &
                       , u_d  = u             &
                       , F    = w             &
                       , F_d  = F_d           )

    ! u_i += dt w
    call MergeArrays(ONE, u_i, dt, w, multi=.true.)

    ! pressure, continuity and diffusion .......................................

    if (n_cpi < 0 .and. present(F)) then
      ! consistent: pressure BC using exact ∂u/∂t at Γ, F providing workspace
      call AssignArray(F, u_i, multi=.true.)              ! F = (ũ - u₀)/dt
      call MergeArrays(1/dt, F, -1/dt, u_0, multi=.true.) !   ≈ F(u, p=0)
      call PressureSolver(problem, flow_op, t, F, p, w, consistent = .true.)
    else
      ! standard: pressure BC using approximate ∂u/∂t at Γ
      call PressureSolver(problem, flow_op, dt, u_i, p, w)
    end if

    ! u_i -= dt (∇p + F_d)
    call WeakGradient(mesh, eop%w, eop%D, p, w)
    call MergeArrays(ONE, u_i, -dt, w  , multi=.true. )
    call MergeArrays(ONE, u_i, -dt, F_d, multi=.true.)

    call DiffusionStep(problem, flow_op, dt, f=u_i, u=u, w=w)

    ! SDC: time derivative .....................................................

    if (present(G) .and. present(F)) then

      ! F = -∇·vu + ν(Δu -∇∇·v) + f
      call TimeDerivative(problem, flow_op, t, u_c = u, u_d = u, F = F)

      ! recompute pressure
      if (present(n_cpi)) then

        call PressureSolver( problem, flow_op, t, F, p, w  &
                           , consistent = .true.           &
                           , i_max      =  n_cpi           )

      end if

      ! F -= ∇p
      call WeakGradient(mesh, eop%w, eop%D, p, w)
      call MergeArrays(ONE, F, -ONE, w, multi=.true.)

      ! G
      call TimeDerivative(problem, flow_op, t, u_c = u_0, u_d = u, p=p, F = G)

    end if

    ! clean-up .................................................................

    !$omp barrier
    !$omp master
    deallocate(F_d, u_i, w)
    !$omp end master

  end associate

end subroutine EulerVC

!===============================================================================
! An adaptation of consistent splitting

!-------------------------------------------------------------------------------
!>

subroutine EulerCS(problem, flow_op, t, dt, u_0, u, F, G, S, n_cpi)

  ! arguments ..................................................................

  class(FlowProblem),   intent(in)    :: problem !< flow problem
  class(FlowOperators), intent(inout) :: flow_op !< flow operators
  real(RNP),            intent(inout) :: t       !< time
  real(RNP),            intent(in)    :: dt      !< time step size
  real(RNP),            intent(in)    :: u_0     !< solution u(t)
  real(RNP),            intent(inout) :: u       !< solution u(t+dt)
  real(RNP),  optional, intent(out)   :: F       !< ∂u/∂t(t+dt)
  real(RNP),  optional, intent(inout) :: G       !< δu/δt(t) → δu/δt(t+dt)
  real(RNP),  optional, intent(in)    :: S       !< ∫∂u/dt over (t,t+dt)
  integer,    optional, intent(in)    :: n_cpi   !< = num consist p-iterations

  dimension :: u_0 (:,:,:,:,:)
  dimension :: u   (:,:,:,:,:)
  dimension :: F   (:,:,:,:,:)
  dimension :: G   (:,:,:,:,:)
  dimension :: S   (:,:,:,:,:)

  ! local variables  ...........................................................

  real(RNP), dimension(:,:,:,:,:), allocatable, save :: u_i, w

  real(RNP) :: t_0
  integer   :: np, ne, nc

  associate( mesh => flow_op % mesh  &
           , x    => flow_op % x     &
           , eop  => flow_op % eop_u &
           , p    => u(:,:,:,:,4)    )

    ! initialization ...........................................................

    t_0 = t
    t   = t + dt

    ! workspace
    !$omp single
    allocate(u_i, mold = u)
    allocate(w  , mold = u)
    !$omp end single

    np = size(u,1)
    ne = size(u,4)
    nc = size(u,5)

    ! boundary conditions ......................................................

    call GetBoundaryValues(problem, mesh, flow_op%bv_x, t, flow_op%bv_u)

    ! predictor ................................................................

    ! u_i = u_0
    call AssignArray(u_i, u_0, multi=.true.)

    ! SDC: u_i += -dt δu/δt(t) + ∫∂u/dt
    if (present(G) .and. present(S)) then
      call MergeArrays(ONE, u_i, -dt, G, multi=.true.)
      call MergeArrays(ONE, u_i, ONE, S, multi=.true.)
    end if

    ! w= F_c + F_d + F_s
    call TimeDerivative( problem, flow_op, t  &
                       , u_c  = u_0           &
                       , p    = u(:,:,:,:,4)  &
                       , F    = w             )

    ! u_i += dt w
    call MergeArrays(ONE, u_i, dt, w, multi=.true.)

    ! diffusion ................................................................

    call DiffusionStep(problem, flow_op, dt, f=u_i, u=u, w=w)

    ! pressure .................................................................

    ! w = -∇·vu + ν(Δu -∇∇·v) + f
    call TimeDerivative(problem, flow_op, t, u_c = u, u_d = u, F = w)
    call PressureSolver(problem, flow_op, t, w, p, w = u_i, consistent = .true.)

    ! SDC terms ................................................................

    if (present(G) .and. present(F)) then

      ! F = -∇·vu + ν(Δu -∇∇·v) + f - ∇p = w - ∇p
      call AssignArray(F, w, multi=.true.)
      call WeakGradient(mesh, eop%w, eop%D, p, w)
      call MergeArrays(ONE, F, -ONE, w, multi=.true.)

      ! G
      call TimeDerivative(problem, flow_op, t, u_c = u_0, u_d = u, p=p, F = G)

    end if

    ! clean-up .................................................................

    !$omp barrier
    !$omp master
    deallocate(u_i, w)
    !$omp end master

  end associate

end subroutine EulerCS

!===============================================================================
! Pressure correction

!-------------------------------------------------------------------------------
!>

subroutine EulerPC(problem, flow_op, t, dt, u_0, u, F, G, S, n_cpi)

  ! arguments ..................................................................

  class(FlowProblem),   intent(in)    :: problem !< flow problem
  class(FlowOperators), intent(inout) :: flow_op !< flow operators
  real(RNP),            intent(inout) :: t       !< time
  real(RNP),            intent(in)    :: dt      !< time step size
  real(RNP),            intent(in)    :: u_0     !< solution u(t)
  real(RNP),            intent(inout) :: u       !< solution u(t+dt)
  real(RNP),  optional, intent(out)   :: F       !< ∂u/∂t(t+dt)
  real(RNP),  optional, intent(inout) :: G       !< δu/δt(t) → δu/δt(t+dt)
  real(RNP),  optional, intent(in)    :: S       !< ∫∂u/dt over (t,t+dt)
  integer,    optional, intent(in)    :: n_cpi   !< = num consist p-iterations

  dimension :: u_0 (:,:,:,:,:)
  dimension :: u   (:,:,:,:,:)
  dimension :: F   (:,:,:,:,:)
  dimension :: G   (:,:,:,:,:)
  dimension :: S   (:,:,:,:,:)

  ! local variables  ...........................................................

  real(RNP), dimension(:,:,:,:,:), allocatable, save :: u_i, w

  real(RNP) :: t_0
  integer   :: np, ne, nc

  associate( mesh => flow_op % mesh  &
           , x    => flow_op % x     &
           , eop  => flow_op % eop_u &
           , p    => u(:,:,:,:,4)    )

    ! initialization ...........................................................

    t_0 = t
    t   = t + dt

    ! workspace
    !$omp single
    allocate(u_i  , mold = u)
    allocate(w    , mold = u)
    !$omp end single

    np = size(u,1)
    ne = size(u,4)
    nc = size(u,5)

    ! boundary conditions ......................................................

    call GetBoundaryValues(problem, mesh, flow_op%bv_x, t, flow_op%bv_u)

    ! predictor ................................................................

    ! ũ = u_0
    call AssignArray(u_i, u_0, multi=.true.)

    ! SDC: ũ += -dt δu/δt(t) + ∫∂u/dt
    if (present(G) .and. present(S)) then
      call MergeArrays(ONE, u_i, -dt, G, multi=.true.)
      call MergeArrays(ONE, u_i, ONE, S, multi=.true.)
    end if

    ! w = F_c + F_d + F_s
    call TimeDerivative(problem, flow_op, t, u_c=u_0, F=w)

    ! ũ += dt w
    call MergeArrays(ONE, u_i, dt, w, multi=.true.)

    ! diffusion ................................................................

    call DiffusionStep(problem, flow_op, dt, f=u_i, u=u, w=w)

    ! pressure correction ......................................................

    call PressureSolver(problem, flow_op, dt, u, p, w)

    call WeakGradient  (mesh, eop%w, eop%D, p, w)             ! w(*,1:3) = ∇p
    call WeakDivergence(mesh, eop%w, eop%D, u, w(:,:,:,:,4))  ! w(*, 4 ) = ∇·ṽ

    ! p -= ν∇·v
    call MergeArrays(ONE, p, -problem%nu_ref(1), w(:,:,:,:,4))

    ! v -= dt ∇p
    call MergeArrays( ONE, u, -dt, w, multi=.true.)

    ! SDC: time derivative .....................................................

    if (present(G) .and. present(F)) then

      ! F = -∇·vu + ν(Δu -∇∇·v) + f
      call TimeDerivative(problem, flow_op, t, u_c = u, u_d = u, F = F)

      ! recompute pressure
      if (present(n_cpi)) then
        if (n_cpi > 0) then

          call PressureSolver( problem, flow_op, t, F, p, w  &
                             , consistent = .true.           &
                             , i_max      =  n_cpi           )

        end if
      end if

      ! F -= ∇p
      call WeakGradient(mesh, eop%w, eop%D, p, w)
      call MergeArrays(ONE, F, -ONE, w, multi=.true.)

      ! G
      call TimeDerivative(problem, flow_op, t, u_c = u_0, u_d = u, p=p, F = G)

    end if

    ! clean-up .................................................................

    !$omp barrier
    !$omp master
    deallocate(u_i, w)
    !$omp end master

  end associate

end subroutine EulerPC

!===============================================================================

end module CART__ISP_Flow__Euler
