!> summary:  IMEX Euler with velocity or pressure correction scheme
!> author:   Joerg Stiller
!> date:     2018/04/10
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### IMEX Euler with velocity or pressure correction scheme
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
  public :: EulerPC

contains

!-------------------------------------------------------------------------------
!> Velocity correction

subroutine EulerVC(problem, flow_op, t, dt, u_0, u, F, G, S)

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
    call SetArray(u_i, u_0, multi=.true.)

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

    call PressureSolver(problem, flow_op, dt, u_i, p, w)

    ! u_i -= dt (∇p + F_d)
    call WeakGradient(mesh, eop%w, eop%D, p, w)
    call MergeArrays(ONE, u_i, -dt, w  , multi=.true. )
    call MergeArrays(ONE, u_i, -dt, F_d, multi=.true.)

    call DiffusionStep(problem, flow_op, dt, f=u_i, u=u, w=w)

    ! SDC: time derivative .....................................................

    if (present(F) .and. present(G)) then

      call TimeDerivative(problem, flow_op, t, u,   u, p, F = F)
      call TimeDerivative(problem, flow_op, t, u_0, u, p, F = G)

    end if

    ! clean-up .................................................................

    !$omp barrier
    !$omp master
    deallocate(F_d, u_i, w)
    !$omp end master

  end associate

end subroutine EulerVC

!-------------------------------------------------------------------------------
!> Pressure correction

subroutine EulerPC(problem, flow_op, t, dt, u_0, u, F, G, S)

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
    call SetArray(u_i, u_0, multi=.true.)

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

      call TimeDerivative(problem, flow_op, t, u,   u, p, F = F)
      call TimeDerivative(problem, flow_op, t, u_0, u, p, F = G)

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
