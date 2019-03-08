!> summary:  ISP flow: pressure step using DG with mixed-order approximation
!> author:   Joerg Stiller
!> date:     2018/04/17
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### ISP flow: pressure step using DG with mixed-order approximation
!===============================================================================

module CART__ISP_Flow__Pressure__MO

  use Kind_Parameters, only: RNP
  use Constants,       only: ONE
  use ISP_Flow_Problem
  use CART__Mesh_Partition
  use CART__Boundary_Variable
  use CART__Weak_Divergence
  use CART__Elliptic_PMG
  use CART__ISP_Flow__Boundary_Values
  use CART__ISP_Flow__Operators

  implicit none
  private

  public :: PressureSolver_MO

  interface PressureSolver_MO
    module procedure PressureSolver_IBC
    module procedure PressureSolver_CBC
  end interface

contains

!-------------------------------------------------------------------------------
!> Mixed-order pressure solver with implied boundary conditions

subroutine PressureSolver_IBC(problem, flow_op, dt, v_i, p, w, i_max)

  ! arguments ..................................................................

  class(FlowProblem),   intent(in)    :: problem         !< flow problem
  class(FlowOperators), intent(inout) :: flow_op         !< flow operators
  real(RNP),            intent(in)    :: dt              !< time-step size
  real(RNP),            intent(in)    :: v_i(:,:,:,:,:)  !< ṽ         @ po_u
  real(RNP),            intent(inout) :: p(:,:,:,:)      !< pressure  @ po_u
  real(RNP),            intent(out)   :: w(:,:,:,:,:)    !< workspace @ po_u
  integer,    optional, intent(in)    :: i_max           !< num MG/CG cycles

  ! local variables  ...........................................................

  real(RNP), allocatable :: q(:,:,:,:)                 ! pressure  @ pq = po_p
  real(RNP), allocatable :: f(:,:,:,:)                 ! RHS for q
  type(BoundaryVariable), allocatable, save :: bv_q(:) ! boundary values for q

  real(RNP) :: r_2
  integer   :: pq, ne, ni
  integer   :: b

  if (present(i_max)) then
    if (i_max <= 0) return
  end if

  associate( mesh   => flow_op % mesh      &
           , eop_u  => flow_op % eop_u     &
           , pop_up => flow_op % pop_up    &
           , pmg    => flow_op % pmg_p     &
           , bc     => problem % bc(:,4)   &
           , div_v  => w(:,:,:,:,1)        &
           )

    ! preliminaries ............................................................

    pq = flow_op % po_p

    !$omp single

    ! pressure and RHS @ pq
    allocate(q(0:pq, 0:pq, 0:pq, ne))
    allocate(f, mold=q)

    ! pressure boundary values @ pq
    allocate(bv_q(mesh % n_boundary))
    do b = 1, mesh % n_boundary
      bv_q(b) = BoundaryVariable(mesh, pq, b, bc(b))
    end do

    !$omp end single

    ! interpolate initial values from p to q
    call flow_op % iop_up % Apply(p, q)

    ! compute divergence of v_i and project to f
    call WeakDivergence(mesh, eop_u%w, eop_u%D, v_i, div_v)  ! div_v = ∇·v*
    call pop_up % Apply(div_v, f)                            ! f = P(div_v)

    ! compute and apply boundary conditions
    call GetImpliedBC(problem, mesh, pop_up%A, dt, flow_op%bv_u, v_i, bv_q)
    call pmg % BcToRHS(bv_q, f)

    ! pressure .................................................................

    ! compute
    call pmg % MG_CG_Solver(q, f, ni, r_2, i_max)

    ! monitoring
    !$omp single
    if (flow_op % monitor_level > 0 .and. mesh % part == 0) then
      print '(4X,A,I4,A,ES9.2)', 'pressure  p:    ni =', ni, ', ‖r‖ =', r_2
    end if
    !$omp end single

    ! interpolate to velocity space: p = I(q)
    call flow_op % iop_pu % Apply(q, p)

    ! clean-up .................................................................

    !$omp barrier
    !$omp master
    deallocate(q, f, bv_q)
    !$omp end master

  end associate

end subroutine PressureSolver_IBC

!-------------------------------------------------------------------------------
!> Implied pressure boundary conditions

subroutine GetImpliedBC(problem, mesh, A, dt, bv_v, v_i, bv_q)

  ! arguments ..................................................................

  class(FlowProblem),      intent(in)    :: problem  !< flow problem
  class(MeshPartition),    intent(in)    :: mesh     !< mesh partition
  real(RNP),               intent(in)    :: A(:,:)   !< 1D projection operator
  real(RNP),               intent(in)    :: dt       !< time step size
  class(BoundaryVariable), intent(in)    :: bv_v(:)  !< v(t) @ boundary
  real(RNP),               intent(in)    :: v_i      !< ṽ
  class(BoundaryVariable), intent(inout) :: bv_q(:)  !< BC for q(t)

  dimension :: v_i(0:,0:,0:,:,:)

  ! local variables  ...........................................................

  real(RNP), pointer, contiguous, save :: v_b(:,:,:,:), dn_q(:,:,:)
  real(RNP), allocatable :: dn_p(:,:), w(:,:)
  real(RNP) :: c

  integer :: pq, pv
  integer :: nb, nq, nv
  integer :: b, e, i, j, k, l

  ! initialization .............................................................

  if (all(problem % bc == 'P')) return

  nb = mesh % n_boundary ! number of boundaries
  nq = size(A,1)         ! number of pressure points per direction
  pq = nq - 1            ! polynomial order of q
  nv = size(A,2)         ! number of velocity points per direction
  pv = nv - 1            ! polynomial order of ṽ

  allocate(dn_p(0:pv, 0:pv), w(nq,nv))

  c = ONE / dt

  do b = 1, nb

    if (problem % bc(b,4) == 'P') cycle

    !$omp single
    v_b  => bv_v(b) % Components(1,3)
    dn_q => bv_q(b) % Component(1)
    !$omp end single

    associate(boundary => mesh % boundary(b))

      !$omp do private(dn_p, w)
      do l = 1, boundary % nf
        e = boundary % face(l) % mesh_element % id

        select case(boundary % face(l) % mesh_element % face)
        case(1)
          i = 0
          do k = 0, pv
          do j = 0, pv
            dn_p(j,k) =  (c*(v_b(j,k,l,1) - v_i(i,j,k,e,1)))
          end do
          end do
        case(2)
          i = pv
          do k = 0, pv
          do j = 0, pv
            dn_p(j,k) = -(c*(v_b(j,k,l,1) - v_i(i,j,k,e,1)))
          end do
          end do
        case(3)
          j = 0
          do k = 0, pv
          do i = 0, pv
            dn_p(i,k) =  (c*(v_b(i,k,l,2) - v_i(i,j,k,e,2)))
          end do
          end do
        case(4)
          j = pv
          do k = 0, pv
          do i = 0, pv
            dn_p(i,k) = -(c*(v_b(i,k,l,2) - v_i(i,j,k,e,2)))
          end do
          end do
        case(5)
          k = 0
          do j = 0, pv
          do i = 0, pv
            dn_p(i,j) =  (c*(v_b(i,j,l,3) - v_i(i,j,k,e,3)))
          end do
          end do
        case(6)
          k = pv
          do j = 0, pv
          do i = 0, pv
            dn_p(i,j) = -(c*(v_b(i,j,l,3) - v_i(i,j,k,e,3)))
          end do
          end do
        end select

        call ProjectToPressureSpace(nq, nv, A, dn_p, dn_q(:,:,l), w)

      end do

    end associate
  end do

end subroutine GetImpliedBC

!-------------------------------------------------------------------------------
!> Mixed-order pressure solver with consistent boundary conditions

subroutine PressureSolver_CBC(problem, flow_op, F_v, t, p, w, i_max)

  ! arguments ..................................................................

  class(FlowProblem),   intent(in)    :: problem         !< flow problem
  class(FlowOperators), intent(inout) :: flow_op         !< flow operators
  real(RNP),            intent(in)    :: F_v(:,:,:,:,:)  !< ∂ṽ/∂t     @ po_u
  real(RNP),            intent(in)    :: t               !< time
  real(RNP),            intent(inout) :: p(:,:,:,:)      !< pressure  @ po_u
  real(RNP),            intent(out)   :: w(:,:,:,:,:)    !< workspace @ po_u
  integer,    optional, intent(in)    :: i_max           !< num MG/CG cycles

  ! local variables  ...........................................................

  real(RNP), allocatable :: q(:,:,:,:)                 ! pressure  @ pq = po_p
  real(RNP), allocatable :: f(:,:,:,:)                 ! RHS for q
  type(BoundaryVariable), allocatable, save :: bv_q(:) ! boundary values for q

  real(RNP) :: r_2
  integer   :: pq, ne, ni
  integer   :: b

  if (present(i_max)) then
    if (i_max <= 0) return
  end if

  associate( mesh   => flow_op % mesh      &
           , eop_u  => flow_op % eop_u     &
           , pop_up => flow_op % pop_up    &
           , pmg    => flow_op % pmg_p     &
           , bc     => problem % bc(:,4)   &
           , div_F  => w(:,:,:,:,1)        &
           )

    ! preliminaries ............................................................

    pq = flow_op % po_p

    !$omp single

    ! pressure and RHS @ pq
    allocate(q(0:pq, 0:pq, 0:pq, ne))
    allocate(f, mold=q)

    ! pressure boundary values @ pq
    allocate(bv_q(mesh % n_boundary))
    do b = 1, mesh % n_boundary
      bv_q(b) = BoundaryVariable(mesh, pq, b, bc(b))
    end do

    !$omp end single

    ! interpolate initial values from p to q
    call flow_op % iop_up % Apply(p, q)

    ! compute divergence of v_i and project to f
    call WeakDivergence(mesh, eop_u%w, eop_u%D, F_v, div_F)  ! div_F = ∇·∂ṽ/∂t
    call pop_up % Apply(div_F, f)                            ! f = P(div_F)

    ! compute and apply boundary conditions
    call GetConsistentBC(problem, mesh, pop_up%A, flow_op%bv_x, t, F_v, bv_q)
    call pmg % BcToRHS(bv_q, f)

    ! pressure .................................................................

    ! compute
    call pmg % MG_CG_Solver(q, f, ni, r_2, i_max)

    ! monitoring
    !$omp single
    if (flow_op % monitor_level > 0 .and. mesh % part == 0) then
      print '(4X,A,I4,A,ES9.2)', 'pressure  p:    ni =', ni, ', ‖r‖ =', r_2
    end if
    !$omp end single

    ! interpolate to velocity space: p = I(q)
    call flow_op % iop_pu % Apply(q, p)

    ! clean-up .................................................................

    !$omp barrier
    !$omp master
    deallocate(q, f, bv_q)
    !$omp end master

  end associate

end subroutine PressureSolver_CBC

!-------------------------------------------------------------------------------
!> Consistent pressure boundary conditions

subroutine GetConsistentBC(problem, mesh, A, bv_x, t, F_v, bv_q)

  ! arguments ..................................................................

  class(FlowProblem),      intent(in)    :: problem  !< flow problem
  class(MeshPartition),    intent(in)    :: mesh     !< mesh partition
  real(RNP),               intent(in)    :: A(:,:)   !< 1D projection operator
  class(BoundaryVariable), intent(in)    :: bv_x(:)  !< boundary points
  real(RNP),               intent(in)    :: t        !< time
  real(RNP),               intent(in)    :: F_v      !< ∂ṽ/∂t
  class(BoundaryVariable), intent(inout) :: bv_q(:)  !< BC for p(t)

  dimension :: F_v(0:,0:,0:,:,:)

  ! local variables  ...........................................................

  type(BoundaryVariable), allocatable, save :: bv_dt_u(:)
  real(RNP), pointer, contiguous, save :: dt_v(:,:,:,:), dn_q(:,:,:)
  real(RNP), allocatable :: dn_p(:,:), w(:,:)

  integer :: pq, pv
  integer :: nb, nq, nv
  integer :: b, e, i, j, k, l

  ! initialization .............................................................

  if (all(problem % bc == 'P')) return

  nb = mesh % n_boundary ! number of boundaries
  nq = size(A,1)         ! number of pressure points per direction
  pq = nq - 1            ! polynomial order of q
  nv = size(A,2)         ! number of velocity points per direction
  pv = nv - 1            ! polynomial order of ṽ

  !$omp single
  allocate(bv_dt_u(nb))
  do b = 1, mesh%n_boundary
    bv_dt_u(b) = BoundaryVariable(mesh, pv, b, problem%bc(b,:))
  end do
  !$omp end single

  call GetBoundaryTimeDerivative(problem, mesh, bv_x, t, bv_dt_u)

  allocate(dn_p(0:pv, 0:pv), w(nq,nv))

  ! ∂p/∂n ......................................................................

  do b = 1, nb

    if (problem % bc(b,4) == 'P') cycle

    !$omp single
    dt_v => bv_dt_u(b) % Components(1,3)
    dn_q => bv_q(b) % Component(1)
    !$omp end single

    associate(boundary => mesh % boundary(b))

      !$omp do private(dn_p, w)
      do l = 1, boundary % nf
        e = boundary % face(l) % mesh_element % id

        select case(boundary % face(l) % mesh_element % face)
        case(1)
          i = 0
          do k = 0, pv
          do j = 0, pv
            dn_p(j,k) =  (dt_v(j,k,l,1) - F_v(i,j,k,e,1))
          end do
          end do
        case(2)
          i = pv
          do k = 0, pv
          do j = 0, pv
            dn_p(j,k) = -(dt_v(j,k,l,1) - F_v(i,j,k,e,1))
          end do
          end do
        case(3)
          j = 0
          do k = 0, pv
          do i = 0, pv
            dn_p(i,k) =  (dt_v(i,k,l,2) - F_v(i,j,k,e,2))
          end do
          end do
        case(4)
          j = pv
          do k = 0, pv
          do i = 0, pv
            dn_p(i,k) = -(dt_v(i,k,l,2) - F_v(i,j,k,e,2))
          end do
          end do
        case(5)
          k = 0
          do j = 0, pv
          do i = 0, pv
            dn_p(i,j) =  (dt_v(i,j,l,3) - F_v(i,j,k,e,3))
          end do
          end do
        case(6)
          k = pv
          do j = 0, pv
          do i = 0, pv
            dn_p(i,j) = -(dt_v(i,j,l,3) - F_v(i,j,k,e,3))
          end do
          end do
        end select

        call ProjectToPressureSpace(nq, nv, A, dn_p, dn_q(:,:,l), w)

      end do

    end associate
  end do

end subroutine GetConsistentBC

!-------------------------------------------------------------------------------
!> Project face data from velocity into pressure space

pure subroutine ProjectToPressureSpace(np, nv, A, v, p, w)
  integer,   intent(in)    :: np       !< number pressure points per direction
  integer,   intent(in)    :: nv       !< number velocity points per direction
  real(RNP), intent(in)    :: A(np,nv) !< 1D projection operator
  real(RNP), intent(in)    :: v(nv,nv) !< variable in velocity space
  real(RNP), intent(out)   :: p(np,np) !< variable in pressure space
  real(RNP), intent(inout) :: w(np,nv) !< workspace for the intermediate result

  real(RNP) :: tmp
  integer   :: i, j, k

  ! projection along direction 1
  do j = 1, nv
  do i = 1, np
    tmp = 1
    do k = 1, nv
      tmp = tmp + A(i,k) * v(k,j)
    end do
    w(i,j) = tmp
  end do
  end do

  ! projection along direction 2
  do j = 1, np
  do i = 1, np
    tmp = 1
    do k = 1, nv
      tmp = tmp + A(j,k) * w(i,k)
    end do
    p(i,j) = tmp
  end do
  end do

end subroutine ProjectToPressureSpace

!===============================================================================

end module CART__ISP_Flow__Pressure__MO

