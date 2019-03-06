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
  use CART__ISP_Flow__Operators

  implicit none
  private

  public :: PressureSolver_MO

contains

!-------------------------------------------------------------------------------
!>  Pressure solver using IP/DG-SEM with mixed-order approximation

subroutine PressureSolver_MO(problem, flow_op, dt, v_i, p, w, i_max)

  ! arguments ..................................................................

  class(FlowProblem),   intent(in)    :: problem        !< flow problem
  class(FlowOperators), intent(inout) :: flow_op        !< flow operators
  real(RNP),            intent(in)    :: dt             !< time-step size
  real(RNP),            intent(in)    :: v_i(:,:,:,:,:) !< v*        @ po_u
  real(RNP),            intent(inout) :: p(:,:,:,:)     !< pressure  @ po_u
  real(RNP),            intent(out)   :: w(:,:,:,:,:)   !< workspace @ po_u
  integer,    optional, intent(in)    :: i_max          !< num MG/CG cycles

  ! local variables  ...........................................................

  real(RNP), allocatable :: q(:,:,:,:)                 ! pressure  @ po_q = po_p
  real(RNP), allocatable :: f(:,:,:,:)                 ! RHS for q
  type(BoundaryVariable), allocatable, save :: bv_q(:) ! boundary values for q

  real(RNP) :: r_2
  integer   :: po_q, ne, ni
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

    po_q = flow_op % po_p

    !$omp single

    ! pressure and RHS @ po_q
    allocate(q(0:po_q, 0:po_q, 0:po_q, ne))
    allocate(f, mold=q)

    ! pressure boundary values @ po_q
    allocate(bv_q(mesh % n_boundary))
    do b = 1, mesh % n_boundary
      bv_q(b) = BoundaryVariable(mesh, po_q, b, bc(b))
    end do

    !$omp end single

    ! interpolate initial values from p to q
    call flow_op % iop_up % Apply(p, q)

    ! compute divergence of v_i and project to f
    call WeakDivergence(mesh, eop_u%w, eop_u%D, v_i, div_v)  ! div_v = ∇·v*
    call pop_up % Apply(div_v, f)                            ! f = P(div_v)

    ! compute and apply boundary conditions
    call GetBC(problem, mesh, pop_up%A, dt, flow_op%bv_u, v_i, bv_q)
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

end subroutine PressureSolver_MO

!-------------------------------------------------------------------------------
!> Implied pressure boundary conditions

subroutine GetBC(problem, mesh, A, dt, bv_v, v_i, bv_q)

  ! arguments ..................................................................

  class(FlowProblem),      intent(in)    :: problem  !< flow problem
  class(MeshPartition),    intent(in)    :: mesh     !< mesh partition
  real(RNP),               intent(in)    :: A(0:,0:) !< 1D projection operator
  real(RNP),               intent(in)    :: dt       !< time step size
  class(BoundaryVariable), intent(in)    :: bv_v(:)  !< v(t) @ boundary
  real(RNP),               intent(in)    :: v_i      !< v*
  class(BoundaryVariable), intent(inout) :: bv_q(:)  !< BC for q(t)

  dimension :: v_i(0:,0:,0:,:,:)

  ! local variables  ...........................................................

  real(RNP), pointer, contiguous, save :: v_b(:,:,:,:), dn_q(:,:,:)
  real(RNP), allocatable :: dn_p(:,:), z(:,:)
  real(RNP) :: c, tmp
  integer   :: nb, po_q, po_v
  integer   :: b, e, i, j, k, l

  ! initialization .............................................................

  if (all(problem % bc == 'P')) return

  po_q = ubound(A,1)       ! polynomial order of q
  po_v = ubound(A,2)       ! polynomial order of v*
  nb   = mesh % n_boundary ! number of boundaries

  allocate(dn_p(0:po_v, 0:po_v), z(0:po_q, 0:po_v))

  c = ONE / dt

  do b = 1, nb

    if (problem % bc(b,4) == 'P') cycle

    !$omp single
    v_b  => bv_v(b) % Components(1,3)
    dn_q => bv_q(b) % Component(1)
    !$omp end single

    associate(boundary => mesh % boundary(b))

      !$omp do private(dn_p, z)
      do l = 1, boundary % nf
        e = boundary % face(l) % mesh_element % id

        ! ∂p/∂n ................................................................

        select case(boundary % face(l) % mesh_element % face)
        case(1)
          i = 0
          do k = 0, po_v
          do j = 0, po_v
            dn_p(j,k) =  (c*(v_b(j,k,l,1) - v_i(i,j,k,e,1)))
          end do
          end do
        case(2)
          i = po_v
          do k = 0, po_v
          do j = 0, po_v
            dn_p(j,k) = -(c*(v_b(j,k,l,1) - v_i(i,j,k,e,1)))
          end do
          end do
        case(3)
          j = 0
          do k = 0, po_v
          do i = 0, po_v
            dn_p(i,k) =  (c*(v_b(i,k,l,2) - v_i(i,j,k,e,2)))
          end do
          end do
        case(4)
          j = po_v
          do k = 0, po_v
          do i = 0, po_v
            dn_p(i,k) = -(c*(v_b(i,k,l,2) - v_i(i,j,k,e,2)))
          end do
          end do
        case(5)
          k = 0
          do j = 0, po_v
          do i = 0, po_v
            dn_p(i,j) =  (c*(v_b(i,j,l,3) - v_i(i,j,k,e,3)))
          end do
          end do
        case(6)
          k = po_v
          do j = 0, po_v
          do i = 0, po_v
            dn_p(i,j) = -(c*(v_b(i,j,l,3) - v_i(i,j,k,e,3)))
          end do
          end do
        end select

        ! ∂q/∂n = P(∂p/∂n) .....................................................

        ! project along direction 1
        do j = 0, po_v
        do i = 0, po_q
          tmp = 0
          do k = 0, po_v
            tmp = tmp + A(i,k) * dn_p(k,j)
          end do
          z(i,j) = tmp
        end do
        end do

        ! project along direction 2
        do j = 0, po_q
        do i = 0, po_q
          tmp = 0
          do k = 0, po_v
            tmp = tmp + A(j,k) * z(i,k)
          end do
          dn_q(i,j,l) = tmp
        end do
        end do

      end do

    end associate
  end do

end subroutine GetBC

!===============================================================================

end module CART__ISP_Flow__Pressure__MO

