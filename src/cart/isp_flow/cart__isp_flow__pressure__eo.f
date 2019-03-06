!> summary:  ISP flow: pressure step using DG with equal-order approximation
!> author:   Joerg Stiller
!> date:     2018/04/17
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### ISP flow: pressure step using DG with equal-order approximation
!===============================================================================

module CART__ISP_Flow__Pressure__EO

  use Kind_Parameters, only: RNP
  use Constants,       only: ONE
  use Array_Assignments
  use TPO_sDDD
  use ISP_Flow_Problem
  use CART__Mesh_Partition
  use CART__Boundary_Variable
  use CART__Weak_Divergence
  use CART__Elliptic_PMG
  use CART__ISP_Flow__Operators

  implicit none
  private

  public :: PressureSolver_EO

contains

!-------------------------------------------------------------------------------
!>  Pressure solver using IP/DG-SEM with equal-order approximation

subroutine PressureSolver_EO(problem, flow_op, dt, v_i, p, w, i_max)

  ! arguments ..................................................................

  class(FlowProblem),   intent(in)    :: problem        !< flow problem
  class(FlowOperators), intent(inout) :: flow_op        !< flow operators
  real(RNP),            intent(in)    :: dt             !< time-step size
  real(RNP),            intent(in)    :: v_i(:,:,:,:,:) !< v*
  real(RNP),            intent(inout) :: p(:,:,:,:)     !< pressure
  real(RNP),            intent(out)   :: w(:,:,:,:,:)   !< workspace
  integer,    optional, intent(in)    :: i_max          !< num MG/CG cycles

  ! local variables  ...........................................................

  type(BoundaryVariable), allocatable, save :: bv_p(:)
  real(RNP) :: g, r_2
  integer   :: np, ne, ni

  if (present(i_max)) then
    if (i_max <= 0) return
  end if

  associate( mesh  => flow_op % mesh      &
           , eop   => flow_op % eop_u     &
           , pmg   => flow_op % pmg_p     &
           , bc    => problem % bc(:,4)   &
           , div_v => w(:,:,:,:,1)        &
           , f     => w(:,:,:,:,2)        &
           )

    ! preliminaries ............................................................

    ! handle for pressure boundary values
    !$omp single
    allocate(bv_p(mesh % n_boundary))
    call flow_op % bv_u % GetHandle(c = 4, handle = bv_p)
    !$omp end single

    np = size(p,1)
    ne = size(p,4)

    g = -product(mesh%dx) / (8 * dt)

    ! RHS ......................................................................

    call WeakDivergence(mesh, eop%w, eop%D, v_i, div_v)   ! div_v = ∇·v*
    call TPO_sDDD_Eval(np, ne, g, eop%w, div_v, f)        ! f = M div_v

    call GetBC(problem, mesh, dt, flow_op%bv_u, v_i, bv_p)
    call pmg % BcToRHS(bv_p, f)

    ! pressure .................................................................

    call pmg % MG_CG_Solver(p, f, ni, r_2, i_max)

    ! monitoring
    !$omp single
    if (flow_op % monitor_level > 0 .and. mesh % part == 0) then
      print '(4X,A,I4,A,ES9.2)', 'pressure  p:    ni =', ni, ', ‖r‖ =', r_2
    end if
    !$omp end single

    ! clean-up .................................................................

    !$omp barrier
    !$omp master
    deallocate(bv_p)
    !$omp end master

  end associate

end subroutine PressureSolver_EO

!-------------------------------------------------------------------------------
!> Implied pressure boundary conditions

subroutine GetBC(problem, mesh, dt, bv_v, v_i, bv_p)

  ! arguments ..................................................................

  class(FlowProblem),      intent(in)    :: problem           !< flow problem
  class(MeshPartition),    intent(in)    :: mesh              !< mesh partition
  real(RNP),               intent(in)    :: dt                !< time step size
  class(BoundaryVariable), intent(in)    :: bv_v(:)           !< v(t) @ boundary
  real(RNP),               intent(in)    :: v_i(0:,0:,0:,:,:) !< v*
  class(BoundaryVariable), intent(inout) :: bv_p(:)           !< BC for p(t)

  ! local variables  ...........................................................

  real(RNP), pointer, contiguous, save :: v_b(:,:,:,:), dn_p(:,:,:)
  real(RNP) :: c
  integer   :: nb, po
  integer   :: b, e, i, j, k, l

  ! initialization .............................................................

  if (all(problem % bc == 'P')) return

  po = ubound(v_i, 1)    ! polynomial order
  nb = mesh % n_boundary ! number of boundaries

  ! ∂p/∂n ......................................................................

  c = ONE / dt

  do b = 1, nb

    if (problem % bc(b,4) == 'P') cycle

    v_b  => bv_v(b) % Components(1,3)
    dn_p => bv_p(b) % Component(1)

    associate(boundary => mesh % boundary(b))

      do l = 1, boundary % nf
        e = boundary % face(l) % mesh_element % id

        select case(boundary % face(l) % mesh_element % face)
        case(1)
          i = 0
          do k = 0, po
          do j = 0, po
            dn_p(j,k,l) =  (c*(v_b(j,k,l,1) - v_i(i,j,k,e,1)))
          end do
          end do
        case(2)
          i = po
          do k = 0, po
          do j = 0, po
            dn_p(j,k,l) = -(c*(v_b(j,k,l,1) - v_i(i,j,k,e,1)))
          end do
          end do
        case(3)
          j = 0
          do k = 0, po
          do i = 0, po
            dn_p(i,k,l) =  (c*(v_b(i,k,l,2) - v_i(i,j,k,e,2)))
          end do
          end do
        case(4)
          j = po
          do k = 0, po
          do i = 0, po
            dn_p(i,k,l) = -(c*(v_b(i,k,l,2) - v_i(i,j,k,e,2)))
          end do
          end do
        case(5)
          k = 0
          do j = 0, po
          do i = 0, po
            dn_p(i,j,l) =  (c*(v_b(i,j,l,3) - v_i(i,j,k,e,3)))
          end do
          end do
        case(6)
          k = po
          do j = 0, po
          do i = 0, po
            dn_p(i,j,l) = -(c*(v_b(i,j,l,3) - v_i(i,j,k,e,3)))
          end do
          end do
        end select

      end do

    end associate
  end do

end subroutine GetBC

!===============================================================================

end module CART__ISP_Flow__Pressure__EO
