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
  use CART__DG_Weak_Divergence
  use CART__Elliptic_PMG
  use CART__ISP_Flow__Boundary_Values
  use CART__ISP_Flow__Operators

  implicit none
  private

  public :: PressureSolver_EO

  interface PressureSolver_EO
    module procedure PressureSolver_IBC
    module procedure PressureSolver_CBC
  end interface

contains

!-------------------------------------------------------------------------------
!> Equal-order pressure solver with implied boundary conditions

subroutine PressureSolver_IBC(problem, flow_op, dt, v_i, p, w, i_max)

  ! arguments ..................................................................

  class(FlowProblem),   intent(in)    :: problem         !< flow problem
  class(FlowOperators), intent(inout) :: flow_op         !< flow operators
  real(RNP),            intent(in)    :: dt              !< time-step size
  real(RNP),            intent(in)    :: v_i(:,:,:,:,:)  !< ṽ
  real(RNP),            intent(inout) :: p(:,:,:,:)      !< pressure
  real(RNP),            intent(out)   :: w(:,:,:,:,:)    !< workspace
  integer,    optional, intent(in)    :: i_max           !< num MG/CG cycles

  ! local variables  ...........................................................

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

    np = size(p,1)
    ne = size(p,4)

    g = -product(mesh%dx) / (8 * dt)

    ! RHS ......................................................................

    call WeakDivergence(mesh, eop%w, eop%D, v_i, div_v)   ! div_v = ∇·ṽ
    call TPO_sDDD_Eval(np, ne, g, eop%w, div_v, f)        ! f = -M div_v

    call GetImpliedBC(problem, mesh, dt, v_i, flow_op%bv_u)
    call pmg % BcToRHS(flow_op%bv_u, 4, f)

    ! pressure .................................................................

    call pmg % MG_CG_Solver(p, f, ni, r_2, i_max)

    ! monitoring
    !$omp single
    if (flow_op % control % monitor > 0 .and. mesh % part == 0) then
      print '(4X,A,I4,A,ES9.2)', 'pressure  p:    ni =', ni, ', ‖r‖ =', r_2
    end if
    !$omp end single

  end associate

end subroutine PressureSolver_IBC

!-------------------------------------------------------------------------------
!> Implied pressure boundary conditions

subroutine GetImpliedBC(problem, mesh, dt, v_i, bv_u)

  ! arguments ..................................................................

  class(FlowProblem),      intent(in)    :: problem           !< flow problem
  class(MeshPartition),    intent(in)    :: mesh              !< mesh partition
  real(RNP),               intent(in)    :: dt                !< time step size
  real(RNP),               intent(in)    :: v_i(0:,0:,0:,:,:) !< ṽ
  class(BoundaryVariable), intent(inout) :: bv_u(:)           !< BC

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

   !$omp single
    v_b  => bv_u(b) % Components(1,3)
    dn_p => bv_u(b) % Component(4)
   !$omp end single

    associate(boundary => mesh % boundary(b))

      !$omp do
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

end subroutine GetImpliedBC

!-------------------------------------------------------------------------------
!>  Equal-order pressure solver with consistent boundary conditions

subroutine PressureSolver_CBC(problem, flow_op, F_v, t, p, w, i_max)

  ! arguments ..................................................................

  class(FlowProblem),   intent(in)    :: problem         !< flow problem
  class(FlowOperators), intent(inout) :: flow_op         !< flow operators
  real(RNP),            intent(in)    :: F_v(:,:,:,:,:)  !< ∂ṽ/∂t
  real(RNP),            intent(in)    :: t               !< time
  real(RNP),            intent(inout) :: p(:,:,:,:)      !< pressure
  real(RNP),            intent(out)   :: w(:,:,:,:,:)    !< workspace
  integer,    optional, intent(in)    :: i_max           !< num MG/CG cycles

  ! local variables  ...........................................................

  real(RNP) :: g, r_2
  integer   :: np, ne, ni

  if (present(i_max)) then
    if (i_max <= 0) return
  end if

  ! initialization .............................................................

  if (all(problem % bc == 'P')) return

  associate( mesh  => flow_op % mesh      &
           , eop   => flow_op % eop_u     &
           , pmg   => flow_op % pmg_p     &
           , bc    => problem % bc(:,4)   &
           , div_F => w(:,:,:,:,1)        &
           , f     => w(:,:,:,:,2)        &
           )

    ! preliminaries ............................................................

    np = size(p,1)
    ne = size(p,4)

    g = -product(mesh%dx) / 8

    ! RHS ......................................................................

    call WeakDivergence(mesh, eop%w, eop%D, F_v, div_F) ! div_F = ∇·∂ṽ/∂t
    call TPO_sDDD_Eval(np, ne, g, eop%w, div_F, f)      ! f = -M div_F

    call GetConsistentBC(problem, mesh, flow_op%bv_x, t, F_v, flow_op%bv_u)
    call pmg % BcToRHS(flow_op%bv_u, 4, f)

    ! pressure .................................................................

    call pmg % MG_CG_Solver(p, f, ni, r_2, i_max)

    ! monitoring
    !$omp single
    if (flow_op % control % monitor > 0 .and. mesh % part == 0) then
      print '(4X,A,I4,A,ES9.2)', 'pressure  p:    ni =', ni, ', ‖r‖ =', r_2
    end if
    !$omp end single

  end associate

end subroutine PressureSolver_CBC

!-------------------------------------------------------------------------------
!> Consistent pressure boundary conditions

subroutine GetConsistentBC(problem, mesh, bv_x, t, F_v, bv_u)

  ! arguments ..................................................................

  class(FlowProblem),      intent(in)    :: problem           !< flow problem
  class(MeshPartition),    intent(in)    :: mesh              !< mesh partition
  class(BoundaryVariable), intent(in)    :: bv_x(:)           !< boundary points
  real(RNP),               intent(in)    :: t                 !< time
  real(RNP),               intent(in)    :: F_v(0:,0:,0:,:,:) !< ∂ṽ/∂t
  class(BoundaryVariable), intent(inout) :: bv_u(:)           !< BC for p(t)

  ! local variables  ...........................................................

  type(BoundaryVariable), allocatable, save :: bv_dt_u(:)
  real(RNP), pointer, contiguous, save :: dt_v(:,:,:,:), dn_p(:,:,:)
  integer :: nb, po
  integer :: b, e, i, j, k, l

  ! initialization .............................................................

  if (all(problem % bc == 'P')) return

  po = ubound(F_v, 1)    ! polynomial order
  nb = mesh % n_boundary ! number of boundaries

  !$omp single
  allocate(bv_dt_u(nb))
  do b = 1, mesh%n_boundary
    bv_dt_u(b) = BoundaryVariable(mesh, po, b, problem%bc(b,:))
  end do
  !$omp end single

  call GetBoundaryTimeDerivative(problem, mesh, bv_x, t, bv_dt_u)

  ! ∂p/∂n ......................................................................

  do b = 1, nb

    if (problem % bc(b,4) == 'P') cycle

    !$omp single
    dt_v => bv_dt_u(b) % Components(1,3)
    dn_p => bv_u(b) % Component(4)
    !$omp end single

    associate(boundary => mesh % boundary(b))

      !$omp do
      do l = 1, boundary % nf
        e = boundary % face(l) % mesh_element % id

        select case(boundary % face(l) % mesh_element % face)
        case(1)
          i = 0
          do k = 0, po
          do j = 0, po
            dn_p(j,k,l) =  (dt_v(j,k,l,1) - F_v(i,j,k,e,1))
          end do
          end do
        case(2)
          i = po
          do k = 0, po
          do j = 0, po
            dn_p(j,k,l) = -(dt_v(j,k,l,1) - F_v(i,j,k,e,1))
          end do
          end do
        case(3)
          j = 0
          do k = 0, po
          do i = 0, po
            dn_p(i,k,l) =  (dt_v(i,k,l,2) - F_v(i,j,k,e,2))
          end do
          end do
        case(4)
          j = po
          do k = 0, po
          do i = 0, po
            dn_p(i,k,l) = -(dt_v(i,k,l,2) - F_v(i,j,k,e,2))
          end do
          end do
        case(5)
          k = 0
          do j = 0, po
          do i = 0, po
            dn_p(i,j,l) =  (dt_v(i,j,l,3) - F_v(i,j,k,e,3))
          end do
          end do
        case(6)
          k = po
          do j = 0, po
          do i = 0, po
            dn_p(i,j,l) = -(dt_v(i,j,l,3) - F_v(i,j,k,e,3))
          end do
          end do
        end select

      end do

    end associate
  end do

  ! clean-up ...................................................................

  !$omp barrier
  !$omp master
  deallocate(bv_dt_u)
  !$omp end master

end subroutine GetConsistentBC

!===============================================================================

end module CART__ISP_Flow__Pressure__EO
