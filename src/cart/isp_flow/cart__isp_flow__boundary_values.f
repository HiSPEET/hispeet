!> summary:  Boundary values for incompressible single-phase flow
!> author:   Joerg Stiller
!> date:     2018/04/16
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!=============================================================================

module CART__ISP_Flow__Boundary_Values
  use Kind_Parameters, only: RNP
  use ISP_Flow_Problem
  use CART__Mesh_Partition
  use CART__Boundary_Variable

  implicit none
  private

  public :: GetBoundaryPoints
  public :: GetBoundaryValues
  public :: GetBoundaryTimeDerivative

contains

!-----------------------------------------------------------------------------
!> Extracts the boundary mesh points into a boundary variable
!>
!> Note that the `inout` intent prevents race conditions because of concurrent
!> deallocation with OpenMP and allows to reuse allocated components of
!> matching size.

subroutine GetBoundaryPoints(mesh, x, bv_x)
  class(MeshPartition),   intent(in)    :: mesh         !< mesh partition
  real(RNP),              intent(in)    :: x(:,:,:,:,:) !< mesh points
  type(BoundaryVariable), intent(inout) :: bv_x(:)      !< boundary points

  integer :: i

  do i = 1, mesh % n_boundary
    call bv_x(i) % Extract(mesh, x, i)
  end do

end subroutine GetBoundaryPoints

!-----------------------------------------------------------------------------
!> Extracts the boundary values of all variables at given time
!>
!> Note that the `inout` intent prevents race conditions because of concurrent
!> deallocation with OpenMP and allows to reuse allocated components of
!> matching size.

subroutine GetBoundaryValues(problem, mesh, bv_x, t, bv_u)
  class(FlowProblem),     intent(in)    :: problem !< flow problem
  class(MeshPartition),   intent(in)    :: mesh    !< mesh partition
  type(BoundaryVariable), intent(in)    :: bv_x(:) !< boundary points
  real(RNP),              intent(in)    :: t       !< time
  type(BoundaryVariable), intent(inout) :: bv_u(:) !< boundary values

  real(RNP), contiguous, pointer, save :: xb(:,:,:,:)
  real(RNP), contiguous, pointer, save :: ub(:,:,:,:)
  integer :: i

  do i = 1, mesh % n_boundary

    if (all(problem % bc(i,:) == 'P')) cycle
    !$omp barrier
    !$omp single
    xb => bv_x(i) % Components()
    ub => bv_u(i) % Components()
    !$omp end single
    call problem % GetBoundaryValues(i, xb, t, ub)

  end do

end subroutine GetBoundaryValues

!-----------------------------------------------------------------------------
!> Extracts the boundary time derivative of all variables at given time
!>
!> Note that the `inout` intent prevents race conditions because of concurrent
!> deallocation with OpenMP and allows to reuse allocated components of
!> matching size.

subroutine GetBoundaryTimeDerivative(problem, mesh, bv_x, t, bv_dt_u)
  class(FlowProblem),     intent(in)    :: problem    !< flow problem
  class(MeshPartition),   intent(in)    :: mesh       !< mesh partition
  type(BoundaryVariable), intent(in)    :: bv_x(:)    !< boundary points
  real(RNP),              intent(in)    :: t          !< time
  type(BoundaryVariable), intent(inout) :: bv_dt_u(:) !< boundary values

  real(RNP), contiguous, pointer, save ::    xb(:,:,:,:)
  real(RNP), contiguous, pointer, save :: dt_ub(:,:,:,:)
  integer :: i

  do i = 1, mesh % n_boundary

    if (all(problem % bc(i,:) == 'P')) cycle
    !$omp single
    xb    => bv_x   (i) % Components()
    dt_ub => bv_dt_u(i) % Components()
    !$omp end single
    call problem % GetBoundaryTimeDerivative(i, xb, t, dt_ub)

  end do

end subroutine GetBoundaryTimeDerivative

!=============================================================================

end module CART__ISP_Flow__Boundary_Values
