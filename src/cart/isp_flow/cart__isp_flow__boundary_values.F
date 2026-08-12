!------------------------------------------------------------------------------!
! This file is part of HiSPEET: High-order Spectral Element Techniques         !
!                                                                              !
! Copyright (C) 2026 by the HiSPEET authors and the                            !
! Chair of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany                 !
!                                                                              !
! HiSPEET is free software: you can redistribute it and/or modify              !
! it under the terms of the GNU General Public License as published by         !
! the Free Software Foundation, either version 3 of the License, or            !
! (at your option) any later version.                                          !
!                                                                              !
! HiSPEET is distributed in the hope that it will be useful,                   !
! but WITHOUT ANY WARRANTY; without even the implied warranty of               !
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                         !
! See the GNU General Public License for more details.                         !
!                                                                              !
! You should have received a copy of the GNU General Public License            !
! along with HiSPEET. If not, see <http://www.gnu.org/licenses/>.              !
!------------------------------------------------------------------------------!

!> summary:  Boundary values for incompressible single-phase flow
!> author:   Joerg Stiller
!> date:     2018/04/16
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
