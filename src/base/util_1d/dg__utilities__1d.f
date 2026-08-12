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

!> summary:  Discontinuous 1D spectral element utilities
!> author:   Johannes Stojanow, Joerg Stiller
!> date:     2019/05/15
!>
!>### Discontinuous 1D spectral element utilities
!>
!> Provides common routines for Lobatto-based nodal DG-SEM, including
!>
!>   *  mesh generation (`DG_GetMeshPoints_1D`)
!>
!> Some routines require conditions for the left and right boundary points,
!> which are passed in the character array `bc(1:2)`. The following boundary
!> types are supported:
!>
!>   *  periodic  (`'P'`)
!>   *  Dirichlet (`'D'`)
!>   *  Neumann   (`'N'`)
!>
!> While the latter two can be combined as appropriate, periodic conditions
!> must be always specified on both sides, i.e. `bc(1:2) ='P'`.
!>
!===============================================================================

module DG__Utilities__1D
  use Kind_Parameters,   only: RNP
  use Constants,         only: HALF
  use Standard_Element_Operators__1D
  implicit none
  private

  public :: DG_GetMeshPoints_1D

contains

  !-----------------------------------------------------------------------------
  !> Computes the mesh points for a given interval [a,b]
  !>
  !> Requires the initialized standard operators `sop` providing the collocation
  !> points `sop%x`. The polynomial order `sop%po` must match the upper bound of
  !> the first dimension in `x`. The second dimension of `x` determines the
  !> number of elements that are generated.
  !>
  !> Works with all nodal bases.

  subroutine DG_GetMeshPoints_1D(sop, a, b, dx, x)
    class(StandardElementOperators_1D), intent(in) :: sop !< standard operators
    real(RNP),             intent(in)  :: a        !< left border
    real(RNP),             intent(in)  :: b        !< right border
    real(RNP),             intent(out) :: dx       !< element length
    real(RNP), contiguous, intent(out) :: x(0:,:)  !< mesh points x(0:po,1:ne)

    integer   :: k, ne
    real(RNP) :: xe

    ne = size(x,2)
    dx = (b - a) / ne

    associate(xi => sop%x)
      !$omp do
      do k = 1, ne
         xe = a + (k - HALF) * dx      ! element midpoint
         x(:,k) = xe + HALF * dx * xi  ! transformed Lobatto points
      end do
    end associate

  end subroutine DG_GetMeshPoints_1D

  !=============================================================================

end module DG__Utilities__1D
