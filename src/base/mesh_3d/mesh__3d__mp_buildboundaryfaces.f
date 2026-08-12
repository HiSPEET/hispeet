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

!> summary:  Identification of boundary faces
!> author:   Joerg Stiller
!> date:     2022/07/11
!===============================================================================

submodule(Mesh__3D) MP_BuildBoundaryFaces
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Generation of boundary faces
  !>
  !> Requires
  !>   - `mesh % n_elem`
  !>   - `mesh % n_boundary`
  !>   - `mesh % boundary` without faces
  !>   - `mesh % element % face % boundary`
  !>
  !> Generates
  !>   - `mesh % boundary % n_face`
  !>   - `mesh % boundary % face`

  module subroutine BuildBoundaryFaces(mesh)
    class(Mesh_3D), intent(inout) :: mesh !< mesh partition

    ! local data ...............................................................

    integer :: f_bound(mesh % n_bound)
    integer :: b, e, i

    ! count local boundary faces ...............................................

    f_bound = 0

    do e = 1, mesh % n_elem
      associate(face => mesh % element(e) % face)
        do i = 1, 6
          b = face(i) % boundary
          if (b > 0) then
            f_bound(b) = f_bound(b) + 1
          end if
        end do
      end associate
    end do

    ! create local boundary faces ..............................................

    do b = 1, mesh % n_bound
      mesh % boundary(b) % n_face = f_bound(b)
      allocate(mesh % boundary(b) % face( f_bound(b) ))
    end do

    if (mesh%n_elem == 0) return

    f_bound = 0

    do e = 1, mesh % n_elem
      associate(face => mesh % element(e) % face)
        do i = 1, 6
          b = face(i) % boundary
          if (b > 0) then
            f_bound(b) = f_bound(b) + 1
            mesh % boundary(b) % face(f_bound(b)) % element_id   = e
            mesh % boundary(b) % face(f_bound(b)) % element_face = i
          end if
        end do
      end associate
    end do

  end subroutine BuildBoundaryFaces

  !=============================================================================

end submodule MP_BuildBoundaryFaces
