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

!> summary:  Cartesian mesh face
!> author:   Joerg Stiller
!> date:     2017/04/03
!>
!>### Cartesian mesh face
!===============================================================================

module CART__Mesh_Face
  implicit none
  private

  public :: MeshFace

  !-----------------------------------------------------------------------------
  !> Cartesian mesh face
  !>
  !> The elements are sorted following to the coordinate direction. For example,
  !> if the normal vector is aligned with the x-direction, element(1) precedes
  !> element(2) on the x-axis. Since the orientation of the mesh faces is known
  !> (see MeshPartition), this convention allows to identify the corresponding
  !> element face: In case of an x-face, element(1) abuts with face 2 and
  !> element(2) with face 1. In case of y- and z-faces these are faces 4 and 3,
  !> 6 and 5, respectively.

  type MeshFace
    integer :: element(2) = 0 !< adjacent element ID, including ghost, 0 if none
  end type MeshFace

!===============================================================================

end module CART__Mesh_Face
