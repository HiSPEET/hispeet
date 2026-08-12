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

!> summary:  3D mesh face type
!> author:   Joerg Stiller
!> date:     2020/11/13
!===============================================================================

module Mesh_Face__3D
  implicit none
  private

  public :: MeshFace_3D

  !-----------------------------------------------------------------------------
  !> Adjacent mesh element data

  type AdjacentElement
    integer :: id   = 0 !< element id
    integer :: face = 0 !< corresponding element face {1:6}
  end type AdjacentElement

  !-----------------------------------------------------------------------------
  !> 3D mesh face
  !>
  !> The main purpose is to access the adjacent elements. The latter are sorted
  !> according to the normal direction of the face. In case of a structured mesh
  !> the latter coincides with the corresponding coordinate direction.

  type MeshFace_3D
    type(AdjacentElement) :: element(2) !< adjacent elements, including ghosts
  end type MeshFace_3D

end module Mesh_Face__3D
