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

!> summary:  3D mesh link
!> author:   Joerg Stiller
!> date:     2020/11/20
!===============================================================================

module Mesh_Link__3D
  use Kind_Parameters, only: IXS
  implicit none
  private

  public :: MeshLink_3D
  public :: MeshElementLink_3D

  !-----------------------------------------------------------------------------
  !> Face linking information

  type MeshFaceLink_3D
    integer :: element_id   = 0  !< local element ID
    integer :: element_face = 0  !< linked element face (1:6)
  end type MeshFaceLink_3D

  !-----------------------------------------------------------------------------
  !> Element linking information
  !>
  !> The structure may refer to a local element or a virtual "ghost" element.
  !> The components `vertex`, `edge` and `face` indicate element components,
  !> for which data needs to be transferred. This specification allows for
  !> implementing partial transfers moving only a selected subset of the given
  !> element data.

  type MeshElementLink_3D
    integer      :: id = 0        !< local/ghost element ID
    integer(IXS) :: vertex(8) = 0 !< set to 1(0) for (un)linked vertices
    integer(IXS) :: edge(12)  = 0 !< set to 1(0) for (un)linked edges
    integer(IXS) :: face(6)   = 0 !< set to 1(0) for (un)linked faces
  end type MeshElementLink_3D

  !-----------------------------------------------------------------------------
  !> Structure for keeping information for coupling mesh partitions
  !>
  !> Provides the following lists:
  !>
  !>   - `face  `: element faces linked with partition `part`
  !>   - `master`: link to local elements with a ghost in partition `part`
  !>   - `ghost `: link to ghost elements with their master in partition `part`

  type MeshLink_3D
    integer :: part     = -1  !< remote partition ID
    integer :: n_face   =  0  !< number of linked faces
    integer :: n_master =  0  !< number of master elements
    integer :: n_ghost  =  0  !< number of ghost elements
    type(MeshFaceLink_3D   ), allocatable :: face(:)   !< linked faces
    type(MeshElementLink_3D), allocatable :: master(:) !< master elements
    type(MeshElementLink_3D), allocatable :: ghost(:)  !< ghost elements
  end type MeshLink_3D

  !=============================================================================

end module Mesh_Link__3D
