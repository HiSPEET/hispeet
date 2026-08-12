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

!> summary:  Map to child mesh
!> author:   Joerg Stiller
!> date:     2020/11/20
!===============================================================================

module Mesh_Map_To_Child__3D
  use XMPI
  implicit none
  private

  public :: MeshMapToChild_3D

  !-----------------------------------------------------------------------------
  !> Map supporting the transfer of data to and from child elements
  !>
  !> The entries in `id_elem` are sorted according 1) the element ID, 2) the
  !> activity of the corresponding child elements. Elements with active children
  !> precede those with inactive ones.
  !>
  !> @note
  !> The communicator `comm` is identical to `comm_world` in the embedding data
  !> structure `mesh`, but is included to allow communication without access to
  !> the latter.

  type MeshMapToChild_3D
    type(MPI_Comm) :: comm !< MPI communicator, usually comm_world
    integer :: proc        !< MPI rank of process keeping the child partition
    integer :: n_elem      !< num elements contributing to child partition
    integer :: n_active    !< num elements with active children
    integer :: n_frozen    !< num elements with frozen children
    integer, allocatable :: id_elem(:) !< IDs of contributing elements
  end type MeshMapToChild_3D

  !=============================================================================

end module Mesh_Map_To_Child__3D
