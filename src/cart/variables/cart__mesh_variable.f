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

!> summary:  Multi-component variable defined on mesh partition
!> author:   Joerg Stiller
!> date:     2016/03/18
!>
!>### Multi-component variable defined on mesh partition
!>
!> @todo
!> This is just a sketch
!> @endtodo
!===============================================================================

module Mesh_Variable

  use Kind_Parameters, only: RNP
  use CART__Mesh_Partition

  implicit none
  private

  public :: MeshVariable

  !-----------------------------------------------------------------------------
  !> Type for handling mesh variables

  type MeshVariable

    type(MeshPartition),   pointer :: mesh => null() !< associated partition
    real(RNP), contiguous, pointer :: val(:,:,:,:,:) => null() !< values
    character(:),          pointer :: name(:) => null()! component names

    logical :: is_original = .false. !< indicates original mesh variable
    logical :: has_ghosts  = .false. !< indicates whether ghosts are included

  contains

  ! New
  ! Set
  ! Merge
  ! GetTrace:  extract trace variable
  ! GetHandle: return mesh variable linked to a section of the given one
  ! GetClone?  get full or partial copy
  ! Access:    get direct access to single or section of components

  end type MeshVariable

!===============================================================================

end module Mesh_Variable
