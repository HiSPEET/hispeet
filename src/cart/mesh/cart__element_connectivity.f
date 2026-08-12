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

!> summary:  Structures for temporary storage of element connectivity
!> author:   Joerg Stiller
!> date:     2017/04/03
!>
!>### Structures for temporary storage of element connectivity
!===============================================================================

module CART__Element_Connectivity
  private

  public :: ElementConnectivity

  !-----------------------------------------------------------------------------
  !> Structure for identifying a neighbor element

  type NeighborElement
    integer :: part = -1  !< host partition
    integer :: id   = -1  !< element ID within host partition
  end type NeighborElement

  !-----------------------------------------------------------------------------
  !> Structure for storing the connectivity of a Cartesian mesh element

  type ElementConnectivity
    type(NeighborElement) :: face(6)    !< face neighbors
    type(NeighborElement) :: edge(12)   !< edge neighbors
    type(NeighborElement) :: vertex(8)  !< vertex neighbors
  end type ElementConnectivity

!===============================================================================

end module CART__Element_Connectivity
