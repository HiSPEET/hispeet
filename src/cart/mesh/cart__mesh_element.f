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

!> summary:  Cartesian mesh element
!> author:   Joerg Stiller
!> date:     2017/04/03
!>
!>### Cartesian mesh boundary
!===============================================================================

module CART__Mesh_Element
  use Kind_Parameters, only: RNP
  implicit none
  private

  public :: MeshElement

  !-----------------------------------------------------------------------------
  !> Structure for keeping element face data

  type ElementFace
    integer :: id       = 0  !< mesh face ID
    integer :: neighbor = 0  !< neighbor element, including virtual one
    integer :: boundary = 0  !< adjacent boundary, if any
  end type ElementFace

  !-----------------------------------------------------------------------------
  !> Structure for keeping element edge data

  type ElementEdge
    integer :: neighbor = 0  !< neighbor element, including virtual one
  end type ElementEdge

  !-----------------------------------------------------------------------------
  !> Structure for keeping element vertex data

  type ElementVertex
    real(RNP) :: x(3)          !< coordinates
    integer   :: neighbor = 0  !< neighbor element, including virtual one
  end type ElementVertex

  !-----------------------------------------------------------------------------
  !> Cartesian mesh element

  type MeshElement
    integer             :: global_id = 0 !< global ID
    type(ElementFace)   :: face(6)       !< face ID, neighbor and boundary info
    type(ElementEdge)   :: edge(12)      !< edge neighbor
    type(ElementVertex) :: vertex(8)     !< vertex coordinates and neighbor
  end type MeshElement

!===============================================================================

end module CART__Mesh_Element
