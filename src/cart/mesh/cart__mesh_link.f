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

!> summary:  Cartesian mesh link
!> author:   Joerg Stiller
!> date:     2017/04/24
!>
!>### Cartesian mesh link
!===============================================================================

module CART__Mesh_Link
  use Kind_Parameters, only: IXS
  implicit none
  private

  public :: MeshLink
  public :: ElementLink

  !-----------------------------------------------------------------------------
  !> Structure for keeping element linking information

  type ElementLink
    integer      :: id = 0        !< local (virtual) element ID
    integer(IXS) :: face(6)   = 0 !< set to 1(0) for (un)linked faces
    integer(IXS) :: edge(12)  = 0 !< set to 1(0) for (un)linked edges
    integer(IXS) :: vertex(8) = 0 !< set to 1(0) for (un)linked vertices
  end type ElementLink

  !-----------------------------------------------------------------------------
  !> Structure for keeping mesh coupling information
  !>
  !> The structure provides the following lists
  !>
  !>    *  in `face` the IDs of faces linked with partition `part`
  !>    *  if `part` is the local partition, in `coupled_face` the
  !>       IDs of the coupled local faces
  !>    *  in `master` the ElementLink data of local elements possessing
  !>       a ghost in partition `part`
  !>    *  in `ghost` the ElementLink data of ghost elements governed by
  !>       a master in partition `part`

  type MeshLink
    integer :: part = -1                        !< remote partition ID
    integer :: nf = 0                           !< number of linked faces
    integer :: nm = 0                           !< number of master elements
    integer :: ng = 0                           !< number of ghost elements
    integer, allocatable :: face(:)             !< list of linked faces
    integer, allocatable :: coupled_face(:)     !< list of coupled local faces
    type(ElementLink), allocatable :: master(:) !< master elements
    type(ElementLink), allocatable :: ghost(:)  !< ghost elements
  end type MeshLink

!===============================================================================

end module CART__Mesh_Link
