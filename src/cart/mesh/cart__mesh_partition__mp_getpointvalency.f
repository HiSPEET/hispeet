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

!> summary:  Evaluation of the valency of mesh points
!> author:   Joerg Stiller
!> date:     2018/10/23
!>
!>### Evaluation of the valency of mesh points
!===============================================================================

submodule(CART__Mesh_Partition) MP_GetPointValency
  implicit none

contains

!-------------------------------------------------------------------------------
!> Computes the valency of mesh points
!>
!> The valency of a point equals the number of mesh elements sharing this point.
!> In the Gauss case, every point belongs to only one element and, hence, has a
!> valency of 1. In the Lobatto case, the points are shared across faces, edges
!> and vertices. As a result, the valency of mesh points can range between 1 and
!> 8.

module subroutine GetPointValency(mesh, nodes, v)
  class(MeshPartition), intent(in)  :: mesh          !< mesh parition
  character,            intent(in)  :: nodes         !< 'G' or 'L'
  integer,              intent(out) :: v(0:,0:,0:,:) !< point valency

  integer :: e, i, j, k, P

  ! initialization .............................................................

  P = ubound(v, 1)

  ! default valency ............................................................

  do e = 1, mesh%ne
    do k = 0, P
    do j = 0, P
    do i = 0, P
      v(i,j,k,e) = 1
    end do
    end do
    end do
  end do

  ! Lobatto point valency ......................................................

  if (nodes == 'L') then

    do e = 1, mesh%ne

      associate(face => mesh % element(e)%face)
        if (face(1) % neighbor > 0)  v( 0, :, :, e)  =  v( 0, :, :, e)  +  1
        if (face(2) % neighbor > 0)  v( P, :, :, e)  =  v( P, :, :, e)  +  1
        if (face(3) % neighbor > 0)  v( :, 0, :, e)  =  v( :, 0, :, e)  +  1
        if (face(4) % neighbor > 0)  v( :, P, :, e)  =  v( :, P, :, e)  +  1
        if (face(5) % neighbor > 0)  v( :, :, 0, e)  =  v( :, :, 0, e)  +  1
        if (face(6) % neighbor > 0)  v( :, :, P, e)  =  v( :, :, P, e)  +  1
      end associate

      associate(edge => mesh % element(e)%edge)
        if (edge( 1) % neighbor > 0)  v( :, 0, 0, e)  =  v( :, 0, 0, e)  +  1
        if (edge( 2) % neighbor > 0)  v( :, P, 0, e)  =  v( :, P, 0, e)  +  1
        if (edge( 3) % neighbor > 0)  v( :, 0, P, e)  =  v( :, 0, P, e)  +  1
        if (edge( 4) % neighbor > 0)  v( :, P, P, e)  =  v( :, P, P, e)  +  1
        if (edge( 5) % neighbor > 0)  v( 0, :, 0, e)  =  v( 0, :, 0, e)  +  1
        if (edge( 6) % neighbor > 0)  v( P, :, 0, e)  =  v( P, :, 0, e)  +  1
        if (edge( 7) % neighbor > 0)  v( 0, :, P, e)  =  v( 0, :, P, e)  +  1
        if (edge( 8) % neighbor > 0)  v( P, :, P, e)  =  v( P, :, P, e)  +  1
        if (edge( 9) % neighbor > 0)  v( 0, 0, :, e)  =  v( 0, 0, :, e)  +  1
        if (edge(10) % neighbor > 0)  v( P, 0, :, e)  =  v( P, 0, :, e)  +  1
        if (edge(11) % neighbor > 0)  v( 0, P, :, e)  =  v( 0, P, :, e)  +  1
        if (edge(12) % neighbor > 0)  v( P, P, :, e)  =  v( P, P, :, e)  +  1
      end associate

      associate(vertex => mesh % element(e)%vertex)
        if (vertex(1) % neighbor > 0)  v( 0, 0, 0, e)  =  v( 0, 0, 0, e)  +  1
        if (vertex(2) % neighbor > 0)  v( P, 0, 0, e)  =  v( P, 0, 0, e)  +  1
        if (vertex(3) % neighbor > 0)  v( 0, P, 0, e)  =  v( 0, P, 0, e)  +  1
        if (vertex(4) % neighbor > 0)  v( P, P, 0, e)  =  v( P, P, 0, e)  +  1
        if (vertex(5) % neighbor > 0)  v( 0, 0, P, e)  =  v( 0, 0, P, e)  +  1
        if (vertex(6) % neighbor > 0)  v( P, 0, P, e)  =  v( P, 0, P, e)  +  1
        if (vertex(7) % neighbor > 0)  v( 0, P, P, e)  =  v( 0, P, P, e)  +  1
        if (vertex(8) % neighbor > 0)  v( P, P, P, e)  =  v( P, P, P, e)  +  1
      end associate

    end do

  end if

end subroutine GetPointValency

!===============================================================================

end submodule MP_GetPointValency
