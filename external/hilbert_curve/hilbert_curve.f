! This file is part of FD4: Four-Dimensional Distributed Dynamic Data structures
!
! Copyright 2010-2016 Matthias Lieber, ZIH, TU Dresden, Germany
! matthias.lieber@tu-dresden.de
!
! FD4 is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! FD4 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with FD4.  If not, see <http://www.gnu.org/licenses/>.

!--------------------------------------------------------------!
!                                                              !
!                       3D Hilbert Curve                       !
!                                                              !
!--------------------------------------------------------------!


!! Simple 3D Hilbert curve subroutines.
!! Based on a table description of the curve. Parts of the table came from http://wiki.tcl.tk/8723.
!!
!! @author Matthias Lieber
!!
!! 2026/01/15 adapted to HiSPEET by Joerg Stiller

module Hilbert_Curve
  implicit none
  private

  public :: Hilbert_I2C, Hilbert_C2I

  !! Integer kind parameter for Hilbert Curve indices and related quantities
  !!
  !! Range of 10^18 is chosen to cope with large curve indices
  integer, parameter, public :: IHC = selected_int_kind(18)

  !! Type for one entry in the 3D Hilbert curve description table.
  !!
  !! Values of the $mv item are:
  !! - 1/-1 means x+1/x-1
  !! - 2/-2 means y+1/y-1
  !! - 3/-3 means z+1/z-1
  type hilbert_table
    integer :: mv(0:6)      !! spatial moves of the curve
    integer :: next(0:7)    !! next states when refinining (octant is sequence number)
    integer :: num2ref(0:7) !! sequence number on curve primitve to reference (spatial) octant
    integer :: ref2num(0:7) !! reference (spatial) octant to sequence number on curve primitve
  end type hilbert_table


  ! (parts of the) table from http://wiki.tcl.tk/8723
  type(hilbert_table), parameter :: tab(0:23) =                                                               &
  [                                                                                                           &
    hilbert_table( [ 1, 3,-1, 2, 1,-3,-1], [ 1, 2, 2,20,20, 4, 4, 5], [0,1,5,4,6,7,3,2], [0,1,7,6,3,2,4,5] ), &
    hilbert_table( [ 3, 2,-3, 1, 3,-2,-3], [ 6, 7, 7,19,19, 9, 9,10], [0,4,6,2,3,7,5,1], [0,7,3,4,1,6,2,5] ), &
    hilbert_table( [ 1, 2,-1, 3, 1,-2,-1], [11, 0, 0,23,23,14,14,15], [0,1,3,2,6,7,5,4], [0,1,3,2,7,6,4,5] ), &
    hilbert_table( [-1,-3, 1, 2,-1, 3, 1], [16,17,17, 7, 7,13,13, 8], [5,4,0,1,3,2,6,7], [2,3,5,4,1,0,6,7] ), &
    hilbert_table( [ 1,-2,-1,-3, 1, 2,-1], [19,14,14,10,10, 0, 0,22], [6,7,5,4,0,1,3,2], [4,5,7,6,3,2,0,1] ), &
    hilbert_table( [ 3,-2,-3,-1, 3, 2,-3], [23, 9, 9,15,15, 7, 7,21], [3,7,5,1,0,4,6,2], [4,3,7,0,5,2,6,1] ), &
    hilbert_table( [ 2, 1,-2, 3, 2,-1,-2], [ 0,11,11,13,13,15,15,14], [0,2,3,1,5,7,6,4], [0,3,1,2,7,4,6,5] ), &
    hilbert_table( [ 3, 1,-3, 2, 3,-1,-3], [ 2, 1, 1, 3, 3, 5, 5, 4], [0,4,5,1,3,7,6,2], [0,3,7,4,1,2,6,5] ), &
    hilbert_table( [-3,-2, 3, 1,-3, 2, 3], [21,18,18,11,11,20,20,23], [6,2,0,4,5,1,3,7], [2,5,1,6,3,4,0,7] ), &
    hilbert_table( [ 3,-1,-3,-2, 3, 1,-3], [13, 5, 5,14,14, 1, 1,17], [3,7,6,2,0,4,5,1], [4,7,3,0,5,6,2,1] ), &
    hilbert_table( [ 2,-1,-2,-3, 2, 1,-2], [ 3,15,15, 4, 4,11,11,12], [5,7,6,4,0,2,3,1], [4,7,5,6,3,0,2,1] ), &
    hilbert_table( [ 2, 3,-2, 1, 2,-3,-2], [ 7, 6, 6, 8, 8,10,10, 9], [0,2,6,4,5,7,3,1], [0,7,1,6,3,4,2,5] ), &
    hilbert_table( [-1, 3, 1,-2,-1,-3, 1], [ 5,13,13,18,18,17,17, 1], [3,2,6,7,5,4,0,1], [6,7,1,0,5,4,2,3] ), &
    hilbert_table( [-1,-2, 1, 3,-1, 2, 1], [22,12,12, 6, 6, 3, 3,19], [3,2,0,1,5,4,6,7], [2,3,1,0,5,4,6,7] ), &
    hilbert_table( [ 1,-3,-1,-2, 1, 3,-1], [ 8, 4, 4, 9, 9, 2, 2,16], [6,7,3,2,0,1,5,4], [4,5,3,2,7,6,0,1] ), &
    hilbert_table( [ 2,-3,-2,-1, 2, 3,-2], [20,10,10, 5, 5, 6, 6,18], [5,7,3,1,0,2,6,4], [4,3,5,2,7,0,6,1] ), &
    hilbert_table( [-3, 2, 3,-1,-3,-2, 3], [10,20,20,22,22,18,18, 6], [5,1,3,7,6,2,0,4], [6,1,5,2,7,0,4,3] ), &
    hilbert_table( [-1, 2, 1,-3,-1,-2, 1], [15, 3, 3,21,21,12,12,11], [5,4,6,7,3,2,0,1], [6,7,5,4,1,0,2,3] ), &
    hilbert_table( [-3, 1, 3,-2,-3,-1, 3], [ 4, 8, 8,12,12,16,16, 2], [6,2,3,7,5,1,0,4], [6,5,1,2,7,4,0,3] ), &
    hilbert_table( [-2,-3, 2, 1,-2, 3, 2], [18,21,21, 1, 1,23,23,20], [6,4,0,2,3,1,5,7], [2,5,3,4,1,6,0,7] ), &
    hilbert_table( [-3,-1, 3, 2,-3, 1, 3], [17,16,16, 0, 0, 8, 8,13], [5,1,0,4,6,2,3,7], [2,1,5,6,3,0,4,7] ), &
    hilbert_table( [-2, 1, 2,-3,-2,-1, 2], [14,19,19,17,17,22,22, 0], [6,4,5,7,3,1,0,2], [6,5,7,4,1,2,0,3] ), &
    hilbert_table( [-2, 3, 2,-1,-2,-3, 2], [ 9,23,23,16,16,21,21, 7], [3,1,5,7,6,4,0,2], [6,1,7,0,5,2,4,3] ), &
    hilbert_table( [-2,-1, 2, 3,-2, 1, 2], [12,22,22, 2, 2,19,19, 3], [3,1,0,2,6,4,5,7], [2,1,3,0,5,6,4,7] )  &
   ]

contains


  !! Get index on curve at given 3D coordinate (coordinates-to-index).
  !!
  !! Has logarithmic complexity: number of loop iterations = $lev.
  !! Due to performance reasons, no checking of arguments is performend!
  pure subroutine Hilbert_C2I(lev, p, idx)

    implicit none

    integer,      intent(in)  :: lev  !! level of the hilbert curve (1=2x2x2, 2=4x4x4, 3=8x8x8, etc)
    integer(IHC), intent(in)  :: p(3) !! 3D coordinates [0 ... 2^$lev - 1]
    integer(IHC), intent(out) :: idx  !! index on the 3D hilbert curve [0 ... (2^$lev)^3 - 1 ]

    integer(IHC) :: len, dim, cp(3)
    integer      :: l, q, s, rq

    cp(:) = p
    idx   = 0
    s     = 0
    dim   = 2**lev

    ! length of curve in base octant
    len = dim**3

    do l = lev,1,-1

      ! dimension and length of refined octants
      dim = dim/2
      len = len/8

      ! determine spatial octant and reduce spatial coordinates
      ! so that they are relative to refined octant
      ! reference octant:
      ! z=0 | z=1
      ! 2 3 | 6 7
      ! 0 1 | 4 5
      rq = 0
      if(cp(1) >= dim) then
        rq    = ior(1,rq)
        cp(1) = cp(1)-dim
      end if
      if(cp(2) >= dim) then
        rq    = ior(2,rq)
        cp(2) = cp(2)-dim
      end if
      if(cp(3) >= dim) then
        rq    = ior(4,rq)
        cp(3) = cp(3)-dim
      end if

      ! get sequence number of octant q in the curve primitive
      q = tab(s)%ref2num(rq)

      ! add offset to index
      idx = idx + q*len

      ! get next state
      s = tab(s)%next(q)

    end do

  end subroutine Hilbert_C2I


  !! Get 3D coordinate at given index on curve (index-to-coordinates).
  !!
  !! Has logarithmic complexity: number of loop iterations = $lev.
  !! Due to performance reasons, no checking of arguments is performend!
  pure subroutine Hilbert_I2C(lev, idx, p)

    implicit none

    integer,      intent(in)  :: lev   !! level of the hilbert curve (1=2x2x2, 2=4x4x4, 3=8x8x8, etc)
    integer(IHC), intent(in)  :: idx   !! index on the 3D hilbert curve [0 ... (2^$lev)^3 - 1 ]
    integer(IHC), intent(out) :: p(3)  !! 3D coordinates [0 ... 2^$lev - 1]

    integer(IHC) :: len, dim, cidx
    integer      :: l, q, s, rq

    p(:) = 0
    s    = 0
    cidx = idx
    dim  = 2**lev

    ! length of curve in base octant
    len = dim**3

    do l = lev,1,-1

      ! dimension of refined octants
      dim = dim/2

      ! determine sequence number of octant in the curve primitive.
      ! this alorithm is faster than doing simply this:
      !   len=len/8 ;  q=cidx/len
      ! or
      !   q=(8*cidx)/len
      if(cidx < len/2) then
        if(cidx < len/4) then
          if(cidx < len/8) then
            q = 0
          else
            q = 1
          end if
        else
          if(cidx < 3*len/8) then
            q = 2
          else
            q = 3
          end if
        end if
      else
        if(cidx < 3*len/4) then
          if(cidx < 5*len/8) then
            q = 4
          else
            q = 5
          end if
        else
          if(cidx < 7*len/8) then
            q = 6
          else
            q = 7
          end if
        end if
      end if

      ! get actual spatial octant, depending on curve primitive
      ! reference octant:
      ! z=0 | z=1
      ! 2 3 | 6 7
      ! 0 1 | 4 5
      rq = tab(s)%num2ref(q)

      ! update coordinates
      if(iand(1,rq) == 1) p(1) = p(1) + dim
      if(iand(2,rq) == 2) p(2) = p(2) + dim
      if(iand(4,rq) == 4) p(3) = p(3) + dim

      ! length of curve in refined octant
      len = len/8

      ! decrease idx by offset in current octant
      cidx = cidx - q*len

      ! next state
      s = tab(s)%next(q)

    end do

  end subroutine Hilbert_I2C

end module Hilbert_Curve
