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

!-------------------------------------------------------------------------------
!> Performs `v = g AᵀxIxI u` for single element
!>
!>   * k: unrolled and jammed with length of 4 (u4)
!>   * j: joined with i and flattened (f)
!>   * i: joined with j and flattened (f)
!>   * p: simple loop (l)
!>   * explicit remainder handling (r)

subroutine SubOp_3i(g, A, u, v)
  !$acc routine vector
  real(RNP), intent(in)  :: g                   !< metric factor
  real(RNP), intent(in)  :: A(__NA__,__NA__)    !< A
  real(RNP), intent(in)  :: u(__NA__**2,__NA__) !< operand
  real(RNP), intent(out) :: v(__NA__**2,__NA__) !< result

  real(RNP) :: tmp0, tmp1
  integer   :: ij, k, p

  !$acc loop collapse(2) vector
  do k = 1, __NA_T2__, 2
  do ij = 1, __NA__**2
    tmp0 = 0
    tmp1 = 0
    do p = 1, __NA__
      tmp0 = tmp0 + A(p,k  ) * u(ij,p)
      tmp1 = tmp1 + A(p,k+1) * u(ij,p)
    end do
    v(ij,k  ) = g * tmp0
    v(ij,k+1) = g * tmp1
  end do
  end do

#if __NA__ == __NA_T2__ + 1

  ! remainder: k = na ........................................................

  !$acc loop vector
  do ij = 1, __NA__**2
    tmp0 = 0
    do p = 1, __NA__
      tmp0 = tmp0 + A(p,__NA__) * u(ij,p)
    end do
    v(ij,__NA__) = g * tmp0
  end do

#endif

end subroutine SubOp_3i
