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
!>   * k: simple loop (l)
!>   * j: joined with i and flattened (f)
!>   * i: joined with j and flattened (f)
!>   * p: simple loop (l)

subroutine SubOp_3i(g, A, u, v)
  !$acc routine vector
  real(RNP), intent(in)  :: g                   !< metric factor
  real(RNP), intent(in)  :: A(__NA__,__NA__)    !< A
  real(RNP), intent(in)  :: u(__NA__**2,__NA__) !< operand
  real(RNP), intent(out) :: v(__NA__**2,__NA__) !< result

  real(RNP) :: tmp
  integer   :: ij, k, p

  !$acc loop collapse(2) vector
  do k = 1, __NA__
  !DIR$ SIMD VECREMAINDER
  do ij = 1, __NA__**2
    tmp = 0
    do p = 1, __NA__
      tmp = tmp + A(p,k) * u(ij,p)
    end do
    v(ij,k) = g * tmp
  end do
  end do

end subroutine SubOp_3i
