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
!> Performs `v = g D AᵀxIxI u` for single element
!>
!>   * k: unrolled and jammed with length of 4 (u4)
!>   * j: joined with i and flattened (f)
!>   * i: joined with j and flattened (f)
!>   * p: simple loop (l)
!>   * explicit remainder handling (r)

subroutine SubOp_3id(g, A, D, u, v)
  !$acc routine vector
  real(RNP), intent(in)  :: g                   !< metric factor
  real(RNP), intent(in)  :: A(__NA__,__NA__)    !< A
  real(RNP), intent(in)  :: D(__NA__**2,__NA__) !< diagonal operator
  real(RNP), intent(in)  :: u(__NA__**2,__NA__) !< operand
  real(RNP), intent(out) :: v(__NA__**2,__NA__) !< result

  real(RNP) :: tmp0, tmp1, tmp2, tmp3
  integer   :: ij, k, p

#if __NA_T4__ > 0

  !$acc loop collapse(2) vector
  do k = 1, __NA_T4__, 4
  do ij = 1, __NA__**2
    tmp0 = 0
    tmp1 = 0
    tmp2 = 0
    tmp3 = 0
    do p = 1, __NA__
      tmp0 = tmp0 + A(p,k  ) * u(ij,p)
      tmp1 = tmp1 + A(p,k+1) * u(ij,p)
      tmp2 = tmp2 + A(p,k+2) * u(ij,p)
      tmp3 = tmp3 + A(p,k+3) * u(ij,p)
    end do
    v(ij,k  ) = g * D(ij,k  ) * tmp0
    v(ij,k+1) = g * D(ij,k+1) * tmp1
    v(ij,k+2) = g * D(ij,k+2) * tmp2
    v(ij,k+3) = g * D(ij,k+3) * tmp3
  end do
  end do

#endif

#if __NA__ == __NA_T4__ + 1

  ! remainder: k = na, p = 1:na ..............................................

  !$acc loop vector
  do ij = 1, __NA__**2
    tmp0 = 0
    do p = 1, __NA__
      tmp0 = tmp0 + A(p,__NA__) * u(ij,p)
    end do
    v(ij,__NA__) = g * D(ij,__NA__) * tmp0
  end do

#elif __NA__ == __NA_T4__ + 2

  ! remainder: k = na-1:na, p = 1:na .........................................

  !$acc loop vector
  do ij = 1, __NA__**2
    tmp0 = 0
    tmp1 = 0
    do p = 1, __NA__
      tmp0 = tmp0 + A(p,__NA__-1) * u(ij,p)
      tmp1 = tmp1 + A(p,__NA__  ) * u(ij,p)
    end do
    v(ij,__NA__-1) = g * D(ij,__NA__-1) * tmp0
    v(ij,__NA__  ) = g * D(ij,__NA__  ) * tmp1
  end do

#elif __NA__ == __NA_T4__ + 3

  !$acc loop vector
  do ij = 1, __NA__**2
    tmp0 = 0
    tmp1 = 0
    tmp2 = 0
    do p = 1, __NA__
      tmp0 = tmp0 + A(p,__NA__-2) * u(ij,p)
      tmp1 = tmp1 + A(p,__NA__-1) * u(ij,p)
      tmp2 = tmp2 + A(p,__NA__  ) * u(ij,p)
    end do
    v(ij,__NA__-2) = g * D(ij,__NA__-2) * tmp0
    v(ij,__NA__-1) = g * D(ij,__NA__-1) * tmp1
    v(ij,__NA__  ) = g * D(ij,__NA__  ) * tmp2
  end do

#endif

end subroutine SubOp_3id
