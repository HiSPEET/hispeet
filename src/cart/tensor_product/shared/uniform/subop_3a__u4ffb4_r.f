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
!> Performs `v += g AᵀxIxI u` for single element
!>
!>   * k: unrolled and jammed with length of 4 (u4)
!>   * j: joined with i and flattened (f)
!>   * i: joined with j and flattened (f)
!>   * p: blocked with length of 4 (b4)
!>   * explicit remainder handling (r)

subroutine SubOp_3a(g, A, u, v)
  !$acc routine vector
  real(RNP), intent(in)    :: g                   !< metric factor
  real(RNP), intent(in)    :: A(__NA__,__NA__)    !< A
  real(RNP), intent(in)    :: u(__NA__**2,__NA__) !< operand
  real(RNP), intent(inout) :: v(__NA__**2,__NA__) !< result

  real(RNP) :: tmp0, tmp1, tmp2, tmp3
  integer   :: ij, k, p, q

#if __NA_T4__ > 0

  !$acc loop independent vector
  do k = 1, __NA_T4__, 4
    !$acc loop seq
    do q = 1, __NA_T4__, 4
      !$acc loop independent vector
      do ij = 1, __NA__**2
        tmp0 = 0
        tmp1 = 0
        tmp2 = 0
        tmp3 = 0
        do p = q, q+3
          tmp0 = tmp0 + A(p,k  ) * u(ij,p)
          tmp1 = tmp1 + A(p,k+1) * u(ij,p)
          tmp2 = tmp2 + A(p,k+2) * u(ij,p)
          tmp3 = tmp3 + A(p,k+3) * u(ij,p)
        end do
        v(ij,k  ) = v(ij,k  ) + g * tmp0
        v(ij,k+1) = v(ij,k+1) + g * tmp1
        v(ij,k+2) = v(ij,k+2) + g * tmp2
        v(ij,k+3) = v(ij,k+3) + g * tmp3
      end do
    end do
  end do

#endif

#if __NA__ == __NA_T4__ + 1

  ! remainder: k = 1:na-1, p = na ............................................

  p = __NA__

  !$acc loop collapse(2) vector
  do k = 1, __NA_T4__, 4
    do ij = 1, __NA__**2
      v(ij,k  ) = v(ij,k  ) + g * A(p,k  ) * u(ij,p)
      v(ij,k+1) = v(ij,k+1) + g * A(p,k+1) * u(ij,p)
      v(ij,k+2) = v(ij,k+2) + g * A(p,k+2) * u(ij,p)
      v(ij,k+3) = v(ij,k+3) + g * A(p,k+3) * u(ij,p)
    end do
  end do

  ! remainder: k = na, p = 1:na ..............................................

  !$acc loop vector
  do ij = 1, __NA__**2
    tmp0 = 0
    do p = 1, __NA__
      tmp0 = tmp0 + A(p,__NA__) * u(ij,p)
    end do
    v(ij,__NA__) = v(ij,__NA__) + g * tmp0
  end do

#elif __NA__ == __NA_T4__ + 2

  ! remainder: k = 1:na-2, p = na-1:na .......................................

  !$acc loop collapse(2) vector
  do k = 1, __NA_T4__, 4
    do ij = 1, __NA__**2
      do p = __NA__-1, __NA__
        v(ij,k  ) = v(ij,k  ) + g * A(p,k  ) * u(ij,p)
        v(ij,k+1) = v(ij,k+1) + g * A(p,k+1) * u(ij,p)
        v(ij,k+2) = v(ij,k+2) + g * A(p,k+2) * u(ij,p)
        v(ij,k+3) = v(ij,k+3) + g * A(p,k+3) * u(ij,p)
      end do
    end do
  end do

  ! remainder: k = na-1:na, p = 1:na .........................................

  !$acc loop vector
  do ij = 1, __NA__**2
    tmp0 = 0
    tmp1 = 0
    do p = 1, __NA__
      tmp0 = tmp0 + A(p,__NA__-1) * u(ij,p)
      tmp1 = tmp1 + A(p,__NA__  ) * u(ij,p)
    end do
    v(ij,__NA__-1) = v(ij,__NA__-1) + g * tmp0
    v(ij,__NA__  ) = v(ij,__NA__  ) + g * tmp1
  end do

#elif __NA__ == __NA_T4__ + 3

  ! remainder: k = 1:na-3, p = na-2:na .......................................

  !$acc loop collapse(2) vector
  do k = 1, __NA_T4__, 4
    do ij = 1, __NA__**2
      do p = __NA__-2, __NA__
        v(ij,k  ) = v(ij,k  ) + g * A(p,k  ) * u(ij,p)
        v(ij,k+1) = v(ij,k+1) + g * A(p,k+1) * u(ij,p)
        v(ij,k+2) = v(ij,k+2) + g * A(p,k+2) * u(ij,p)
        v(ij,k+3) = v(ij,k+3) + g * A(p,k+3) * u(ij,p)
      end do
    end do
  end do

  ! remainder: j = na-2:na, p = 1:na .........................................

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
    v(ij,__NA__-2) = v(ij,__NA__-2) + g * tmp0
    v(ij,__NA__-1) = v(ij,__NA__-1) + g * tmp1
    v(ij,__NA__  ) = v(ij,__NA__  ) + g * tmp2
  end do

#endif

end subroutine SubOp_3a
