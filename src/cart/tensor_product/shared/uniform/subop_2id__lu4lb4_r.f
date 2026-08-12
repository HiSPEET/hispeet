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
!> Performs `v = g D IxAᵀxI u`
!>
!>   * k: simple loop (l)
!>   * j: unrolled and jammed with length of 4 (u4)
!>   * i: simple loop (l)
!>   * p: blocked with length of 4 (b4)
!>   * explicit remainder handling (r)

subroutine SubOp_2id(g, A, D, u, v)
  !$acc routine vector
  real(RNP), intent(in)  :: g                        !< metric factor
  real(RNP), intent(in)  :: A(__NA__,__NA__)         !< A
  real(RNP), intent(in)  :: D(__NA__,__NA__,__NA__)  !< diagonal operator
  real(RNP), intent(in)  :: u(__NA__,__NA__,__NA__)  !< operand
  real(RNP), intent(out) :: v(__NA__,__NA__,__NA__)  !< result

  real(RNP) :: tmp0, tmp1, tmp2, tmp3
  integer   :: i, j, k, p, q

  !$acc loop collapse(3) vector
  do k = 1, __NA__
  do j = 1, __NA__
  do i = 1, __NA__
    v(i,j,k) = 0
  end do
  end do
  end do

#if __NA_T4__ > 0

  !$acc loop collapse(2) independent vector
  do k = 1, __NA__
  do j = 1, __NA_T4__, 4
    !$acc loop seq
    do q = 1, __NA_T4__, 4
      !$acc loop independent vector
      do i = 1, __NA__
        tmp0 = 0
        tmp1 = 0
        tmp2 = 0
        tmp3 = 0
        do p = q, q+3
          tmp0 = tmp0 + A(p,j  ) * u(i,p,k)
          tmp1 = tmp1 + A(p,j+1) * u(i,p,k)
          tmp2 = tmp2 + A(p,j+2) * u(i,p,k)
          tmp3 = tmp3 + A(p,j+3) * u(i,p,k)
        end do
        v(i,j  ,k) = v(i,j  ,k) + g * tmp0
        v(i,j+1,k) = v(i,j+1,k) + g * tmp1
        v(i,j+2,k) = v(i,j+2,k) + g * tmp2
        v(i,j+3,k) = v(i,j+3,k) + g * tmp3
      end do
    end do
  end do
  end do

#endif

#if __NA__ == __NA_T4__ + 1

  !$acc loop independent vector
  do k = 1, __NA__

    ! remainder: j = 1: na-1, p = na .........................................

    p = __NA__

    !$acc loop collapse(2) vector
    do j = 1, __NA_T4__, 4
    do i = 1, __NA__
      v(i,j  ,k) = v(i,j  ,k) + g * A(p,j  ) * u(i,p,k)
      v(i,j+1,k) = v(i,j+1,k) + g * A(p,j+1) * u(i,p,k)
      v(i,j+2,k) = v(i,j+2,k) + g * A(p,j+2) * u(i,p,k)
      v(i,j+3,k) = v(i,j+3,k) + g * A(p,j+3) * u(i,p,k)
    end do
    end do

    ! remainder: j = na, p = 1:na  ...........................................

    !$acc loop vector
    do i = 1, __NA__
      tmp0 = 0
      do p = 1, __NA__
        tmp0 = tmp0 + A(p,__NA__) * u(i,p,k)
      end do
      v(i,__NA__,k) = v(i,__NA__,k) + g * tmp0
    end do

  end do

#elif __NA__ == __NA_T4__ + 2

  !$acc loop independent vector
  do k = 1, __NA__

    ! remainder: j = 1:na-2, p = na-1:na .....................................

    !$acc loop collapse(2) vector
    do j = 1, __NA_T4__, 4
    do i = 1, __NA__
      do p = __NA__-1, __NA__
        v(i,j  ,k) = v(i,j  ,k) + g * A(p,j  ) * u(i,p,k)
        v(i,j+1,k) = v(i,j+1,k) + g * A(p,j+1) * u(i,p,k)
        v(i,j+2,k) = v(i,j+2,k) + g * A(p,j+2) * u(i,p,k)
        v(i,j+3,k) = v(i,j+3,k) + g * A(p,j+3) * u(i,p,k)
      end do
    end do
    end do

    ! remainder: j = na-1:na, p = 1:na .......................................

    !$acc loop vector
    do i = 1, __NA__
      tmp0 = 0
      tmp1 = 0
      do p = 1, __NA__
        tmp0 = tmp0 + A(p,__NA__-1) * u(i,p,k)
        tmp1 = tmp1 + A(p,__NA__  ) * u(i,p,k)
      end do
      v(i,__NA__-1,k) = v(i,__NA__-1,k) + g * tmp0
      v(i,__NA__  ,k) = v(i,__NA__  ,k) + g * tmp1
    end do

  end do

#elif __NA__ == __NA_T4__ + 3

  !$acc loop independent vector
  do k = 1, __NA__

    ! remainder: j = 1:na-3, p = na-2:na .....................................

    !$acc loop collapse(2) vector
    do j = 1, __NA_T4__, 4
    do i = 1, __NA__
      do p = __NA__-2, __NA__
        v(i,j  ,k) = v(i,j  ,k) + g * A(p,j  ) * u(i,p,k)
        v(i,j+1,k) = v(i,j+1,k) + g * A(p,j+1) * u(i,p,k)
        v(i,j+2,k) = v(i,j+2,k) + g * A(p,j+2) * u(i,p,k)
        v(i,j+3,k) = v(i,j+3,k) + g * A(p,j+3) * u(i,p,k)
      end do
    end do
    end do

    ! remainder: j = na-2:na, p = 1:na .......................................

    !$acc loop vector
    do i = 1, __NA__
      tmp0 = 0
      tmp1 = 0
      tmp2 = 0
      do p = 1, __NA__
        tmp0 = tmp0 + A(p,__NA__-2) * u(i,p,k)
        tmp1 = tmp1 + A(p,__NA__-1) * u(i,p,k)
        tmp2 = tmp2 + A(p,__NA__  ) * u(i,p,k)
      end do
      v(i,__NA__-2,k) = v(i,__NA__-2,k) + g * tmp0
      v(i,__NA__-1,k) = v(i,__NA__-1,k) + g * tmp1
      v(i,__NA__  ,k) = v(i,__NA__  ,k) + g * tmp2
    end do

  end do

#endif

  ! v = D v ....................................................................

  !$acc loop collapse(3) vector
  do k = 1, __NA__
  do j = 1, __NA__
  do i = 1, __NA__
    v(i,j,k) = D(i,j,k) * v(i,j,k)
  end do
  end do
  end do

end subroutine SubOp_2id
