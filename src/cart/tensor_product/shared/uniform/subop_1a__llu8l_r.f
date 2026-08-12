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
!> Performs `v += g IxIxAᵀ u`
!>
!>   * k: simple loop (l)
!>   * j: simple loop (l)
!>   * i: unrolled and jammed with length of 8 (u8)
!>   * p: simple loop (l)
!>   * explicit remainder handling (r)

subroutine SubOp_1a(g, A, u, v)
  !$acc routine vector
  real(RNP), intent(in)    :: g                        !< metric factor
  real(RNP), intent(in)    :: A(__NA__,__NA__)         !< A
  real(RNP), intent(in)    :: u(__NA__,__NA__,__NA__)  !< operand
  real(RNP), intent(inout) :: v(__NA__,__NA__,__NA__)  !< result

  real(RNP) :: tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7
  integer   :: i, j, k, p

#if __NA_T8__ > 0

  !$acc loop collapse(3) vector
  do k = 1, __NA__
  do j = 1, __NA__
  do i = 1, __NA_T8__, 8
    tmp0 = 0
    tmp1 = 0
    tmp2 = 0
    tmp3 = 0
    tmp4 = 0
    tmp5 = 0
    tmp6 = 0
    tmp7 = 0
    do p = 1, __NA__
      tmp0 = tmp0 + A(p,i  ) * u(p,j,k)
      tmp1 = tmp1 + A(p,i+1) * u(p,j,k)
      tmp2 = tmp2 + A(p,i+2) * u(p,j,k)
      tmp3 = tmp3 + A(p,i+3) * u(p,j,k)
      tmp4 = tmp4 + A(p,i+4) * u(p,j,k)
      tmp5 = tmp5 + A(p,i+5) * u(p,j,k)
      tmp6 = tmp6 + A(p,i+6) * u(p,j,k)
      tmp7 = tmp7 + A(p,i+7) * u(p,j,k)
    end do
    v(i  ,j,k) = v(i  ,j,k) + g * tmp0
    v(i+1,j,k) = v(i+1,j,k) + g * tmp1
    v(i+2,j,k) = v(i+2,j,k) + g * tmp2
    v(i+3,j,k) = v(i+3,j,k) + g * tmp3
    v(i+4,j,k) = v(i+4,j,k) + g * tmp4
    v(i+5,j,k) = v(i+5,j,k) + g * tmp5
    v(i+6,j,k) = v(i+6,j,k) + g * tmp6
    v(i+7,j,k) = v(i+7,j,k) + g * tmp7
  end do
  end do
  end do

#endif

#if __NA__ == __NA_T8__ + 1

  ! remainder: i = na ........................................................

  i = __NA_T8__ + 1

  !$acc loop collapse(2) vector
  do k = 1, __NA__
  do j = 1, __NA__
    tmp0 = 0
    do p = 1, __NA__
      tmp0 = tmp0 + A(p,i  ) * u(p,j,k)
    end do
    v(i,j,k) = v(i,j,k) + g * tmp0
  end do
  end do

#elif __NA__ == __NA_T8__ + 2

  ! remainder: i = na-1:na ...................................................

  i = __NA_T8__ + 1

  !$acc loop collapse(2) vector
  do k = 1, __NA__
  do j = 1, __NA__
    tmp0 = 0
    tmp1 = 0
    do p = 1, __NA__
      tmp0 = tmp0 + A(p,i  ) * u(p,j,k)
      tmp1 = tmp1 + A(p,i+1) * u(p,j,k)
    end do
    v(i  ,j,k) = v(i  ,j,k) + g * tmp0
    v(i+1,j,k) = v(i+1,j,k) + g * tmp1
  end do
  end do

#elif __NA__ == __NA_T8__ + 3

  ! remainder: i = na-2:na ...................................................

  i = __NA_T8__ + 1

  !$acc loop collapse(2) vector
  do k = 1, __NA__
  do i = 1, __NA__
    tmp0 = 0
    tmp1 = 0
    tmp2 = 0
    do p = 1, __NA__
      tmp0 = tmp0 + A(p,i  ) * u(p,j,k)
      tmp1 = tmp1 + A(p,i+1) * u(p,j,k)
      tmp2 = tmp2 + A(p,i+2) * u(p,j,k)
    end do
    v(i  ,j,k) = v(i  ,j,k) + g * tmp0
    v(i+1,j,k) = v(i+1,j,k) + g * tmp1
    v(i+2,j,k) = v(i+2,j,k) + g * tmp2
  end do
  end do

#elif __NA__ == __NA_T8__ + 4

  ! remainder: i = na-3:na ...................................................

  i = __NA_T8__ + 1

  !$acc loop collapse(2) vector
  do k = 1, __NA__
  do i = 1, __NA__
    tmp0 = 0
    tmp1 = 0
    tmp2 = 0
    tmp3 = 0
    do p = 1, __NA__
      tmp0 = tmp0 + A(p,i  ) * u(p,j,k)
      tmp1 = tmp1 + A(p,i+1) * u(p,j,k)
      tmp2 = tmp2 + A(p,i+2) * u(p,j,k)
      tmp3 = tmp3 + A(p,i+3) * u(p,j,k)
    end do
    v(i  ,j,k) = v(i  ,j,k) + g * tmp0
    v(i+1,j,k) = v(i+1,j,k) + g * tmp1
    v(i+2,j,k) = v(i+2,j,k) + g * tmp2
    v(i+3,j,k) = v(i+3,j,k) + g * tmp3
  end do
  end do

#elif __NA__ == __NA_T8__ + 5

  ! remainder: i = na-4:na ...................................................

  i = __NA_T8__ + 1

  !$acc loop collapse(2) vector
  do k = 1, __NA__
  do i = 1, __NA__
    tmp0 = 0
    tmp1 = 0
    tmp2 = 0
    tmp3 = 0
    tmp4 = 0
    do p = 1, __NA__
      tmp0 = tmp0 + A(p,i  ) * u(p,j,k)
      tmp1 = tmp1 + A(p,i+1) * u(p,j,k)
      tmp2 = tmp2 + A(p,i+2) * u(p,j,k)
      tmp3 = tmp3 + A(p,i+3) * u(p,j,k)
      tmp4 = tmp4 + A(p,i+4) * u(p,j,k)
    end do
    v(i  ,j,k) = v(i  ,j,k) + g * tmp0
    v(i+1,j,k) = v(i+1,j,k) + g * tmp1
    v(i+2,j,k) = v(i+2,j,k) + g * tmp2
    v(i+3,j,k) = v(i+3,j,k) + g * tmp3
    v(i+4,j,k) = v(i+4,j,k) + g * tmp4
  end do
  end do

#elif __NA__ == __NA_T8__ + 6

  ! remainder: i = na-5:na ...................................................

  i = __NA_T8__ + 1

  !$acc loop collapse(2) vector
  do k = 1, __NA__
  do i = 1, __NA__
    tmp0 = 0
    tmp1 = 0
    tmp2 = 0
    tmp3 = 0
    tmp4 = 0
    tmp5 = 0
    do p = 1, __NA__
      tmp0 = tmp0 + A(p,i  ) * u(p,j,k)
      tmp1 = tmp1 + A(p,i+1) * u(p,j,k)
      tmp2 = tmp2 + A(p,i+2) * u(p,j,k)
      tmp3 = tmp3 + A(p,i+3) * u(p,j,k)
      tmp4 = tmp4 + A(p,i+4) * u(p,j,k)
      tmp5 = tmp5 + A(p,i+5) * u(p,j,k)
    end doi
    v(i  ,j,k) = v(i  ,j,k) + g * tmp0
    v(i+1,j,k) = v(i+1,j,k) + g * tmp1
    v(i+2,j,k) = v(i+2,j,k) + g * tmp2
    v(i+3,j,k) = v(i+3,j,k) + g * tmp3
    v(i+4,j,k) = v(i+4,j,k) + g * tmp4
    v(i+5,j,k) = v(i+5,j,k) + g * tmp5
  end do
  end do

#elif __NA__ == __NA_T8__ + 7

  ! remainder: i = na-6:na ...................................................

  i = __NA_T8__ + 1

  !$acc loop collapse(2) vector
  do k = 1, __NA__
  do i = 1, __NA__
    tmp0 = 0
    tmp1 = 0
    tmp2 = 0
    tmp3 = 0
    tmp4 = 0
    tmp5 = 0
    tmp6 = 0
    do p = 1, __NA__
      tmp0 = tmp0 + A(p,i  ) * u(p,j,k)
      tmp1 = tmp1 + A(p,i+1) * u(p,j,k)
      tmp2 = tmp2 + A(p,i+2) * u(p,j,k)
      tmp3 = tmp3 + A(p,i+3) * u(p,j,k)
      tmp4 = tmp4 + A(p,i+4) * u(p,j,k)
      tmp5 = tmp5 + A(p,i+5) * u(p,j,k)
      tmp6 = tmp6 + A(p,i+6) * u(p,j,k)
    end do
    v(i  ,j,k) = v(i  ,j,k) + g * tmp0
    v(i+1,j,k) = v(i+1,j,k) + g * tmp1
    v(i+2,j,k) = v(i+2,j,k) + g * tmp2
    v(i+3,j,k) = v(i+3,j,k) + g * tmp3
    v(i+4,j,k) = v(i+4,j,k) + g * tmp4
    v(i+5,j,k) = v(i+5,j,k) + g * tmp5
    v(i+6,j,k) = v(i+6,j,k) + g * tmp6
  end do
  end do

#endif

end subroutine SubOp_1a
