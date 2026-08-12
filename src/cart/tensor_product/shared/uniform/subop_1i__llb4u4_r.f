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
!> Performs `v = g IxIxAᵀ u`
!>
!>   * k: simple loop (l)
!>   * j: simple loop (l)
!>   * i: blocked with length of 4 (b4)
!>   * p: unrolled and jammed with length of 4 (u4)
!>   * explicit remainder handling (r)

subroutine SubOp_1i(g, A, u, v)
  !$acc routine vector
  real(RNP), intent(in)  :: g                       !< metric factor
  real(RNP), intent(in)  :: A(__NA__,__NA__)        !< A
  real(RNP), intent(in)  :: u(__NA__,__NA__,__NA__) !< operand
  real(RNP), intent(out) :: v(__NA__,__NA__,__NA__) !< result

  real(RNP) :: tmp(0:3)
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

  !$acc loop independent vector
  do k = 1, __NA__
    !$acc loop seq
    do q = 1, __NA_T4__, 4
      !$acc loop collapse(2) independent vector private(tmp)
      do j = 1, __NA__
      do i = 1, __NA_T4__, 4

        tmp =       A(q  ,i:i+3) * u(q  ,j,k)
        tmp = tmp + A(q+1,i:i+3) * u(q+1,j,k)
        tmp = tmp + A(q+2,i:i+3) * u(q+2,j,k)
        tmp = tmp + A(q+3,i:i+3) * u(q+3,j,k)

        v(i:i+3,j,k) = v(i:i+3,j,k) + g * tmp

      end do
      end do
    end do
  end do

#endif

#if __NA__ == __NA_T4__ + 1

  ! 1:na/4*4, na ---------------------------------------------------------------

  q = __NA_T4__ + 1

  !$acc loop collapse(3) independent vector private(tmp)
  do k = 1, __NA__
    do j = 1, __NA__
    do i = 1, __NA_T4__, 4

      v(i:i+3,j,k) = v(i:i+3,j,k) + g * A(q,i:i+3) * u(q,j,k)

    end do
    end do
  end do

  ! na, 1:na -------------------------------------------------------------------

  i = __NA_T4__ + 1

  !$acc loop collapse(2) independent vector
  do k = 1, __NA__
  do j = 1, __NA__
    do p = 1, __NA__

      v(i,j,k) = v(i,j,k) + g * A(p,i) * u(p,j,k)

    end do
  end do
  end do


#elif __NA__ == __NA_T4__ + 2

  ! 1:na/4*4, na-2:na ----------------------------------------------------------

  q = __NA_T4__ + 1

  !$acc loop collapse(3) independent vector private(tmp)
  do k = 1, __NA__
    do j = 1, __NA__
    do i = 1, __NA_T4__, 4

      tmp =       A(q  ,i:i+3) * u(q  ,j,k)
      tmp = tmp + A(q+1,i:i+3) * u(q+1,j,k)

      v(i:i+3,j,k) = v(i:i+3,j,k) + g * tmp

    end do
    end do
  end do

  ! na-2:na, 1:na --------------------------------------------------------------

  i = __NA_T4__ + 1

  !$acc loop collapse(2) independent vector
  do k = 1, __NA__
  do j = 1, __NA__
    do p = 1, __NA__

      v(i:i+1,j,k) = v(i:i+1,j,k) + g * A(p,i:i+1) * u(p,j,k)

    end do
  end do
  end do

#elif __NA__ == __NA_T4__ + 3

  ! 1:na/4*4, na-2:na ----------------------------------------------------------

  q = __NA_T4__ + 1

  !$acc loop collapse(3) independent vector private(tmp)
  do k = 1, __NA__
    do j = 1, __NA__
    do i = 1, __NA_T4__, 4

      tmp =       A(q  ,i:i+3) * u(q  ,j,k)
      tmp = tmp + A(q+1,i:i+3) * u(q+1,j,k)
      tmp = tmp + A(q+2,i:i+3) * u(q+2,j,k)

      v(i:i+3,j,k) = v(i:i+3,j,k) + g * tmp

    end do
    end do
  end do

  ! na-2:na, 1:na --------------------------------------------------------------

  i = __NA_T4__ + 1

  !$acc loop collapse(2) independent vector
  do k = 1, __NA__
  do j = 1, __NA__
    do p = 1, __NA__

      v(i:i+2,j,k) = v(i:i+2,j,k) + g * A(p,i:i+2) * u(p,j,k)

    end do
  end do
  end do

#endif

end subroutine SubOp_1i
