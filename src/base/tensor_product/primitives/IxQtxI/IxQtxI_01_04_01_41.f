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

!-----------------------------------------------------------------------------
!> Computes  v = α I⊗Qᵀ⊗I u + β v
!>
!>   * k:  block 1,  unroll 1
!>   * j:  block 1,  unroll 4
!>   * i:  block 1,  unroll 1
!>   * p:  block 4,  unroll 1
!>   * explicit remainder handling

#define _NQ_T4_ (_NQ_ / 4) * 4

subroutine PROC(IxQtxI__,_NQ_)(Q, alpha, beta, u, v)
  !$acc routine vector
  real(RWP), intent(in)    :: Q(_NQ_,_NQ_)      !< square matrix
  real(RWP), intent(in)    :: alpha             !< factor α
  real(RWP), intent(in)    :: beta              !< factor β
  real(RWP), intent(in)    :: u(_NQ_,_NQ_,_NQ_) !< operand
  real(RWP), intent(inout) :: v(_NQ_,_NQ_,_NQ_) !< result

  integer, parameter :: NQ_T4 = _NQ_T4_

  real(RWP) :: tmp0, tmp1, tmp2, tmp3
  integer   :: i, j, k, p   ! loop counters
  integer   :: pb           ! block counters

  !$acc loop collapse(3)
  do k = 1, _NQ_
  do j = 1, _NQ_
  do i = 1, _NQ_
    v(i,j,k) = beta * v(i,j,k)
  end do
  end do
  end do

#if _NQ_T4_ > 0

  !$acc loop independent vector
  do k = 1, _NQ_
  do j = 1, NQ_T4, 4
    !$acc loop seq
    do pb = 1, NQ_T4, 4
      !$acc loop collapse(2) independent vector
      do i = 1, _NQ_
        tmp0 = 0
        tmp1 = 0
        tmp2 = 0
        tmp3 = 0
        do p = pb, pb+3
          tmp0 = tmp0 + Q(p,j  ) * u(i,p,k)
          tmp1 = tmp1 + Q(p,j+1) * u(i,p,k)
          tmp2 = tmp2 + Q(p,j+2) * u(i,p,k)
          tmp3 = tmp3 + Q(p,j+3) * u(i,p,k)
        end do
        v(i,j  ,k) = alpha * tmp0 + v(i,j  ,k)
        v(i,j+1,k) = alpha * tmp1 + v(i,j+1,k)
        v(i,j+2,k) = alpha * tmp2 + v(i,j+2,k)
        v(i,j+3,k) = alpha * tmp3 + v(i,j+3,k)
      end do
    end do
  end do
  end do

#endif

#if _NQ_ == _NQ_T4_ + 1

  !$acc loop independent vector
  do k = 1, _NQ_

    ! remainder: j = 1: nq-1, p = nq .........................................

    p = _NQ_

    !$acc loop collapse(2) vector
    do j = 1, NQ_T4, 4
    do i = 1, _NQ_
      v(i,j  ,k) = alpha * Q(p,j  ) * u(i,p,k) + v(i,j  ,k)
      v(i,j+1,k) = alpha * Q(p,j+1) * u(i,p,k) + v(i,j+1,k)
      v(i,j+2,k) = alpha * Q(p,j+2) * u(i,p,k) + v(i,j+2,k)
      v(i,j+3,k) = alpha * Q(p,j+3) * u(i,p,k) + v(i,j+3,k)
    end do
    end do

    ! remainder: j = nq, p = 1:nq  ...........................................

    !$acc loop vector
    do i = 1, _NQ_
      tmp0 = 0
      do p = 1, _NQ_
        tmp0 = tmp0 + Q(p,_NQ_) * u(i,p,k)
      end do
      v(i,_NQ_,k) = alpha * tmp0 + v(i,_NQ_,k)
    end do

  end do

#elif _NQ_ == _NQ_T4_ + 2

  !$acc loop independent vector
  do k = 1, _NQ_

    ! remainder: j = 1:nq-2, p = nq-1:nq .....................................

    !$acc loop collapse(2) vector
    do j = 1, NQ_T4, 4
    do i = 1, _NQ_
      do p = _NQ_-1, _NQ_
        v(i,j  ,k) = alpha * Q(p,j  ) * u(i,p,k) + v(i,j  ,k)
        v(i,j+1,k) = alpha * Q(p,j+1) * u(i,p,k) + v(i,j+1,k)
        v(i,j+2,k) = alpha * Q(p,j+2) * u(i,p,k) + v(i,j+2,k)
        v(i,j+3,k) = alpha * Q(p,j+3) * u(i,p,k) + v(i,j+3,k)
      end do
    end do
    end do

    ! remainder: j = nq-1:nq, p = 1:nq .......................................

    !$acc loop vector
    do i = 1, _NQ_
      tmp0 = 0
      tmp1 = 0
      do p = 1, _NQ_
        tmp0 = tmp0 + Q(p,_NQ_-1) * u(i,p,k)
        tmp1 = tmp1 + Q(p,_NQ_  ) * u(i,p,k)
      end do
      v(i,_NQ_-1,k) = alpha * tmp0 + v(i,_NQ_-1,k)
      v(i,_NQ_  ,k) = alpha * tmp1 + v(i,_NQ_  ,k)
    end do

  end do

#elif _NQ_ == _NQ_T4_ + 3

  !$acc loop independent vector
  do k = 1, _NQ_

    ! remainder: j = 1:nq-3, p = nq-2:nq .....................................

    !$acc loop collapse(2) vector
    do j = 1, NQ_T4, 4
    do i = 1, _NQ_
      do p = _NQ_-2, _NQ_
        v(i,j  ,k) = alpha * Q(p,j  ) * u(i,p,k) + v(i,j  ,k)
        v(i,j+1,k) = alpha * Q(p,j+1) * u(i,p,k) + v(i,j+1,k)
        v(i,j+2,k) = alpha * Q(p,j+2) * u(i,p,k) + v(i,j+2,k)
        v(i,j+3,k) = alpha * Q(p,j+3) * u(i,p,k) + v(i,j+3,k)
      end do
    end do
    end do

    ! remainder: j = nq-2:nq, p = 1:nq .......................................

    !$acc loop vector
    do i = 1, _NQ_
      tmp0 = 0
      tmp1 = 0
      tmp2 = 0
      do p = 1, _NQ_
        tmp0 = tmp0 + Q(p,_NQ_-2) * u(i,p,k)
        tmp1 = tmp1 + Q(p,_NQ_-1) * u(i,p,k)
        tmp2 = tmp2 + Q(p,_NQ_  ) * u(i,p,k)
      end do
      v(i,_NQ_-2,k) = alpha * tmp0 + v(i,_NQ_-2,k)
      v(i,_NQ_-1,k) = alpha * tmp1 + v(i,_NQ_-1,k)
      v(i,_NQ_  ,k) = alpha * tmp2 + v(i,_NQ_  ,k)
    end do

  end do

#endif

end subroutine PROC(IxQtxI__,_NQ_)

#undef _NQ_T4
