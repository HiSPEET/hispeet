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
!> Computes  v = α I⊗I⊗Qᵀ u + β v
!>
!>   * k:  block 1,  unroll 1
!>   * j:  block 1,  unroll 1
!>   * i:  block 4,  unroll 1
!>   * p:  block 4,  unroll 4
!>   * explicit remainder handling

#define _NQ_T4_ (_NQ_ / 4) * 4

subroutine PROC(IxIxQt__,_NQ_)(Q, alpha, beta, u, v)
  !$acc routine vector
  real(RWP), intent(in)    :: Q(_NQ_,_NQ_)      !< square matrix
  real(RWP), intent(in)    :: alpha             !< factor α
  real(RWP), intent(in)    :: beta              !< factor β
  real(RWP), intent(in)    :: u(_NQ_,_NQ_,_NQ_) !< operand
  real(RWP), intent(inout) :: v(_NQ_,_NQ_,_NQ_) !< result

  integer, parameter :: NQ_T4 = _NQ_T4_

  real(RWP) :: tmp
  integer   :: i, j, k, p   ! loop counters
  integer   :: ib, pb       ! block counters

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
    !$acc loop seq
    do pb = 1, NQ_T4, 4
      !$acc loop collapse(2) independent vector private(tmp)
      do j = 1, _NQ_
      do ib = 1, NQ_T4, 4
        do i = ib, ib+3

          tmp =       Q(pb  ,i) * u(pb  ,j,k)
          tmp = tmp + Q(pb+1,i) * u(pb+1,j,k)
          tmp = tmp + Q(pb+2,i) * u(pb+2,j,k)
          tmp = tmp + Q(pb+3,i) * u(pb+3,j,k)

          v(i,j,k) = alpha * tmp + v(i,j,k)

        end do
      end do
      end do
    end do
  end do

#endif

#if _NQ_ == _NQ_T4_ + 1

  ! 1:nq/4*4, nq -----------------------------------------------------------

  pb = NQ_T4 + 1

  !$acc loop collapse(3) independent vector private(tmp)
  do k = 1, _NQ_
    do j = 1, _NQ_
    do ib = 1, NQ_T4, 4
      do i = ib, ib+3

        v(i,j,k) = alpha * Q(pb,i) * u(pb,j,k) + v(i,j,k)

      end do
    end do
    end do
  end do

  ! nq, 1:nq ---------------------------------------------------------------

  ib = NQ_T4 + 1
  i  = ib

  !$acc loop collapse(2) independent vector
  do k = 1, _NQ_
    do j = 1, _NQ_
      do p = 1, _NQ_

        v(i,j,k) = alpha * Q(p,i) * u(p,j,k) + v(i,j,k)

      end do
    end do
  end do

#elif _NQ_ == _NQ_T4_ + 2

  ! 1:nq/4*4, nq-2:nq ----------------------------------------------------

  pb = NQ_T4 + 1

  !$acc loop collapse(3) independent vector private(tmp)
  do k = 1, _NQ_
    do j = 1, _NQ_
    do ib = 1, NQ_T4, 4
      do i = ib, ib+3

        tmp =       Q(pb  ,i) * u(pb  ,j,k)
        tmp = tmp + Q(pb+1,i) * u(pb+1,j,k)

        v(i,j,k) = alpha * tmp + v(i,j,k)

      end do
    end do
    end do
  end do

  ! nq-2:nq, 1:nq --------------------------------------------------------

  ib = NQ_T4 + 1

  !$acc loop collapse(2) independent vector
  do k = 1, _NQ_
    do j = 1, _NQ_
    do i = ib, ib+1
      do p = 1, _NQ_

        v(i,j,k) = alpha * Q(p,i) * u(p,j,k) + v(i,j,k)

      end do
    end do
    end do
  end do

#elif _NQ_ == _NQ_T4_ + 3

  ! 1:nq/4*4, nq-2:nq ----------------------------------------------------

  pb = NQ_T4 + 1

  !$acc loop collapse(3) independent vector private(tmp)
  do k = 1, _NQ_
    do j = 1, _NQ_
    do ib = 1, NQ_T4, 4
      do i = ib, ib+3

        tmp =       Q(pb  ,i) * u(pb  ,j,k)
        tmp = tmp + Q(pb+1,i) * u(pb+1,j,k)
        tmp = tmp + Q(pb+2,i) * u(pb+2,j,k)

        v(i,j,k) = alpha * tmp + v(i,j,k)

      end do
    end do
    end do
  end do

  ! nq-2:nq, 1:nq --------------------------------------------------------

  ib = NQ_T4 + 1

  !$acc loop collapse(2) independent vector
  do k = 1, _NQ_
    do j = 1, _NQ_
    do i = ib, ib+2
      do p = 1, _NQ_

        v(i,j,k) = alpha * Q(p,i) * u(p,j,k) + v(i,j,k)

      end do
    end do
    end do
  end do

#endif

end subroutine PROC(IxIxQt__,_NQ_)

#undef _NQ_T4_
