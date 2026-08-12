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
!> Computes  v = α I⊗I⊗Aᵀ u + β v
!>
!>   * k:  block 1,  unroll 1
!>   * j:  block 1,  unroll 1
!>   * i:  block 4,  unroll 1
!>   * p:  block 1,  unroll 4
!>   * using Intel SIMD directive

#define _NA1_T4_ (_NA1_ / 4) * 4
#define _NA2_T4_ (_NA2_ / 4) * 4

subroutine PROC(IxIxAt__,_NA1_,_NA2_)(nb, nc, A, alpha, beta, u, v)
  !$acc routine vector
  integer,   intent(in)    :: nb              !< 2nd dimension of u,v
  integer,   intent(in)    :: nc              !< 3rd dimension of u,v
  real(RWP), intent(in)    :: A(_NA1_,_NA2_)  !< rectangular matrix
  real(RWP), intent(in)    :: alpha           !< factor α
  real(RWP), intent(in)    :: beta            !< factor β
  real(RWP), intent(in)    :: u(_NA1_,nb,nc)  !< operand
  real(RWP), intent(inout) :: v(_NA2_,nb,nc)  !< result

  real(RWP) :: tmp(0:3)
  integer   :: i, j, k, p

  if (beta == 0) then
    do k = 1, nc
    do j = 1, nb
    do i = 1, _NA2_
     v(i,j,k) = 0
    end do
    end do
    end do
  else
    do k = 1, nc
    do j = 1, nb
    do i = 1, _NA2_
     v(i,j,k) = beta * v(i,j,k)
    end do
    end do
    end do
  end if

#if _NA1_T4_ > 0

  !$acc loop independent vector
  do k = 1, nc
    !$acc loop seq
    do p = 1, _NA1_T4_, 4
      !$acc loop collapse(2) independent vector private(tmp)
      do j = 1, nb
      do i = 1, _NA2_T4_, 4

        tmp =       A(p  ,i:i+3) * u(p  ,j,k)
        tmp = tmp + A(p+1,i:i+3) * u(p+1,j,k)
        tmp = tmp + A(p+2,i:i+3) * u(p+2,j,k)
        tmp = tmp + A(p+3,i:i+3) * u(p+3,j,k)

        v(i:i+3,j,k) = alpha * tmp + v(i:i+3,j,k)

      end do
      end do
    end do
  end do

#endif

#if _NA1_ == _NA1_T4_ + 1

  ! 1:na1/4*4, na2 -------------------------------------------------------------

  p = _NA1_T4_ + 1

  !$acc loop collapse(3) independent vector
  do k = 1, nc
    do j = 1, nb
    do i = 1, _NA2_T4_, 4

      v(i:i+3,j,k) = v(i:i+3,j,k) + alpha * A(p,i:i+3) * u(p,j,k)

    end do
    end do
  end do

#elif _NA1_ == _NA1_T4_ + 2

  ! 1:na1/4*4, na2-1:na2 -------------------------------------------------------

  p = _NA1_T4_ + 1

  !$acc loop collapse(3) independent vector private(tmp)
  do k = 1, nc
    do j = 1, nb
    do i = 1, _NA2_T4_, 4

      tmp =       A(p  ,i:i+3) * u(p  ,j,k)
      tmp = tmp + A(p+1,i:i+3) * u(p+1,j,k)

      v(i:i+3,j,k) = alpha * tmp + v(i:i+3,j,k)

    end do
    end do
  end do

#elif _NA1_ == _NA1_T4_ + 3

  ! 1:na1/4*4, na2-2:na2 -------------------------------------------------------

  p = _NA1_T4_ + 1

  !$acc loop collapse(3) independent vector private(tmp)
  do k = 1, nc
    do j = 1, nb
    do i = 1, _NA2_T4_, 4

      tmp =       A(p  ,i:i+3) * u(p  ,j,k)
      tmp = tmp + A(p+1,i:i+3) * u(p+1,j,k)
      tmp = tmp + A(p+2,i:i+3) * u(p+2,j,k)

      v(i:i+3,j,k) = alpha * tmp + v(i:i+3,j,k)

    end do
    end do
  end do

#endif

#if _NA2_ == _NA2_T4_ + 1

  ! na1, 1:na2 -----------------------------------------------------------------

  i = _NA2_T4_ + 1

  !$acc loop collapse(2) independent vector
  do k = 1, nc
  do j = 1, nb
    do p = 1, _NA1_

      v(i,j,k) = v(i,j,k) + alpha * A(p,i) * u(p,j,k)

    end do
  end do
  end do

#elif _NA2_ == _NA2_T4_ + 2

  ! na1-1:na1, 1:na2 -----------------------------------------------------------

  i = _NA2_T4_ + 1

  !$acc loop collapse(2) independent vector
  do k = 1, nc
  do j = 1, nb
    do p = 1, _NA1_

      v(i:i+1,j,k) = v(i:i+1,j,k) + alpha * A(p,i:i+1) * u(p,j,k)

    end do
  end do
  end do

#elif _NA2_ == _NA2_T4_ + 3

  ! na1-2:na1, 1:na2 -----------------------------------------------------------

  i = _NA1_T4_ + 1

  !$acc loop collapse(2) independent vector
  do k = 1, nc
  do j = 1, nb
    do p = 1, _NA1_

      v(i:i+2,j,k) = v(i:i+2,j,k) + alpha * A(p,i:i+2) * u(p,j,k)

    end do
  end do
  end do

#endif

end subroutine PROC(IxIxAt__,_NA1_,_NA2_)

#undef _NA1_T4_
#undef _NA2_T4_
