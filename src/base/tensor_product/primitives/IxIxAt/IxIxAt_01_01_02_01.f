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

!> Computes  v = α I⊗I⊗Aᵀ u + β v
!>
!>   * k:  block 1,  unroll 1
!>   * j:  block 1,  unroll 1
!>   * i:  block 1,  unroll 2
!>   * p:  block 1,  unroll 1
!>   * using Intel SIMD directive

#define _NA2_T2_ (_NA2_ / 2) * 2

subroutine PROC(IxIxAt__,_NA1_,_NA2_)(nb, nc, A, alpha, beta, u, v)
  !$acc routine vector
  integer,   intent(in)    :: nb              !< 2nd dimension of u,v
  integer,   intent(in)    :: nc              !< 3rd dimension of u,v
  real(RWP), intent(in)    :: A(_NA1_,_NA2_)  !< rectangular matrix
  real(RWP), intent(in)    :: alpha           !< factor α
  real(RWP), intent(in)    :: beta            !< factor β
  real(RWP), intent(in)    :: u(_NA1_,nb,nc)  !< operand
  real(RWP), intent(inout) :: v(_NA2_,nb,nc)  !< result

  real(RWP) :: tmp0, tmp1
  integer   :: i, j, k, p

#if _NA2_T2_ > 0

  do k = 1, nc
  do j = 1, nb
    do i = 1, _NA2_T2_, 2
     tmp0 = 0
     tmp1 = 0
      do p = 1, _NA1_
        tmp0 = tmp0 + A(p,i  ) * u(p,j,k)
        tmp1 = tmp1 + A(p,i+1) * u(p,j,k)
      end do
      v(i  ,j,k) = alpha * tmp0 + beta * v(i  ,j,k)
      v(i+1,j,k) = alpha * tmp1 + beta * v(i+1,j,k)
    end do
  end do
  end do

#endif

#if _NA2_ == _NA2_T2_ + 1

  i = _NA2_

  !$acc loop collapse(2) vector
  do k = 1, nc
    do j = 1, nb
      tmp0 = 0
      do p = 1, _NA1_
        tmp0 = tmp0 + A(p,i) * u(p,j,k)
      end do
      v(i,j,k) = alpha * tmp0 + beta * v(i,j,k)
    end do
  end do

#endif

end subroutine PROC(IxIxAt__,_NA1_,_NA2_)

#undef _NA2_T2_
