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
!> Computes  v = α Qᵀ⊗I⊗I u + β v
!>
!>   * k:  block 1,  unroll 1
!>   * j:  block 1,  unroll 1
!>   * i:  block 1,  unroll 1
!>   * p:  block 1,  unroll 1
!>   * using Intel SIMD directive

subroutine PROC(QtxIxI__,_NQ_)(Q, alpha, beta, u, v)
  !$acc routine vector
  real(RWP), intent(in)    :: Q(_NQ_,_NQ_)      !< square matrix
  real(RWP), intent(in)    :: alpha             !< factor α
  real(RWP), intent(in)    :: beta              !< factor β
  real(RWP), intent(in)    :: u(_NQ_*_NQ_,_NQ_) !< operand
  real(RWP), intent(inout) :: v(_NQ_*_NQ_,_NQ_) !< result

  real(RWP) :: tmp
  integer   :: ij, k, p

  !$acc loop collapse(2) vector
  do k = 1, _NQ_
    !DIR$ SIMD VECREMAINDER
    do ij = 1, _NQ_**2
      tmp = 0
      do p = 1, _NQ_
        tmp = tmp + Q(p,k) * u(ij,p)
      end do
      v(ij,k) = alpha * tmp + beta * v(ij,k)
    end do
  end do

end subroutine PROC(QtxIxI__,_NQ_)
