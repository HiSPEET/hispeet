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
!> Computes  v = α Cᵀ⊗I⊗I u + β v
!>
!>   * k:  block 1,  unroll 1
!>   * j:  joined with i
!>   * i:  joined with j
!>   * p:  block 1,  unroll 1
!>   * using Intel SIMD directive

subroutine PROC(CtxIxI__,_NC1_,_NC2_)(na_nb, C, alpha, beta, u, v)
  !$acc routine vector
  integer,   intent(in)    :: na_nb           !< size of dimensions 1+2 of u,v
  real(RWP), intent(in)    :: C(_NC1_,_NC2_)  !< rectangular matrix
  real(RWP), intent(in)    :: alpha           !< factor α
  real(RWP), intent(in)    :: beta            !< factor β
  real(RWP), intent(in)    :: u(na_nb,_NC1_)  !< operand
  real(RWP), intent(inout) :: v(na_nb,_NC2_)  !< result

  real(RWP) :: tmp
  integer   :: ij, k, p

  !$acc loop collapse(2) vector
  do k = 1, _NC2_
    !DIR$ SIMD
    do ij = 1, na_nb
      tmp = 0
      do p = 1, _NC1_
        tmp = tmp + C(p,k) * u(ij,p)
      end do
      v(ij,k) = alpha * tmp + beta * v(ij,k)
    end do
  end do

end subroutine PROC(CtxIxI__,_NC1_,_NC2_)
