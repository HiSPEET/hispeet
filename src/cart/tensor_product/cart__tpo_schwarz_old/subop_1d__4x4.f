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
!> Performs `v =(d₀IxIxI + d₁IxIxV₁ + d₂IxIxV₂ + d₃IxIxV₃)⁻¹ (IxIxA') u`
!> for ns = 4
!>
!>   * k: simple loop (l)
!>   * j: simple loop (l)
!>   * i: blocked with length of 4 (b4)
!>   * p: unrolled and jammed with length of 4 (u4)
!>   * no remainder handling

subroutine SubOp_1d(d0, d1, d2, d3, V1, V2, V3, A, u, v)
  !$acc routine vector
  real(RNP), intent(in)  :: d0        !< λ Δx₁ Δx₂ Δx₃
  real(RNP), intent(in)  :: d1        !< ν Δx₂ Δx₃ / Δx₁
  real(RNP), intent(in)  :: d2        !< ν Δx₃ Δx₁ / Δx₂
  real(RNP), intent(in)  :: d3        !< ν Δx₁ Δx₂ / Δx₃
  real(RNP), intent(in)  :: V1(4)     !< 1D eigenvalues in direction 1
  real(RNP), intent(in)  :: V2(4)     !< 1D eigenvalues in direction 2
  real(RNP), intent(in)  :: V3(4)     !< 1D eigenvalues in direction 3
  real(RNP), intent(in)  :: A(4,4)    !< A
  real(RNP), intent(in)  :: u(4,4,4)  !< operand
  real(RNP), intent(out) :: v(4,4,4)  !< result

  real(RNP) :: tmp(0:3)
  integer   :: j, k

  !$acc loop collapse(2) independent vector private(tmp)
  do k = 1, 4
    do j = 1, 4

      tmp =       A(1,1:4) * u(1,j,k)
      tmp = tmp + A(2,1:4) * u(2,j,k)
      tmp = tmp + A(3,1:4) * u(3,j,k)
      tmp = tmp + A(4,1:4) * u(4,j,k)

      v(1:4,j,k) = tmp / (d0 + d1*V1(1:4) + d2*V2(j) + d3*V3(k))

    end do
  end do

end subroutine SubOp_1d
