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
!>
!>   * k: simple loop (l)
!>   * j: simple loop (l)
!>   * i: unrolled and jammed with length of 4 (u4)
!>   * p: simple loop (l)
!>   * explicit remainder handling (r)

subroutine SubOp_1d(d0, d1, d2, d3, V1, V2, V3, A, u, v)
  !$acc routine vector
  real(RNP), intent(in)  :: d0                       !< λ Δx₁ Δx₂ Δx₃
  real(RNP), intent(in)  :: d1                       !< ν Δx₂ Δx₃ / Δx₁
  real(RNP), intent(in)  :: d2                       !< ν Δx₃ Δx₁ / Δx₂
  real(RNP), intent(in)  :: d3                       !< ν Δx₁ Δx₂ / Δx₃
  real(RNP), intent(in)  :: V1(__NA__)               !< 1D eigenvalues in dir 1
  real(RNP), intent(in)  :: V2(__NA__)               !< 1D eigenvalues in dir 2
  real(RNP), intent(in)  :: V3(__NA__)               !< 1D eigenvalues in dir 3
  real(RNP), intent(in)  :: A(__NA__,__NA__)         !< A
  real(RNP), intent(in)  :: u(__NA__,__NA__,__NA__)  !< operand
  real(RNP), intent(out) :: v(__NA__,__NA__,__NA__)  !< result

  real(RNP) :: tmp0, tmp1, tmp2, tmp3
  integer   :: i, j, k, p

#if __NA_T4__ > 0

  !$acc loop collapse(3) vector
  do k = 1, __NA__
  do j = 1, __NA__
  do i = 1, __NA_T4__, 4
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
    v(i  ,j,k) = tmp0
    v(i+1,j,k) = tmp1
    v(i+2,j,k) = tmp2
    v(i+3,j,k) = tmp3
  end do
  end do
  end do

#endif

#if __NA__ == __NA_T4__ + 1

  ! remainder: i = np ........................................................

  !$acc loop collapse(2) vector
  do k = 1, __NA__
  do j = 1, __NA__
    tmp0 = 0
    do p = 1, __NA__
      tmp0 = tmp0 + A(p,__NA__) * u(p,j,k)
    end do
    v(__NA__,j,k) = tmp0
  end do
  end do

#elif __NA__ == __NA_T4__ + 2

  ! remainder: i = np-1:np ...................................................

  !$acc loop collapse(2) vector
  do k = 1, __NA__
  do j = 1, __NA__
    tmp0 = 0
    tmp1 = 0
    do p = 1, __NA__
      tmp0 = tmp0 + A(p,__NA__-1) * u(p,j,k)
      tmp1 = tmp1 + A(p,__NA__  ) * u(p,j,k)
    end do
    v(__NA__-1,j,k) = tmp0
    v(__NA__  ,j,k) = tmp1
  end do
  end do

#elif __NA__ == __NA_T4__ + 3

  ! remainder: i = np-2:np ...................................................

  !$acc loop collapse(2) vector
  do k = 1, __NA__
  do j = 1, __NA__
    tmp0 = 0
    tmp1 = 0
    tmp2 = 0
    do p = 1, __NA__
      tmp0 = tmp0 + A(p,__NA__-2) * u(p,j,k)
      tmp1 = tmp1 + A(p,__NA__-1) * u(p,j,k)
      tmp2 = tmp2 + A(p,__NA__  ) * u(p,j,k)
    end do
    v(__NA__-2,j,k) = tmp0
    v(__NA__-1,j,k) = tmp1
    v(__NA__  ,j,k) = tmp2
  end do
  end do

#endif

  !$acc loop collapse(3) vector
  do k = 1, __NA__
  do j = 1, __NA__
  do i = 1, __NA__
    v(i,j,k) = v(i,j,k) / (d0 + d1*V1(i) + d2*V2(j) + d3*V3(k))
  end do
  end do
  end do

end subroutine SubOp_1d

