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

!> summary:  Cubic and quintic Hermite polynomials and interpolants
!> author:   Joerg Stiller
!> date:     2014/03/12
!>
!>### Cubic and quintic Hermite polynomials and interpolants
!===============================================================================

module Hermite_Interpolation
  use Kind_Parameters, only: RNP
  use Constants,       only: HALF
  implicit none
  private

  public :: CubicHermitePolynomial
  public :: CubicHermiteInterpolant

  public :: QuinticHermitePolynomial
  public :: QuinticHermiteInterpolant

  !----------------------------------------------------------------------------
  !> Cubic Hermite interpolant h(x) = sum(i=0:3) c(i) H³_i(x)
  !>
  !> Coefficients:
  !>
  !>     c(0) = h (0)
  !>     c(1) = h'(0)
  !>     c(2) = h (1)
  !>     c(3) = h'(1)

  interface CubicHermiteInterpolant
    module procedure CubicHermiteInterpolant00  ! scalar coeff., single pos.
    module procedure CubicHermiteInterpolant01  ! scalar coeff., vector pos.
  end interface

  !----------------------------------------------------------------------------
  !> Quintic Hermite interpolant h(x) = sum(i=0:5) c(i) H⁵_i(x)
  !>
  !> Coefficients:
  !>
  !>     c(0) = h (0)
  !>     c(1) = h'(0)
  !>     c(2) = h"(0)
  !>     c(3) = h (1)
  !>     c(4) = h'(1)
  !>     c(5) = h"(1)

  interface QuinticHermiteInterpolant
    module procedure QuinticHermiteInterpolant00  ! scalar coeff., single pos.
    module procedure QuinticHermiteInterpolant01  ! scalar coeff., vector pos.
  end interface

contains

!===============================================================================
! Cubic case

!-------------------------------------------------------------------------------
!> Cubic Hermite polynomials H³_i(x), i = 0 ... 3
!>
!> The cubic Hermite polynomials satisfy the conditions
!>
!>     H₀(0) = 1,  H₀(1) = 0,  H₀'(0) = 0,  H₀'(1) = 0
!>     H₁(0) = 0,  H₁(1) = 1,  H₁'(0) = 0,  H₁'(1) = 0
!>     H₂(0) = 0,  H₂(1) = 0,  H₂'(0) = 1,  H₂'(1) = 0
!>     H₃(0) = 0,  H₃(1) = 0,  H₃'(0) = 0,  H₃'(1) = 1

pure real(RNP) function CubicHermitePolynomial(i, x) result(h)
  integer,   intent(in) :: i  !< ID ranging from 0 to 3
  real(RNP), intent(in) :: x  !< 0 <= x <= 1

  select case(i)
  case(0)
     h = 1 - 3*x**2 + 2*x**3
  case(1)
     h = 3*x**2 - 2*x**3
  case(2)
     h = x - 2*x**2 + x**3
  case(3)
     h = -x**2 + x**3
  case default
     h = 0
  end select

end function CubicHermitePolynomial

!-------------------------------------------------------------------------------
!> Cubic Hermite interpolant for scalar coefficients and single argument

pure real(RNP) function CubicHermiteInterpolant00(c, x) result(h)
  real(RNP), intent(in) :: c(0:3) !< coefficients
  real(RNP), intent(in) :: x      !< argument

  h = c(0) * (1 - 3*x**2 + 2*x**3)  &
    + c(1) * (3*x**2 - 2*x**3)      &
    + c(2) * (x - 2*x**2 + x**3)    &
    + c(3) * (-x**2 + x**3)

end function CubicHermiteInterpolant00

!-------------------------------------------------------------------------------
!> Cubic Hermite interpolant for scalar coefficients and multiple arguments

pure function CubicHermiteInterpolant01(c, x) result(h)
  real(RNP), intent(in) :: c(0:3) !< coefficients
  real(RNP), intent(in) :: x(:)   !< arguments

  real(RNP), dimension(size(x)) :: h

  h = c(0) * (1 - 3*x**2 + 2*x**3)  &
    + c(1) * (3*x**2 - 2*x**3)      &
    + c(2) * (x - 2*x**2 + x**3)    &
    + c(3) * (-x**2 + x**3)

end function CubicHermiteInterpolant01

!===============================================================================
! Quintic case

!-------------------------------------------------------------------------------
!> Quintic Hermite polynomials H⁵_i(x), i = 0 ... 5
!>
!> The quintic Hermite polynomials satisfy the conditions
!>
!>     H₀(0) = 1,  H₀(1) = 0,  H₀'(0) = 0,  H₀'(1) = 0,  H₀"(0) = 0,  H₀"(1) = 0
!>     H₁(0) = 0,  H₁(1) = 1,  H₁'(0) = 0,  H₁'(1) = 0,  H₁"(0) = 0,  H₁"(1) = 0
!>     H₂(0) = 0,  H₂(1) = 0,  H₂'(0) = 1,  H₂'(1) = 0,  H₂"(0) = 0,  H₂"(1) = 0
!>     H₃(0) = 0,  H₃(1) = 0,  H₃'(0) = 0,  H₃'(1) = 1,  H₃"(0) = 0,  H₃"(1) = 0
!>     H₄(0) = 0,  H₄(1) = 0,  H₄'(0) = 0,  H₄'(1) = 0,  H₄"(0) = 1,  H₄"(1) = 0
!>     H₅(0) = 0,  H₅(1) = 0,  H₅'(0) = 0,  H₅'(1) = 0,  H₅"(0) = 0,  H₅"(1) = 1

pure real(RNP) function QuinticHermitePolynomial(i, x) result(h)
  integer,   intent(in) :: i  !< ID ranging from 0 to 5
  real(RNP), intent(in) :: x  !< 0 <= x <= 1

  select case(i)
  case(0)
     h = 1 - 10*x**3 + 15*x**4 - 6*x**5
  case(1)
     h = 10*x**3 - 15*x**4 + 6*x**5
  case(2)
     h = x - 6*x**3 + 8*x**4 - 3*x**5
  case(3)
     h = -4*x**3 + 7*x**4 - 3*x**5
  case(4)
     h = HALF * (x**2 - 3*x**3 + 3*x**4 - x**5)
  case(5)
     h = HALF * (x**3 - 2*x**4 + x**5)
  case default
     h = 0
  end select

end function QuinticHermitePolynomial

!-------------------------------------------------------------------------------
!> Quintic Hermite interpolant for scalar coefficients and single argument

pure real(RNP) function QuinticHermiteInterpolant00(c, x) result(h)
  real(RNP), intent(in) :: c(0:5) !< coefficients
  real(RNP), intent(in) :: x      !< argument

  h = c(0) * (1 - 10*x**3 + 15*x**4 - 6*x**5)   &
    + c(1) * (10*x**3 - 15*x**4 + 6*x**5)       &
    + c(2) * (x - 6*x**3 + 8*x**4 - 3*x**5)     &
    + c(3) * (-4*x**3 + 7*x**4 - 3*x**5)        &
    + c(4) * (x**2 - 3*x**3 + 3*x**4 - x**5)/2  &
    + c(5) * (x**3 - 2*x**4 + x**5)/2

end function QuinticHermiteInterpolant00

!-------------------------------------------------------------------------------
!> Quintic Hermite interpolant for scalar coefficients and multiple arguments

pure function QuinticHermiteInterpolant01(c, x) result(h)
  real(RNP), intent(in) :: c(0:5) !< coefficients
  real(RNP), intent(in) :: x(:)   !< arguments
  real(RNP), dimension(size(x)) :: h

  h = c(0) * (1 - 10*x**3 + 15*x**4 - 6*x**5)   &
    + c(1) * (10*x**3 - 15*x**4 + 6*x**5)       &
    + c(2) * (x - 6*x**3 + 8*x**4 - 3*x**5)     &
    + c(3) * (-4*x**3 + 7*x**4 - 3*x**5)        &
    + c(4) * (x**2 - 3*x**3 + 3*x**4 - x**5)/2  &
    + c(5) * (x**3 - 2*x**4 + x**5)/2

end function QuinticHermiteInterpolant01

!===============================================================================

end module Hermite_Interpolation
