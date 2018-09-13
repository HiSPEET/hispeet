!> summary:  Cubic and quintic Hermite polynomials and interpolants
!> author:   Joerg Stiller
!> date:     2014/03/12
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
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
  !> Cubic Hermite interpolant h(t) = sum(i=0:3) c(i) H³_i(t)
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
  !> Quintic Hermite interpolant h(t) = sum(i=0:5) c(i) H⁵_i(t)
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
!> Cubic Hermite polynomials H³_i(t), i = 0 ... 3
!>
!> The cubic Hermite polynomials satisfy the conditions
!>
!>     H₀(0) = 1,  H₀(1) = 0,  H₀'(0) = 0,  H₀'(1) = 0
!>     H₁(0) = 0,  H₁(1) = 1,  H₁'(0) = 0,  H₁'(1) = 0
!>     H₂(0) = 0,  H₂(1) = 0,  H₂'(0) = 1,  H₂'(1) = 0
!>     H₃(0) = 0,  H₃(1) = 0,  H₃'(0) = 0,  H₃'(1) = 1

pure real(RNP) function CubicHermitePolynomial(i, t) result(h)
  integer,   intent(in) :: i  !< ID ranging from 0 to 3
  real(RNP), intent(in) :: t  !< 0 <= t <= 1

  select case(i)
  case(0)
     h = 1 - 3*t**2 + 2*t**3
  case(1)
     h = 3*t**2 - 2*t**3
  case(2)
     h = t - 2*t**2 + t**3
  case(3)
     h = -t**2 + t**3
  case default
     h = 0
  end select

end function CubicHermitePolynomial

!-------------------------------------------------------------------------------
!> Cubic Hermite interpolant for scalar coefficients and single argument

pure real(RNP) function CubicHermiteInterpolant00(c, t) result(h)
  real(RNP), intent(in) :: c(0:3) !< coefficients
  real(RNP), intent(in) :: t      !< argument

  h = c(0) * (1 - 3*t**2 + 2*t**3)  &
    + c(1) * (3*t**2 - 2*t**3)      &
    + c(2) * (t - 2*t**2 + t**3)    &
    + c(3) * (-t**2 + t**3)

end function CubicHermiteInterpolant00

!-------------------------------------------------------------------------------
!> Cubic Hermite interpolant for scalar coefficients and multiple arguments

pure function CubicHermiteInterpolant01(c, t) result(h)
  real(RNP), intent(in) :: c(0:3) !< coefficients
  real(RNP), intent(in) :: t(:)   !< arguments

  real(RNP), dimension(size(t)) :: h

  h = c(0) * (1 - 3*t**2 + 2*t**3)  &
    + c(1) * (3*t**2 - 2*t**3)      &
    + c(2) * (t - 2*t**2 + t**3)    &
    + c(3) * (-t**2 + t**3)

end function CubicHermiteInterpolant01

!===============================================================================
! Quintic case

!-------------------------------------------------------------------------------
!> Quintic Hermite polynomials H⁵_i(t), i = 0 ... 5
!>
!> The quintic Hermite polynomials satisfy the conditions
!>
!>     H₀(0) = 1,  H₀(1) = 0,  H₀'(0) = 0,  H₀'(1) = 0,  H₀"(0) = 0,  H₀"(1) = 0
!>     H₁(0) = 0,  H₁(1) = 1,  H₁'(0) = 0,  H₁'(1) = 0,  H₁"(0) = 0,  H₁"(1) = 0
!>     H₂(0) = 0,  H₂(1) = 0,  H₂'(0) = 1,  H₂'(1) = 0,  H₂"(0) = 0,  H₂"(1) = 0
!>     H₃(0) = 0,  H₃(1) = 0,  H₃'(0) = 0,  H₃'(1) = 1,  H₃"(0) = 0,  H₃"(1) = 0
!>     H₄(0) = 0,  H₄(1) = 0,  H₄'(0) = 0,  H₄'(1) = 0,  H₄"(0) = 1,  H₄"(1) = 0
!>     H₅(0) = 0,  H₅(1) = 0,  H₅'(0) = 0,  H₅'(1) = 0,  H₅"(0) = 0,  H₅"(1) = 1

pure real(RNP) function QuinticHermitePolynomial(i, t) result(h)
  integer,   intent(in) :: i  !< ID ranging from 0 to 5
  real(RNP), intent(in) :: t  !< 0 <= t <= 1

  select case(i)
  case(0)
     h = 1 - 10*t**3 + 15*t**4 - 6*t**5
  case(1)
     h = 10*t**3 - 15*t**4 + 6*t**5
  case(2)
     h = t - 6*t**3 + 8*t**4 - 3*t**5
  case(3)
     h = -4*t**3 + 7*t**4 - 3*t**5
  case(4)
     h = HALF * (t**2 - 3*t**3 + 3*t**4 - t**5)
  case(5)
     h = HALF * (t**3 - 2*t**4 + t**5)
  case default
     h = 0
  end select

end function QuinticHermitePolynomial

!-------------------------------------------------------------------------------
!> Quintic Hermite interpolant for scalar coefficients and single argument

pure real(RNP) function QuinticHermiteInterpolant00(c, t) result(h)
  real(RNP), intent(in) :: c(0:5) !< coefficients
  real(RNP), intent(in) :: t      !< argument

  h = c(0) * (1 - 10*t**3 + 15*t**4 - 6*t**5)   &
    + c(1) * (10*t**3 - 15*t**4 + 6*t**5)       &
    + c(2) * (t - 6*t**3 + 8*t**4 - 3*t**5)     &
    + c(3) * (-4*t**3 + 7*t**4 - 3*t**5)        &
    + c(4) * (t**2 - 3*t**3 + 3*t**4 - t**5)/2  &
    + c(5) * (t**3 - 2*t**4 + t**5)/2

end function QuinticHermiteInterpolant00

!-------------------------------------------------------------------------------
!> Quintic Hermite interpolant for scalar coefficients and multiple arguments

pure function QuinticHermiteInterpolant01(c, t) result(h)
  real(RNP), intent(in) :: c(0:5) !< coefficients
  real(RNP), intent(in) :: t(:)   !< arguments
  real(RNP), dimension(size(t)) :: h

  h = c(0) * (1 - 10*t**3 + 15*t**4 - 6*t**5)   &
    + c(1) * (10*t**3 - 15*t**4 + 6*t**5)       &
    + c(2) * (t - 6*t**3 + 8*t**4 - 3*t**5)     &
    + c(3) * (-4*t**3 + 7*t**4 - 3*t**5)        &
    + c(4) * (t**2 - 3*t**3 + 3*t**4 - t**5)/2  &
    + c(5) * (t**3 - 2*t**4 + t**5)/2

end function QuinticHermiteInterpolant01

!===============================================================================

end module Hermite_Interpolation
