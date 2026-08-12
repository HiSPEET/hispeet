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

!> summary:  Solution of linear, quadratic and cubic equations
!> author:   Joerg Stiller
!> date:     2016/08/15
!>
!>### Solution of linear, quadratic and cubic equations
!===============================================================================

module Algebraic_Equations
  use Kind_Parameters, only: RNP
  use Constants
  implicit none
  private

  public :: SolveLinearEquation
  public :: SolveQuadraticEquation
  public :: SolveCubicEquation

contains

!-------------------------------------------------------------------------------
!> Solution of a linear equation
!>
!> Determines the solution of
!>
!>       a + b·x = 0
!>
!> such that
!>
!>       xmin <= x <= xmax
!>
!> with optional bounds xmin and xmax.

pure subroutine SolveLinearEquation(a, b, xmin, xmax, x, n)
  real(RNP),           intent(in)  ::  a, b  !< coefficients
  real(RNP),           intent(out) ::  x     !< solution
  integer,             intent(out) ::  n     !< number of solutions
  real(RNP), optional, intent(in)  ::  xmin  !< optional lower bound
  real(RNP), optional, intent(in)  ::  xmax  !< optional upper bound

  real(RNP) ::  s, smin, smax

  x = 0
  n = 0
  if (abs(b) > 0) then
     s = -sign(ONE,a) * sign(ONE,b) * min(abs(a)/abs(b),huge(s)) ! = -a/b
     smin = -huge(s);  if (present(xmin)) smin = xmin
     smax =  huge(s);  if (present(xmax)) smax = xmax
     if (smin <= s .and. s <= smax) then
        x = s
        n = 1
     end if
  end if
end subroutine SolveLinearEquation

!-------------------------------------------------------------------------------
!> Solution of a quadratic equation
!>
!> Determines the real solution(s) of
!>
!>       a + b·x + c·x² = 0
!>
!> such that
!>
!>       xmin <= x <= xmax
!>
!> with optional bounds xmin and xmax.

pure subroutine SolveQuadraticEquation(a, b, c, xmin, xmax, x, n)
  real(RNP), intent(in)  ::  a, b, c  !< coefficients
  real(RNP), intent(out) ::  x(2)     !< solution(s),  x(1) ≥ x(2)
  integer,   intent(out) ::  n        !< number of solutions
  real(RNP), optional, intent(in)  ::  xmin  !< optional lower bound
  real(RNP), optional, intent(in)  ::  xmax  !< optional upper bound

  real(RNP) ::  s(2), smin, smax
  real(RNP) ::  d, p, q
  integer   ::  i, m

  n = 0
  x = 0
  s = 0

  if (c**2 > (a**2 + b**2)/1E12) then
     p = b / (2*c)
     q = a / c
     d = p*p - q
     if (d > 0) then
        d = sqrt(d)
        s = [ d, -d ] - p
        m = 2
     else if (d > -epsilon(x)) then
        s = -p
        m = 1
     else
        s = 0
        m = 0
     end if
  else if (abs(b) > 0) then
     s = -sign(ONE,a) * sign(ONE,b) * min(abs(a)/abs(b),huge(s)) ! = -a/b
     m = 1
  else
     s = 0
     m = 0
  end if

  smin = -huge(s);  if (present(xmin)) smin = xmin
  smax =  huge(s);  if (present(xmax)) smax = xmax
  do i = 1, m
     if (smin <= s(i) .and. s(i) <= smax) then
        n = n + 1
        x(n) = s(i)
     end if
  end do

end subroutine SolveQuadraticEquation

!-------------------------------------------------------------------------------
!> Solution of a cubic equation
!>
!> Determines the real solution(s) of
!>
!>       a + b·x + c·x² + x³ = 0
!>
!> such that
!>
!>       xmin <= x <= xmax
!>
!> with optional bounds xmin and xmax.

pure subroutine SolveCubicEquation(a, b, c, xmin, xmax, x, n)
  real(RNP), intent(in)  ::  a, b, c         !< coefficients
  real(RNP), intent(out) ::  x(3)            !< solution(s)
  integer,   intent(out) ::  n               !< number of solutions
  real(RNP), optional, intent(in)  ::  xmin  !< optional lower bound
  real(RNP), optional, intent(in)  ::  xmax  !< optional upper bound

  real(RNP) ::  s(3), smin, smax
  real(RNP) ::  d, p, q, beta, z
  integer   ::  i, m

  n = 0
  x = 0
  s = 0

  p = (b - (c**2)/3)/3
  q = ((2*c**3)/27 - b*c/3 + a)/2
  d = p**3 + q**2
  if (p < 0) then
     p = sign(ONE,q) * sqrt(-p)
     z = q / p**3
     if (d <= 0) then
        m = 3
        z = min(z, ONE)
        beta =  acos(z) / 3
        s(1) = -2*p*cos(beta)
        s(2) =  2*p*cos(beta + real(PI/3, RNP) )
        s(3) =  2*p*cos(beta - real(PI/3, RNP) )
     else
        m = 1
        beta =  log(z + sqrt(max(z*z - ONE, ZERO))) / 3 ! = acosh(z)/3
        s(1) = -2*p*cosh(beta)
     end if
  else if (p > 0) then
     p = sign(ONE,q) * sqrt(p)
     z = q / p**3
     m = 1
     beta =  log(z + sqrt(z*z + ONE)) / 3
     s(1) = -2*p*sinh(beta)
  else ! p == 0
     m = 1
     if (q <= 0) then
        s(1) = (-2*q)**THIRD
     else
        s(1) = -(2*q)**THIRD
     end if
  end if
  s(1:m) = s(1:m) - c/3

  smin = -huge(s);  if (present(xmin)) smin = xmin
  smax =  huge(s);  if (present(xmax)) smax = xmax
  do i = 1, m
     if (smin <= s(i) .and. s(i) <= smax) then
        n = n + 1
        x(n) = s(i)
     end if
  end do

end subroutine SolveCubicEquation

!===============================================================================

end module Algebraic_Equations
