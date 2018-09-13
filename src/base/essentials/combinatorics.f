!> summary:  Routines for combinatorics
!> author:   Joerg Stiller, Immo Huismann
!> date:     2014/05/09
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Routines for combinatorics
!>
!> Currently provides an implementation for the factorial and for binomials.
!===============================================================================

module Combinatorics
  use Kind_Parameters, only: RNP, RHP
  implicit none
  private

  public ::  Factorial
  public ::  Binomial

contains

!-------------------------------------------------------------------------------
!> Computes the factorial of n

pure elemental real(RNP) function Factorial(n) result(fac)
  integer, intent(in) ::  n  !< natural number

  ! first few factorials
  integer, parameter ::  factorials(0:12) =  &
    [ 1, 1, 2, 6, 24, 120, 720, 5040, 40320, &
      362880, 3628800, 39916800, 479001600   ]

  ! select either a known value, or use the Gamma function to compute the value
  select case(n)
  case( lbound(factorials,1) : ubound(factorials,1) )
    fac = factorials(n)
  case default
    fac = real(gamma(real(n+1,RHP)), RNP)
  end select

end function Factorial

!-------------------------------------------------------------------------------
!> Returns the binomial (m,n) = m! / (n!(m-n)!)

pure elemental real(RNP) function Binomial(m, n) result(b)
  integer, intent(in) ::  m, n

  integer ::  i, k, l

  if (0 <= n .and. n <= m) then
    k = min(n, m-n)
    l = max(n, m-n)
    b = 1
    do i = l+1, m
      b = b * i
    end do
    b = b / Factorial(k)
  else
    b = 0
  end if

end function Binomial

!===============================================================================

end module Combinatorics
