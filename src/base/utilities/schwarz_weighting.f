!> author:   Joerg Stiller
!> date:     2016/04/06
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### One-dimensional weighting methods for Schwarz methods
!>
!> The following weighting methods are supported
!>
!>    *  0:  arithmetic average
!>    *  1:  linear transition
!>    *  3:  cubic transition
!>    *  5:  quintic transition
!>    *  7:  7th order transition
!>    *  9:  jump (top hat)
!>    *  default: w = 1, i.e. unweighted additive Schwarz
!>
!===============================================================================

module Schwarz_Weighting
  use Kind_Parameters, only: RNP
  use Constants,       only: ZERO, ONE, HALF
  implicit none
  private

  public :: WeightDistribution

  interface WeightDistribution
    module procedure WeightDistribution_V
    module procedure WeightDistribution_C
  end interface

contains

!-------------------------------------------------------------------------------
!> Weight distribution for varying overlap

pure subroutine WeightDistribution_V(xi, no, method, w)
  real(RNP), intent(in)  :: xi(0:)  !< collocation points (0:po)
  integer,   intent(in)  :: no(2)   !< left/right overlap
  integer,   intent(in)  :: method  !< weighting method {-1,0,1,3,5,7,9}
  real(RNP), intent(out) :: w(:)    !< weight distribution

  integer   :: i, j, o, np
  real(RNP) :: delta(2)

  ! number of element points
  np = size(xi)

  ! return w = -1 if compatibility check fails
  if (size(w) /= np + sum(no)) then
    w = -1
    return
  end if

  ! left and right overlap width in standard coordinates
  delta(1) = OverlapWidth(xi, no(1))
  delta(2) = OverlapWidth(xi, no(2))

  ! initialize subdomain point index
  i = 0

  ! preceding element
  o =  np - no(1) - 1  ! element point offset
  do j = 1, no(1)
    i = i + 1
    w(i) = Weight(xi(o+i), -1, method, delta)
  end do

  ! this element
  o = o - np
  do j = 1, np
    i = i + 1
    w(i) = Weight(xi(o+i), 0, method, delta)
  end do

  ! following element
  o = o - np
  do j = 1, no(2)
    i = i + 1
    w(i) = Weight(xi(o+i), 1, method, delta)
  end do

end subroutine WeightDistribution_V

!-------------------------------------------------------------------------------
!> Weight distribution for constant overlap

pure subroutine WeightDistribution_C(xi, no, method, w)
  real(RNP), intent(in)  :: xi(0:)  !< collocation points (0:po)
  integer,   intent(in)  :: no      !< left/right overlap
  integer,   intent(in)  :: method  !< weighting method {-1,0,1,3,5,7,9}
  real(RNP), intent(out) :: w(:)    !< weight distribution

  integer :: no_v(2)

  no_v = no

  call WeightDistribution_V(xi, no_v, method, w)

end subroutine WeightDistribution_C

!-------------------------------------------------------------------------------
!> Overlap width in standard coordinates

pure real(RNP) function OverlapWidth(xi, no) result(delta)
  real(RNP), intent(in) :: xi(0:)  !< collocation points (0:po)
  integer,   intent(in) :: no      !< overlap

  if (no < 0) then
    delta = 0
  else if (no < size(xi)) then
    delta = 1 + xi(no)
  else
    delta = 2
  end if

end function OverlapWidth

!-------------------------------------------------------------------------------
!> Weighting function

pure real(RNP) function Weight(x, s, method, delta) result(w)
  real(RNP), intent(in) :: x  !< position in standard interval [-1,1]
  integer,   intent(in) :: s  !< element shift: -1|0|1
  integer,   intent(in) :: method    !< index of transition function
  real(RNP), intent(in) :: delta(2)  !< left/right overlap width

  real(RNP) :: xi

  ! default
  w = 1

  select case(method)

  case(0) ! aritmetic average

    if ( x < -(1 + epsilon(x)) + delta(1) )  w = w + 1
    if ( x >  (1 + epsilon(x)) - delta(2) )  w = w + 1
    w = 1 / w

  case(1,3,5,7) ! smooth transition

    select case(s)
    case(0)
      xi = x
    case default
      xi = x + 2*s
    end select

    if (delta(1) > ZERO)  w = w - Psi(-(xi+1)/delta(1), method)
    if (delta(2) > ZERO)  w = w - Psi( (xi-1)/delta(2), method)

  case(9:) ! top-hat distribution

    if (abs(s) > 0)  w = 0

  end select

end function Weight

!-------------------------------------------------------------------------------
!> Transition function for the overlap region

pure real(RNP) function Psi(x, n)
  real(RNP), intent(in) :: x  !< reference coordinate
  integer,   intent(in) :: n  !< exponent, 1|3|5|7 : linear|cubic|quintic|7th

  if (abs(x) < 1) then
    Psi = HALF + HALF * Phi(x, n)
  else if (x < 0) then
    Psi = ZERO
  else
    Psi = ONE
  end if

end function Psi

!-------------------------------------------------------------------------------
!> Shape function for the overlap region

pure real(RNP) function Phi(x, n)
  real(RNP), intent(in) :: x  !< reference coordinate
  integer,   intent(in) :: n  !< exponent, 1|3|5|7 : linear|cubic|quintic|7th

  select case(n)
  case(1) ! linear
    Phi = x
  case(3) ! cubic
    Phi = (3*x - x**3) / 2
  case(5) ! quintic
    Phi = (15*x - 10*x**3 + 3*x**5) / 8
  case(7) ! degree of 7
    Phi = (35*x - 35*x**3 + 21*x**5 - 5*x**7) / 16
  case default
    Phi = 0
  end select

end function Phi

!===============================================================================

end module Schwarz_Weighting
