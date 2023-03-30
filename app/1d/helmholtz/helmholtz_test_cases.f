!> summary:  Test cases for the 1D Helmholtz equation
!> author:   Joerg Stiller
!> date:     2019/01/23
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Test cases for the 1D Helmholtz equation
!>
!> Available test cases
!>
!>     1)  u = sin(PI x)
!>     2)  u = abs(x)^3
!>
!===============================================================================

module Helmholtz_Test_Cases
  use Kind_Parameters, only: RNP
  use Constants,       only: PI
  implicit none
  private

  public :: SetTestCase
  public :: u_exact
  public :: du_exact
  public :: ddu_exact

  integer :: test = 1

contains

!-------------------------------------------------------------------------------
!> Set test case

subroutine SetTestCase(test_case)
  integer, intent(in) :: test_case ! test chosen {1,2}

  test = max(1, min(2, test_case))

end subroutine SetTestCase

!-------------------------------------------------------------------------------
!> Exact solution

elemental real(RNP) function u_exact(x) result(u)
  real(RNP), intent(in) :: x !< point coordinate

  select case(test)
  case(1)
    u = sin(PI*x)
  case(2)
    u = abs(x)**3
  case default
    u = 0
  end select

end function u_exact

!-------------------------------------------------------------------------------
!> Exact first derivative

elemental real(RNP) function du_exact(x) result(du)
  real(RNP), intent(in) :: x !< point coordinate

  select case(test)
  case(1)
    du = PI * cos(PI*x)
  case(2)
    du = 3 * abs(x) * x
  case default
    du = 0
  end select

end function du_exact

!-------------------------------------------------------------------------------
!> Exact second derivative

elemental real(RNP) function ddu_exact(x) result(ddu)
  real(RNP), intent(in) :: x !< point coordinate

  select case(test)
  case(1)
    ddu = -PI**2 * sin(PI*x)
  case(2)
    ddu = 6 * abs(x)
  case default
    ddu = 0
  end select

end function ddu_exact

!===============================================================================

end module Helmholtz_Test_Cases
