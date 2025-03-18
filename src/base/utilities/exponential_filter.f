!> summary:  Exponential filter
!> author:   Joerg Stiller
!> date:     2025/03/18
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Exponential_Filter
  use Kind_Parameters, only: RNP
  use Constants,       only: ONE
  implicit none
  private

  public :: ExponentialFilter

contains

  !-----------------------------------------------------------------------------
  !> Exponential filter function

  pure real(RNP) function ExponentialFilter(theta, p) result(sigma)
    real(RNP), intent(in) :: theta !< non-dimensional wave number in [0,1]
    integer,   intent(in) :: p     !< order of the filter ≥ 1

    real(RNP), parameter :: alpha = log(epsilon(ONE))

    if (abs(theta) < epsilon(ONE)) then
      sigma = ONE
    else
      sigma = exp(alpha * abs(theta)**p)
    end if

  end function ExponentialFilter

  !=============================================================================

end module Exponential_Filter
