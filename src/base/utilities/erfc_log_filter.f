!> summary:  Erfc-Log filter of Boyd
!> author:   Joerg Stiller
!> date:     2025/03/17
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> References
!>   1.  John P Boyd: The Erfc-Log Filter and the Asymptotics of the Euler and
!>       Vandeven Sequence Accelerations, In Ilin & Ridgway Scott: Proc Third
!>       Int Conf Spectral & High Order Meth, Houston J Math 267-276 (1996)
!===============================================================================

module Erfc_Log_Filter
  use Kind_Parameters, only: RNP
  use Constants,       only: ZERO, ONE, HALF
  implicit none
  private

  public :: ErfcLogFilter

contains

  !-----------------------------------------------------------------------------
  !> Erfc-Log filter function of Boyd (1996)

  pure real(RNP) function ErfcLogFilter(theta, p) result(sigma)
    real(RNP), intent(in) :: theta !< non-dimensional wave number in [0,1]
    real(RNP), intent(in) :: p     !< order of the filter ≥ 1

    real(RNP) :: x

    if (abs(theta) < epsilon(ONE)) then
      ! σ(0) = 1
      sigma = ONE
    else if (ONE - abs(theta) < epsilon(ONE)) then
      ! σ(θ) = 0 if |θ| ≥ 1
      sigma = ZERO
    else if (abs(abs(theta) - HALF) < epsilon(ONE)) then
      sigma = HALF
    else
      x = 2 * theta - 1
      sigma = HALF * erfc(x * sqrt( -p/(x*x) * log(1 - x*x) ))
    end if

  end function ErfcLogFilter

  !=============================================================================

end module Erfc_Log_Filter
