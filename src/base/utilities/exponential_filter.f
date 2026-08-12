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

!> summary:  Exponential filter
!> author:   Joerg Stiller
!> date:     2025/03/18
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
    real(RNP), intent(in) :: p     !< order of the filter ≥ 1

    real(RNP), parameter :: alpha = log(epsilon(ONE))

    if (abs(theta) < epsilon(ONE)) then
      sigma = ONE
    else
      sigma = exp(alpha * abs(theta)**p)
    end if

  end function ExponentialFilter

  !=============================================================================

end module Exponential_Filter
