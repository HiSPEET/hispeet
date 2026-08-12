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

!> summary:  Definition of common constants
!> author:   Joerg Stiller
!> date:     2013/10/17
!>
!>### Definition of common constants
!===============================================================================

module Constants
  use Kind_Parameters
  implicit none
  private

  !-----------------------------------------------------------------------------
  ! numeric constants

  real(RNP), parameter, public ::  ZERO  = 0
  real(RNP), parameter, public ::  ONE   = 1
  real(RNP), parameter, public ::  TWO   = 2
  real(RNP), parameter, public ::  THREE = 3
  real(RNP), parameter, public ::  FOUR  = 4
  real(RNP), parameter, public ::  HALF  = 0.5_RNP
  real(RNP), parameter, public ::  THIRD = ONE/3

  real(RHP), parameter, public ::  PI = 3.1415926535897932384626433832795029_RHP

  !==============================================================================

end module Constants
