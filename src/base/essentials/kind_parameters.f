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

!> summary:  Definition of intrinsic type kind parameters
!> author:   Joerg Stiller
!> date:     2014/03/12
!>
!>### Definition of kind parameters for intrinsic types
!===============================================================================

module Kind_Parameters
  use, intrinsic :: ISO_Fortran_Env
  implicit none
  public

  !-----------------------------------------------------------------------------
  ! integer kinds

  integer, parameter :: IXS = selected_int_kind( 2) !<  8 Bit, range > 10^2
  integer, parameter :: IXL = selected_int_kind(18) !< 64 Bit, range > 10^18

  !-----------------------------------------------------------------------------
  ! real kinds

  integer, parameter :: RSP = kind(1E0) !< real single precision
  integer, parameter :: RDP = kind(1D0) !< real double precision

  !> real high precision: 128 bit, if available
  integer, parameter :: RHP = merge(REAL128, RDP, REAL128 > 0)

  !> real normal precision
  integer, parameter :: RNP = RDP

  !=============================================================================

end module Kind_Parameters
