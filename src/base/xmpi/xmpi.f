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

!> summary:  Extended Fortran binding to MPI
!> author:   Joerg Stiller
!> date:     2014/11/06, revised 2017/04/10
!===============================================================================

module XMPI

  use MPI_Binding

  use XMPI__Character
  use XMPI__Integer
  use XMPI__Integer_IXS
  use XMPI__Integer_IXL
  use XMPI__Logical
  use XMPI__Real_RSP
  use XMPI__Real_RDP

  logical, private :: initialized = .false.

contains

!-------------------------------------------------------------------------------
!> Initializes the extended MPI Fortran binding

subroutine XMPI_Init()

  if (initialized) then
    return
  else
    call Init_MPI_Binding
    initialized = .true.
  end if

end subroutine XMPI_Init

!===============================================================================

end module XMPI
