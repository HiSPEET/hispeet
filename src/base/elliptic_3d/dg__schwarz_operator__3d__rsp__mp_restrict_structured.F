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

!> summary:  Restriction a structured mesh variable to Schwarz subdomains
!> author:   Joerg Stiller
!> date:     2022/03/01
!===============================================================================

submodule (DG__Schwarz_Operator__3D) MP_Restrict_Structured_RSP

  integer, parameter :: RWP = RSP

contains

#define  Restrict_Structured Restrict_Structured_RSP
#include "dg__schwarz_operator__3d__rwp__mp_restrict_structured.F"

end submodule MP_Restrict_Structured_RSP
