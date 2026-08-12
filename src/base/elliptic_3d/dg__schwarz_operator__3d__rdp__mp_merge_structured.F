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

!> summary:  Merge subdomain contributions into structured mesh variable
!> author:   Joerg Stiller
!> date:     2022/03/07
!===============================================================================

submodule (DG__Schwarz_Operator__3D) MP_Merge_Structured_RDP

  integer, parameter :: RWP = RDP

contains

#define  Merge_Structured Merge_Structured_RDP
#include "dg__schwarz_operator__3d__rwp__mp_merge_structured.F"
#undef   Merge_Structured

end submodule MP_Merge_Structured_RDP
