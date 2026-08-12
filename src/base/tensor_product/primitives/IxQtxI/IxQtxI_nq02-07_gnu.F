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

!> summary:   Manually optimized IxQtxI routines for GNU compilers, NQ = 2..7
!> author:    Joerg Stiller, Erik Pfister
!> date:      2019/10/28
!===============================================================================
!
! available variants
! * IxQtxI_01_01_01_01.f
! * IxQtxI_01_04_01_01.f
! * IxQtxI_01_04_01_41.f
! * IxQtxI_01_44_01_44.f

#define _NQ_  2
#include "IxQtxI_01_44_01_44.f"
#undef _NQ_
#define _NQ_  3
#include "IxQtxI_01_44_01_44.f"
#undef _NQ_
#define _NQ_  4
#include "IxQtxI_01_44_01_44.f"
#undef _NQ_
#define _NQ_  5
#include "IxQtxI_01_01_01_01.f"
#undef _NQ_
#define _NQ_  6
#include "IxQtxI_01_44_01_44.f"
#undef _NQ_
#define _NQ_  7
#include "IxQtxI_01_44_01_44.f"
#undef _NQ_

