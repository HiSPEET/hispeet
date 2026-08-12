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

!> summary:   Manually optimized IxIxQt routines for Intel compilers, NQ = 2..33
!> author:    Joerg Stiller, Erik Pfister
!> date:      2019/10/28
!===============================================================================
!
! available variants
! * IxIxQt_01_01_01_01.f
! * IxIxQt_01_01_41_44.f
! * IxIxQt_01_01_41_44_array.f
! * IxIxQt_01_01_44_44.f
! * IxIxQt_01_01_88_01.f

#define _NQ_  2
#include "IxIxQt_01_01_01_01.f"
#undef _NQ_
#define _NQ_  3
#include "IxIxQt_01_01_44_44.f"
#undef _NQ_
#define _NQ_  4
#include "IxIxQt_01_01_01_01.f"
#undef _NQ_
#define _NQ_  5
#include "IxIxQt_01_01_01_01.f"
#undef _NQ_
#define _NQ_  6
#include "IxIxQt_01_01_41_44.f"
#undef _NQ_
#define _NQ_  7
#include "IxIxQt_01_01_41_44.f"
#undef _NQ_
#define _NQ_  8
#include "IxIxQt_01_01_41_44_array.f"
#undef _NQ_
#define _NQ_  9
#include "IxIxQt_01_01_41_44_array.f"
#undef _NQ_
#define _NQ_ 10
#include "IxIxQt_01_01_41_44_array.f"
#undef _NQ_
#define _NQ_ 11
#include "IxIxQt_01_01_41_44_array.f"
#undef _NQ_
#define _NQ_ 12
#include "IxIxQt_01_01_41_44_array.f"
#undef _NQ_
#define _NQ_ 13
#include "IxIxQt_01_01_41_44_array.f"
#undef _NQ_
#define _NQ_ 14
#include "IxIxQt_01_01_41_44_array.f"
#undef _NQ_
#define _NQ_ 15
#include "IxIxQt_01_01_41_44_array.f"
#undef _NQ_
#define _NQ_ 16
#include "IxIxQt_01_01_41_44_array.f"
#undef _NQ_
#define _NQ_ 17
#include "IxIxQt_01_01_41_44_array.f"
#undef _NQ_
#define _NQ_ 18
#include "IxIxQt_01_01_41_44_array.f"
#undef _NQ_
#define _NQ_ 19
#include "IxIxQt_01_01_41_44_array.f"
#undef _NQ_
#define _NQ_ 20
#include "IxIxQt_01_01_41_44_array.f"
#undef _NQ_
#define _NQ_ 21
#include "IxIxQt_01_01_41_44_array.f"
#undef _NQ_
#define _NQ_ 22
#include "IxIxQt_01_01_41_44_array.f"
#undef _NQ_
#define _NQ_ 23
#include "IxIxQt_01_01_41_44_array.f"
#undef _NQ_
#define _NQ_ 24
#include "IxIxQt_01_01_41_44_array.f"
#undef _NQ_
#define _NQ_ 25
#include "IxIxQt_01_01_41_44_array.f"
#undef _NQ_
#define _NQ_ 26
#include "IxIxQt_01_01_41_44_array.f"
#undef _NQ_
#define _NQ_ 27
#include "IxIxQt_01_01_41_44_array.f"
#undef _NQ_
#define _NQ_ 28
#include "IxIxQt_01_01_41_44_array.f"
#undef _NQ_
#define _NQ_ 29
#include "IxIxQt_01_01_41_44_array.f"
#undef _NQ_
#define _NQ_ 30
#include "IxIxQt_01_01_41_44_array.f"
#undef _NQ_
#define _NQ_ 31
#include "IxIxQt_01_01_41_44_array.f"
#undef _NQ_
#define _NQ_ 32
#include "IxIxQt_01_01_41_44_array.f"
#undef _NQ_
#define _NQ_ 33
#include "IxIxQt_01_01_41_44_array.f"
#undef _NQ_
