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

!> summary:   Manually optimized CtxIxI routines for Intel compilers, NC = 2..33
!> author:    Joerg Stiller, Erik Pfister
!> date:      2020/05/07
!===============================================================================
!
! available variants
! CtxIxI_01_jj_01.f
! CtxIxI_02_jj_01.f
! CtxIxI_04_jj_01.f
! CtxIxI_08_jj_01.f
! CtxIxI_04_jj_41.f

! 1:1 operators ----------------------------------------------------------------

#define _NC1_  2
#define _NC2_  2
#include "CtxIxI_01_jj_01.f"
#undef _NC1_
#undef _NC2_

#define _NC1_  3
#define _NC2_  3
#include "CtxIxI_01_jj_01.f"
#undef _NC1_
#undef _NC2_

#define _NC1_  4
#define _NC2_  4
#include "CtxIxI_01_jj_01.f"
#undef _NC1_
#undef _NC2_

#define _NC1_  5
#define _NC2_  5
#include "CtxIxI_01_jj_01.f"
#undef _NC1_
#undef _NC2_

#define _NC1_  6
#define _NC2_  6
#include "CtxIxI_01_jj_01.f"
#undef _NC1_
#undef _NC2_

#define _NC1_  7
#define _NC2_  7
#include "CtxIxI_01_jj_01.f"
#undef _NC1_
#undef _NC2_

#define _NC1_  8
#define _NC2_  8
#include "CtxIxI_01_jj_01.f"
#undef _NC1_
#undef _NC2_

#define _NC1_  9
#define _NC2_  9
#include "CtxIxI_01_jj_01.f"
#undef _NC1_
#undef _NC2_

#define _NC1_  10
#define _NC2_  10
#include "CtxIxI_01_jj_01.f"
#undef _NC1_
#undef _NC2_

#define _NC1_  11
#define _NC2_  11
#include "CtxIxI_01_jj_01.f"
#undef _NC1_
#undef _NC2_

#define _NC1_  12
#define _NC2_  12
#include "CtxIxI_01_jj_01.f"
#undef _NC1_
#undef _NC2_

#define _NC1_  13
#define _NC2_  13
#include "CtxIxI_01_jj_01.f"
#undef _NC1_
#undef _NC2_

#define _NC1_  14
#define _NC2_  14
#include "CtxIxI_01_jj_01.f"
#undef _NC1_
#undef _NC2_

#define _NC1_  15
#define _NC2_  15
#include "CtxIxI_01_jj_01.f"
#undef _NC1_
#undef _NC2_

#define _NC1_  16
#define _NC2_  16
#include "CtxIxI_01_jj_01.f"
#undef _NC1_
#undef _NC2_

#define _NC1_  17
#define _NC2_  17
#include "CtxIxI_01_jj_01.f"
#undef _NC1_
#undef _NC2_

#define _NC1_  18
#define _NC2_  18
#include "CtxIxI_01_jj_01.f"
#undef _NC1_
#undef _NC2_

#define _NC1_  19
#define _NC2_  19
#include "CtxIxI_01_jj_01.f"
#undef _NC1_
#undef _NC2_

#define _NC1_  20
#define _NC2_  20
#include "CtxIxI_01_jj_01.f"
#undef _NC1_
#undef _NC2_

#define _NC1_  21
#define _NC2_  21
#include "CtxIxI_01_jj_01.f"
#undef _NC1_
#undef _NC2_

#define _NC1_  22
#define _NC2_  22
#include "CtxIxI_01_jj_01.f"
#undef _NC1_
#undef _NC2_

#define _NC1_  23
#define _NC2_  23
#include "CtxIxI_01_jj_01.f"
#undef _NC1_
#undef _NC2_

#define _NC1_  24
#define _NC2_  24
#include "CtxIxI_01_jj_01.f"
#undef _NC1_
#undef _NC2_

#define _NC1_  25
#define _NC2_  25
#include "CtxIxI_01_jj_01.f"
#undef _NC1_
#undef _NC2_

#define _NC1_  26
#define _NC2_  26
#include "CtxIxI_01_jj_01.f"
#undef _NC1_
#undef _NC2_

#define _NC1_  27
#define _NC2_  27
#include "CtxIxI_01_jj_01.f"
#undef _NC1_
#undef _NC2_

#define _NC1_  28
#define _NC2_  28
#include "CtxIxI_01_jj_01.f"
#undef _NC1_
#undef _NC2_

#define _NC1_  29
#define _NC2_  29
#include "CtxIxI_01_jj_01.f"
#undef _NC1_
#undef _NC2_

#define _NC1_  30
#define _NC2_  30
#include "CtxIxI_01_jj_01.f"
#undef _NC1_
#undef _NC2_

#define _NC1_  31
#define _NC2_  31
#include "CtxIxI_01_jj_01.f"
#undef _NC1_
#undef _NC2_

#define _NC1_  32
#define _NC2_  32
#include "CtxIxI_01_jj_01.f"
#undef _NC1_
#undef _NC2_

#define _NC1_  33
#define _NC2_  33
#include "CtxIxI_01_jj_01.f"
#undef _NC1_
#undef _NC2_


! 1:2 operators ----------------------------------------------------------------

#define _NC1_  2
#define _NC2_  3
#include "CtxIxI_01_jj_01.f"
#undef _NC1_
#undef _NC2_

#define _NC1_  2
#define _NC2_  4
#include "CtxIxI_01_jj_01.f"
#undef _NC1_
#undef _NC2_

#define _NC1_  3
#define _NC2_  4
#include "CtxIxI_01_jj_01.f"
#undef _NC1_
#undef _NC2_

#define _NC1_  3
#define _NC2_  5
#include "CtxIxI_01_jj_01.f"
#undef _NC1_
#undef _NC2_

#define _NC1_  4
#define _NC2_  8
#include "CtxIxI_01_jj_01.f"
#undef _NC1_
#undef _NC2_

#define _NC1_  5
#define _NC2_  9
#include "CtxIxI_01_jj_01.f"
#undef _NC1_
#undef _NC2_

#define _NC1_  8
#define _NC2_  16
#include "CtxIxI_01_jj_01.f"
#undef _NC1_
#undef _NC2_

#define _NC1_  9
#define _NC2_  17
#include "CtxIxI_01_jj_01.f"
#undef _NC1_
#undef _NC2_

#define _NC1_  16
#define _NC2_  30
#include "CtxIxI_01_jj_01.f"
#undef _NC1_
#undef _NC2_

#define _NC1_  17
#define _NC2_  33
#include "CtxIxI_01_jj_01.f"
#undef _NC1_
#undef _NC2_


! 2:1 operators ----------------------------------------------------------------

#define _NC1_  3
#define _NC2_  2
#include "CtxIxI_01_jj_01.f"
#undef _NC1_
#undef _NC2_

#define _NC1_  4
#define _NC2_  2
#include "CtxIxI_01_jj_01.f"
#undef _NC1_
#undef _NC2_

#define _NC1_  4
#define _NC2_  3
#include "CtxIxI_01_jj_01.f"
#undef _NC1_
#undef _NC2_

#define _NC1_  5
#define _NC2_  3
#include "CtxIxI_01_jj_01.f"
#undef _NC1_
#undef _NC2_

#define _NC1_  8
#define _NC2_  4
#include "CtxIxI_01_jj_01.f"
#undef _NC1_
#undef _NC2_

#define _NC1_  9
#define _NC2_  5
#include "CtxIxI_01_jj_01.f"
#undef _NC1_
#undef _NC2_

#define _NC1_  16
#define _NC2_  8
#include "CtxIxI_01_jj_01.f"
#undef _NC1_
#undef _NC2_

#define _NC1_  17
#define _NC2_  9
#include "CtxIxI_01_jj_01.f"
#undef _NC1_
#undef _NC2_

#define _NC1_  30
#define _NC2_  16
#include "CtxIxI_01_jj_01.f"
#undef _NC1_
#undef _NC2_

#define _NC1_  33
#define _NC2_  17
#include "CtxIxI_01_jj_01.f"
#undef _NC1_
#undef _NC2_
