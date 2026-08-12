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

!> summary:   Manually optimized IxIxAt routines for Intel compilers, NA = 2..7
!> author:    Joerg Stiller, Erik Pfister
!> date:      2020/05/07
!===============================================================================
!
! available variants
! * IxIxAt_01_01_01_01.f


! 1:1 operators ----------------------------------------------------------------

#define _NA1_  2
#define _NA2_  2
#include "IxIxAt_01_01_01_01.f"
#undef _NA1_
#undef _NA2_

#define _NA1_  3
#define _NA2_  3
#include "IxIxAt_01_01_01_01.f"
#undef _NA1_
#undef _NA2_

#define _NA1_  4
#define _NA2_  4
#include "IxIxAt_01_01_01_01.f"
#undef _NA1_
#undef _NA2_

#define _NA1_  5
#define _NA2_  5
#include "IxIxAt_01_01_01_01.f"
#undef _NA1_
#undef _NA2_

#define _NA1_  6
#define _NA2_  6
#include "IxIxAt_01_01_01_01.f"
#undef _NA1_
#undef _NA2_

#define _NA1_  7
#define _NA2_  7
#include "IxIxAt_01_01_01_01.f"
#undef _NA1_
#undef _NA2_


! 1:2 operators ----------------------------------------------------------------

#define _NA1_  2
#define _NA2_  3
#include "IxIxAt_01_01_01_01.f"
#undef _NA1_
#undef _NA2_

#define _NA1_  2
#define _NA2_  4
#include "IxIxAt_01_01_01_01.f"
#undef _NA1_
#undef _NA2_

#define _NA1_  3
#define _NA2_  4
#include "IxIxAt_01_01_01_01.f"
#undef _NA1_
#undef _NA2_

#define _NA1_  3
#define _NA2_  5
#include "IxIxAt_01_01_01_01.f"
#undef _NA1_
#undef _NA2_

! 2:1 operators ----------------------------------------------------------------

#define _NA1_  3
#define _NA2_  2
#include "IxIxAt_01_01_01_01.f"
#undef _NA1_
#undef _NA2_

#define _NA1_  4
#define _NA2_  2
#include "IxIxAt_01_01_01_01.f"
#undef _NA1_
#undef _NA2_

#define _NA1_  4
#define _NA2_  3
#include "IxIxAt_01_01_01_01.f"
#undef _NA1_
#undef _NA2_

#define _NA1_  5
#define _NA2_  3
#include "IxIxAt_01_01_01_01.f"
#undef _NA1_
#undef _NA2_
