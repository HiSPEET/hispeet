!> summary:   Manually optimized CtxIxI routines for Intel compilers, NC = 2..7
!> author:    Joerg Stiller, Erik Pfister
!> date:      2020/05/07
!> license:   Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================
!
! available variants
! * CtxIxI_01_jj_01.f


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
