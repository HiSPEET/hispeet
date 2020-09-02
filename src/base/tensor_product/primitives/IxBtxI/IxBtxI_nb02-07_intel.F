!> summary:   Manually optimized IxBtxI routines for Intel compilers, NB = 2..7
!> author:    Joerg Stiller, Erik Pfister
!> date:      2020/05/07
!> license:   Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================
!
! available variants
! * IxBtxI_01_01_01_01.f


! 1:1 operators ----------------------------------------------------------------

#define _NB1_  2
#define _NB2_  2
#include "IxBtxI_01_01_01_01.f"
#undef _NB1_
#undef _NB2_

#define _NB1_  3
#define _NB2_  3
#include "IxBtxI_01_01_01_01.f"
#undef _NB1_
#undef _NB2_

#define _NB1_  4
#define _NB2_  4
#include "IxBtxI_01_01_01_01.f"
#undef _NB1_
#undef _NB2_

#define _NB1_  5
#define _NB2_  5
#include "IxBtxI_01_01_01_01.f"
#undef _NB1_
#undef _NB2_

#define _NB1_  6
#define _NB2_  6
#include "IxBtxI_01_01_01_01.f"
#undef _NB1_
#undef _NB2_

#define _NB1_  7
#define _NB2_  7
#include "IxBtxI_01_01_01_01.f"
#undef _NB1_
#undef _NB2_


! 1:2 operators ----------------------------------------------------------------

#define _NB1_  2
#define _NB2_  3
#include "IxBtxI_01_01_01_01.f"
#undef _NB1_
#undef _NB2_

#define _NB1_  2
#define _NB2_  4
#include "IxBtxI_01_01_01_01.f"
#undef _NB1_
#undef _NB2_

#define _NB1_  3
#define _NB2_  4
#include "IxBtxI_01_01_01_01.f"
#undef _NB1_
#undef _NB2_

#define _NB1_  3
#define _NB2_  5
#include "IxBtxI_01_01_01_01.f"
#undef _NB1_
#undef _NB2_

! 2:1 operators ----------------------------------------------------------------

#define _NB1_  3
#define _NB2_  2
#include "IxBtxI_01_01_01_01.f"
#undef _NB1_
#undef _NB2_

#define _NB1_  4
#define _NB2_  2
#include "IxBtxI_01_01_01_01.f"
#undef _NB1_
#undef _NB2_

#define _NB1_  4
#define _NB2_  3
#include "IxBtxI_01_01_01_01.f"
#undef _NB1_
#undef _NB2_

#define _NB1_  5
#define _NB2_  3
#include "IxBtxI_01_01_01_01.f"
#undef _NB1_
#undef _NB2_
