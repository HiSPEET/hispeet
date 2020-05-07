!> summary:   Manually optimized IxIxAt routines for Intel compilers, NA = 2..7
!> author:    Joerg Stiller, Erik Pfister
!> date:      2020/05/07
!> license:   Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
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
