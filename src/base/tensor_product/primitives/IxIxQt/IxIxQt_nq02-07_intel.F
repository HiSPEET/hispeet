!> summary:   Manually optimized IxIxQt routines for Intel compilers, NQ = 2..6
!> author:    Joerg Stiller, Erik Pfister
!> date:      2019/10/28
!> license:   Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
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
