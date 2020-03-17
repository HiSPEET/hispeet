!> summary:   Manually optimized IxQtxI routines for GNU compilers, NQ = 2..7
!> author:    Joerg Stiller, Erik Pfister
!> date:      2019/10/28
!> license:   Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
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

