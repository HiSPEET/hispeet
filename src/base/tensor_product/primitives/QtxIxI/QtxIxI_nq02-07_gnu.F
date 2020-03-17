!> summary:   Manually optimized QtxIxI routines for GNU compilers, NQ = 2..7
!> author:    Joerg Stiller, Erik Pfister
!> date:      2019/10/28
!> license:   Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================
!
! available variants
! * QtxIxI_01_01_01_01.f
! * QtxIxI_04_01_01_41.f
! * QtxIxI_44_01_01_44.f

#define _NQ_  2
#include "QtxIxI_44_01_01_44.f"
#undef _NQ_
#define _NQ_  3
#include "QtxIxI_44_01_01_44.f"
#undef _NQ_
#define _NQ_  4
#include "QtxIxI_01_01_01_01.f"
#undef _NQ_
#define _NQ_  5
#include "QtxIxI_04_01_01_41.f"
#undef _NQ_
#define _NQ_  6
#include "QtxIxI_01_01_01_01.f"
#undef _NQ_
#define _NQ_  7
#include "QtxIxI_04_01_01_41.f"
#undef _NQ_

