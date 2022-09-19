!> summary:  Merge subdomain core contributions into mesh variable
!> author:   Joerg Stiller
!> date:     2022/03/07
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule (DG__Schwarz_Operator__3D) MP_Merge_Core_RSP

  integer, parameter :: RWP = RSP

contains

#define  Merge_Core Merge_Core_RSP
#include "dg__schwarz_operator__3d__rwp__mp_merge_core.F"
#undef   Merge_Core

end submodule MP_Merge_Core_RSP
