!> summary:  Merge subdomain contributions into structured mesh variable
!> author:   Joerg Stiller
!> date:     2022/03/07
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule (DG__Schwarz_Operator__3D) MP_Merge_Structured_RSP

  integer, parameter :: RWP = RSP

contains

#define  Merge_Structured Merge_Structured_RSP
#include "dg__schwarz_operator__3d__rwp__mp_merge_structured.F"
#undef   Merge_Structured

end submodule MP_Merge_Structured_RSP
