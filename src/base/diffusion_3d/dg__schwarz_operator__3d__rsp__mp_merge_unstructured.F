!> summary:  Merge subdomain contributions into unstructured mesh variable
!> author:   Joerg Stiller
!> date:     2022/03/07
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule (DG__Schwarz_Operator__3D) MP_Merge_Unstructured_RSP

  integer, parameter :: RWP = RSP

contains

#define  Merge_Unstructured Merge_Unstructured_RSP
#include "dg__schwarz_operator__3d__rwp__mp_merge_unstructured.F"
#undef   Merge_Unstructured

end submodule MP_Merge_Unstructured_RSP
