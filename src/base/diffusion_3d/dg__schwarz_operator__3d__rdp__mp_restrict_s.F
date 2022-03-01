!> summary:  Restriction a structured mesh variable to Schwarz subdomains
!> author:   Joerg Stiller
!> date:     2022/03/01
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule (DG__Schwarz_Operator__3D__RDP) MP_Restrict_S
contains

#define  DG_SchwarzOperator_3D_RWP DG_SchwarzOperator_3D_RDP
#include "dg__schwarz_operator__3d__rwp__mp_restrict_s.F"
#undef   DG_SchwarzOperator_3D_RWP

end submodule MP_Restrict_S
