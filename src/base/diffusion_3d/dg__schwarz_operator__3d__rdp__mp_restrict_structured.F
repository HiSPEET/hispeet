!> summary:  Restriction a structured mesh variable to Schwarz subdomains
!> author:   Joerg Stiller
!> date:     2022/03/01
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule (DG__Schwarz_Operator__3D) MP_Restrict_Structured_RDP

  integer, parameter :: RWP = RDP

contains

#define  Restrict_Structured Restrict_Structured_RDP
#include "dg__schwarz_operator__3d__rwp__mp_restrict_structured.F"

end submodule MP_Restrict_Structured_RDP
