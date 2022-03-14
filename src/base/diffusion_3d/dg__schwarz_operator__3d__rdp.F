!> summary:  DG Schwarz operator for elliptic equations: double precision
!> author:   Joerg Stiller
!> date:     2022/02/20
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================


submodule (DG__Schwarz_Operator__3D) DG__Schwarz_Operator__3D__RDP
  use Mesh_Element__3D
  use Mesh_Element_Indexing__3D

  integer, parameter :: RWP = RDP

#define  RestrictResidual_RWP RestrictResidual_RDP
#define  MergeCorrections_RWP MergeCorrections_RDP

#include "dg__schwarz_operator__3d__rwp.F"

#undef RestrictResidual_RWP
#undef MergeCorrections_RWP

end submodule DG__Schwarz_Operator__3D__RDP
