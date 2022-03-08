!> summary:  DG Schwarz operator for elliptic equations: double precision
!> author:   Joerg Stiller
!> date:     2022/02/20
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================


module DG__Schwarz_Operator__3D__RDP
  use Kind_Parameters
  use Constants
  use Schwarz_Weighting
  use Mesh__3D
  use Mesh_Element__3D
  use Mesh_Element_Indexing__3D
  use Element_Transfer_Buffer__3D
  use DG__Element_Operators__1D
  use DG__Schwarz_Operator__3D__Base
  implicit none
  private

  integer, parameter :: RWP = RDP

#define  DG_SchwarzOperator_3D_RWP DG_SchwarzOperator_3D_RDP
#include "dg__schwarz_operator__3d__rwp.F"
#undef   DG_SchwarzOperator_3D_RWP

end module DG__Schwarz_Operator__3D__RDP
