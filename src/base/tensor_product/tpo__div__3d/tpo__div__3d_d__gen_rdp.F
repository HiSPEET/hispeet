!> summary:  Generic divergence of a vector field: 3D Cartesian curvilinear
!> author:   Jerome Michel, Joerg Stiller
!> date:     2021/08/02
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Div__3D_D__Gen_RDP
  use Kind_Parameters, only: RWP => RDP
  implicit none
  private

  public :: TPO_Div_D_Gen

  interface TPO_Div_D_Gen
    module procedure TPO_Div_D_Gen_RWP
  end interface

contains

#include "tpo__div__3d_d__gen_kernel.f"

end module TPO__Div__3D_D__Gen_RDP
