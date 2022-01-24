!> summary:  3d generic curvilinear element convection operator for INS (D)
!> author:   Joerg Stiller
!> date:     2022/01/21
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__INS_Convection__3D_D__Gen_RDP
  use Kind_Parameters, only: RWP => RDP
  implicit none
  private

  public :: TPO_INS_Convection_D_Gen

  interface TPO_INS_Convection_D_Gen
    module procedure TPO_INS_Convection_D_Gen_RWP
  end interface

contains

#include "tpo__ins_convection__3d_d__gen_kernel.f"

end module TPO__INS_Convection__3D_D__Gen_RDP
