!> summary:  3d generic curvilinear element stress operator (DLVI)
!> author:   Joerg Stiller
!> date:     2024/09/30
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Stress__3D_DLVI__Gen_RDP
  use Kind_Parameters, only: RWP => RDP
  implicit none
  private

  public :: TPO_Stress_DLVI_Gen

  interface TPO_Stress_DLVI_Gen
    module procedure TPO_Stress_DLVI_Gen_RWP
  end interface

contains

#include "tpo__stress__3d_dlvi__gen_kernel.f"

end module TPO__Stress__3D_DLVI__Gen_RDP
