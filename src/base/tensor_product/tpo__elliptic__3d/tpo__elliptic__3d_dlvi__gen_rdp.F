!> summary:  3d curvilinear generic element elliptic operator (DLVI)
!> author:   Joerg Stiller
!> date:     2024/08/12
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Elliptic__3D_DLVI__Gen_RDP
  use Kind_Parameters, only: RWP => RDP
  implicit none
  private

  public :: TPO_Elliptic_DLVI_Gen

  interface TPO_Elliptic_DLVI_Gen
    module procedure TPO_Elliptic_DLVI_Gen_RWP
  end interface

contains

#include "tpo__elliptic__3d_dlvi__gen_kernel.f"

end module TPO__Elliptic__3D_DLVI__Gen_RDP
