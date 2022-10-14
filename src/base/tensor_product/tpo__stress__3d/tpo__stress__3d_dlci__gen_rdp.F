!> summary:  3d generic curvilinear element stress operator (DLCI)
!> author:   Joerg Stiller
!> date:     2022/01/02
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Stress__3D_DLCI__Gen_RDP
  use Kind_Parameters, only: RWP => RDP
  implicit none
  private

  public :: TPO_Stress_DLCI_Gen

  interface TPO_Stress_DLCI_Gen
    module procedure TPO_Stress_DLCI_Gen_RWP
  end interface

contains

#include "tpo__stress__3d_dlci__gen_kernel.f"

end module TPO__Stress__3D_DLCI__Gen_RDP
