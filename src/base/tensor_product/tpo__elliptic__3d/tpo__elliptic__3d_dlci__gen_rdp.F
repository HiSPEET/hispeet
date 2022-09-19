!> summary:  3d curvilinear generic element elliptic operator (DLCI)
!> author:   Jerome Michel, Joerg Stiller
!> date:     2021/07/06
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Elliptic__3D_DLCI__Gen_RDP
  use Kind_Parameters, only: RWP => RDP
  implicit none
  private

  public :: TPO_Elliptic_DLCI_Gen

  interface TPO_Elliptic_DLCI_Gen
    module procedure TPO_Elliptic_DLCI_Gen_RWP
  end interface

contains

#include "tpo__elliptic__3d_dlci__gen_kernel.f"

end module TPO__Elliptic__3D_DLCI__Gen_RDP
