!> summary:  3d curvilinear generic element diffusion operator (DLCI)
!> author:   Jerome Michel, Joerg Stiller
!> date:     2021/07/06
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Diffusion__3D_DLCI__Gen_RDP
  use Kind_Parameters, only: RWP => RDP
  implicit none
  private

  public :: TPO_Diffusion_DLCI_Gen

  interface TPO_Diffusion_DLCI_Gen
    module procedure TPO_Diffusion_DLCI_Gen_RWP
  end interface

contains

#include "tpo__diffusion__3d_dlci__gen_kernel.f"

end module TPO__Diffusion__3D_DLCI__Gen_RDP
