!> summary:   3d generic element diffusion operator (RLCI) for real(RDP)
!> author:    Jörg Stiller
!> date:      2019/12/06
!> license:   Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Diffusion__3D_RLCI__Gen_RDP
  use Kind_Parameters, only: RWP => RDP
  implicit none
  private

  public :: TPO_Diffusion_RLCI_Gen

  interface TPO_Diffusion_RLCI_Gen
    module procedure TPO_Diffusion_RLCI_Gen_RWP
  end interface

contains

#include "tpo__diffusion__3d_rlci__gen_kernel.f"

end module TPO__Diffusion__3D_RLCI__Gen_RDP
