!> summary:  3d generic element diffusion operator (RLVI)
!> author:   Karl Schoppmann, Joerg Stiller
!> date:     2018/11/20
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Diffusion__3D_RLVI__Gen_RDP
  use Kind_Parameters, only: RWP => RDP
  implicit none
  private

  public :: TPO_Diffusion_RLVI_Gen

  interface TPO_Diffusion_RLVI_Gen
    module procedure TPO_Diffusion_RLVI_Gen_RWP
  end interface

contains

#include "tpo__diffusion__3d_rlvi__gen_kernel.f"

end module TPO__Diffusion__3D_RLVI__Gen_RDP
