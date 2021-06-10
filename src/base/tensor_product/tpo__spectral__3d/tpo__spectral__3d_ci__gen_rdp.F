!> summary:  3d generic isotropic spectral operator
!> author:   Joerg Stiller
!> date:     2020/05/28
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Spectral__3D_CI__Gen_RDP
  use Kind_Parameters, only: RWP => RDP
  implicit none
  private

  public :: TPO_Spectral_CI_Gen

  interface TPO_Spectral_CI_Gen
    module procedure TPO_Spectral_CI_Gen_RWP
  end interface

contains

#include "tpo__spectral__3d_ci__gen_kernel.f"

end module TPO__Spectral__3D_CI__Gen_RDP
