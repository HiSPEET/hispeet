!> summary:  3d generic anisotropic spectral operator
!> author:   Joerg Stiller
!> date:     2020/06/02
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Spectral__3D_CA__Gen_RDP
  use Kind_Parameters, only: RWP => RDP
  implicit none
  private

  public :: TPO_Spectral_CA_Gen

  interface TPO_Spectral_CA_Gen
    module procedure TPO_Spectral_CA_Gen_RWP
  end interface

contains

#include "tpo__spectral__3d_ca__hand_kernel.f"

end module TPO__Spectral__3D_CA__Gen_RDP
