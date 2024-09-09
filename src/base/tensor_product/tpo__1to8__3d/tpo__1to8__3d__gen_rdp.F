!> summary:   3d generic 1:8 interpolation operator for real(RDP)
!> author:    Jörg Stiller
!> date:      2024/07/29
!> license:   Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__1To8__3D__Gen_RDP
  use Kind_Parameters, only: RWP => RDP
  implicit none
  private

  public :: TPO_1To8_Gen

  interface TPO_1To8_Gen
    module procedure TPO_1To8_Gen_RWP
  end interface

contains

#include "tpo__1to8__3d__gen_kernel-321.f"

end module TPO__1To8__3D__Gen_RDP
