!> summary:   3d generic 8:1 projection operator for real(RDP)
!> author:    Jörg Stiller
!> date:      2024/08/07
!> license:   Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__8To1__3D__Gen_RDP
  use Kind_Parameters, only: RWP => RDP
  implicit none
  private

  public :: TPO_8To1_Gen

  interface TPO_8To1_Gen
    module procedure TPO_8To1_Gen_RWP
  end interface

contains

#include "tpo__8to1__3d__gen_kernel.f"

end module TPO__8To1__3D__Gen_RDP
