!> summary:  3d generic element elliptic operator (RLVI)
!> author:   Karl Schoppmann, Joerg Stiller
!> date:     2018/11/20
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Elliptic__3D_RLVI__Gen_RDP
  use Kind_Parameters, only: RWP => RDP
  implicit none
  private

  public :: TPO_Elliptic_RLVI_Gen

  interface TPO_Elliptic_RLVI_Gen
    module procedure TPO_Elliptic_RLVI_Gen_RWP
  end interface

contains

#include "tpo__elliptic__3d_rlvi__gen_kernel.f"

end module TPO__Elliptic__3D_RLVI__Gen_RDP
