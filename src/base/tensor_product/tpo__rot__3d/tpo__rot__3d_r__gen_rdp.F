!> summary:  Generic rotation of a vector field: 3D Cartesian equidistant
!> author:    Joerg Stiller
!> date:      2021/06/09
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, German
!===============================================================================

module TPO__Rot__3D_R__Gen_RDP
  use Kind_Parameters, only: RWP => RDP
  implicit none
  private

  public :: TPO_Rot_R_Gen

  interface TPO_Rot_R_Gen
    module procedure TPO_Rot_R_Gen_RWP
  end interface

contains

#include "tpo__rot__3d_r__gen_kernel.f"

end module TPO__Rot__3D_R__Gen_RDP
