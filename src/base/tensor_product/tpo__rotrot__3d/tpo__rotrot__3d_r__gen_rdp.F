!> summary:  Generic rotrot of a vector field: 3D Cartesian equidistant
!> author:   Joerg Stiller, Erik Pfister
!> date:     2017/01/19
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__RotRot__3D_R__Gen_RWP
  use Kind_Parameters, only: RWP => RDP
  implicit none
  private

  public :: TPO_RotRot_R_Gen

  interface TPO_RotRot_R_Gen
    module procedure TPO_RotRot_R_Gen_RWP
  end interface

contains

#include "tpo__rotrot__3d_r__gen_kernel.f"

end module TPO__RotRot__3D_R__Gen_RWP
