!> summary:  Generic divergence of a vector field: 3D Cartesian equidistant
!> author:   Joerg Stiller, Erik Pfister
!> date:     2017/01/19
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Div__3D_R__Gen_RDP
  use Kind_Parameters, only: RWP => RDP
  implicit none
  private

  public :: TPO_Div_R_Gen

  interface TPO_Div_R_Gen
    module procedure TPO_Div_R_Gen_RWP
  end interface

contains

#include "tpo__div__3d__gen_kernel.f"

end module TPO__Div__3D_R__Gen_RDP
