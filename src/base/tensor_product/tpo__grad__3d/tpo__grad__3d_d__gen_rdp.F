!> summary:  Generic gradient of a vector field: 3D Cartesian equidistant
!> author:   Jerome Michel, Joerg Stiller
!> date:     2017/01/19
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Grad__3D_D__Gen_RDP
  use Kind_Parameters, only: RWP => RDP
  implicit none
  private

  public :: TPO_Grad_D_Gen

  interface TPO_Grad_D_Gen
    module procedure TPO_Grad_D_Gen_RWP
  end interface

contains

#include "tpo__grad__3d_d__gen_kernel.f"

end module TPO__Grad__3D_D__Gen_RDP
