!> summary:  3D generic isotropic averaging operator for real(RDP)
!> author:   Joerg Stiller
!> date:     2022/10/02
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Average__3D_I__Gen_RDP
  use Kind_Parameters, only: RWP => RDP
  implicit none
  private

  public :: TPO_Average_I_Gen

  interface TPO_Average_I_Gen
    module procedure TPO_Average_I_Gen_RWP
  end interface

contains

#include "tpo__average__3d_i__gen_kernel.f"

end module TPO__Average__3D_I__Gen_RDP

