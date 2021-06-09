!> summary:  3D generic scaled constant isotropic diagonal operator for real(RDP)
!> author:   Joerg Stiller
!> date:     2020/06/04
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Diagonal__3D_CI__Gen_RDP
  use Kind_Parameters, only: RWP => RDP
  implicit none
  private

  public :: TPO_Diagonal_CI_Gen

  interface TPO_Diagonal_CI_Gen
    module procedure TPO_Diagonal_CI_Gen_RWP
  end interface

contains

#include "tpo__diagonal__3d_ci__gen_kernel.f"

end module TPO__Diagonal__3D_CI__Gen_RDP

