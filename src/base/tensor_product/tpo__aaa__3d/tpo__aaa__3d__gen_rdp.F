!> summary:   3d generic AxAxA operator for real(RDP)
!> author:    Jörg Stiller
!> date:      2020/05/07
!> license:   Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__AAA__3D__Gen_RDP
  use Kind_Parameters, only: RWP => RDP
  implicit none
  private

  public :: TPO_AAA_Gen

  interface TPO_AAA_Gen
    module procedure TPO_AAA_Gen_RWP
  end interface

contains

#include "tpo__aaa__3d__gen_kernel.f"

end module TPO__AAA__3D__Gen_RDP
