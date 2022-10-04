!> summary:  3d generic isotropic Schwarz operator
!> author:   Joerg Stiller
!> date:     2020/06/02
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Schwarz__3D_I__Gen_RSP
  use Kind_Parameters, only: RWP => RSP
  implicit none
  private

  public :: TPO_Schwarz_I_Gen

  interface TPO_Schwarz_I_Gen
    module procedure TPO_Schwarz_I_Gen_RWP
  end interface

contains

#include "tpo__schwarz__3d_i__gen_kernel.f"

end module TPO__Schwarz__3D_I__Gen_RSP
