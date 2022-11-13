!> summary:  3d curvilinear element stress operator with LIBXSMM (DLCI)
!> author:   Joerg Stiller
!> date:     2022/02/04
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Stress__3D_DLCI__XSMM_RDP

#ifdef __LIBXSMM__

  use Kind_Parameters, only: RWP => RDP
  use LIBXSMM
  implicit none
  private

  public :: TPO_Stress_DLCI_XSMM

  interface TPO_Stress_DLCI_XSMM
    module procedure TPO_Stress_DLCI_XSMM_RWP
  end interface

contains

#ifdef __INS_PENALTY__
#include "tpo__stress__3d_dlci__xsmm_kernel-ins_penalty.f"
#else
#include "tpo__stress__3d_dlci__xsmm_kernel.f"
#endif

#endif

end module TPO__Stress__3D_DLCI__XSMM_RDP
