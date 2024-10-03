!> summary:  3d curvilinear element stress operator with LIBXSMM (DLVI)
!> author:   Joerg Stiller
!> date:     2024/09/30
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Stress__3D_DLVI__XSMM_RDP

#ifdef __LIBXSMM__

  use Kind_Parameters, only: RWP => RDP
  use LIBXSMM
  implicit none
  private

  public :: TPO_Stress_DLVI_XSMM

  interface TPO_Stress_DLVI_XSMM
    module procedure TPO_Stress_DLVI_XSMM_RWP
  end interface

contains

#include "tpo__stress__3d_dlvi__xsmm_kernel.f"

#endif

end module TPO__Stress__3D_DLVI__XSMM_RDP
