!------------------------------------------------------------------------------!
! This file is part of HiSPEET: High-order Spectral Element Techniques         !
!                                                                              !
! Copyright (C) 2026 by the HiSPEET authors and the                            !
! Chair of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany                 !
!                                                                              !
! HiSPEET is free software: you can redistribute it and/or modify              !
! it under the terms of the GNU General Public License as published by         !
! the Free Software Foundation, either version 3 of the License, or            !
! (at your option) any later version.                                          !
!                                                                              !
! HiSPEET is distributed in the hope that it will be useful,                   !
! but WITHOUT ANY WARRANTY; without even the implied warranty of               !
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                         !
! See the GNU General Public License for more details.                         !
!                                                                              !
! You should have received a copy of the GNU General Public License            !
! along with HiSPEET. If not, see <http://www.gnu.org/licenses/>.              !
!------------------------------------------------------------------------------!

!> summary:  3d curvilinear element stress operator with LIBXSMM (DLCI)
!> author:   Joerg Stiller
!> date:     2022/02/04
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
