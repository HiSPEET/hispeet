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

!> summary:  3d curvilinear element stress operator (DLVI)
!> author:   Joerg Stiller
!> date:     2024/09/30
!===============================================================================

module TPO__Stress__3D_DLVI
  use Kind_Parameters, only: RDP
  use TPO__Stress__3D_DLVI__Gen
  use TPO__Stress__3D_DLVI__XSMM
  implicit none
  private

  public :: TPO_Stress

  interface TPO_Stress
    module procedure TPO_Stress_DLVI_RDP
  end interface

contains

  !-----------------------------------------------------------------------------
  !> 3D element stress operator -- variable viscosity

  subroutine TPO_Stress_DLVI_RDP(Ms, Ds, Jd, Ji, n, zeta, eta, v, fv, vb, sb)

    real(RDP), contiguous, intent(in) :: Ms(:)
    !< 1D diagonal standard mass matrix (np)

    real(RDP), contiguous, intent(in)  :: Ds(:,:)
    !< 1D standard diff matrix (np,np)

    real(RDP), contiguous, intent(in)  :: Jd(:,:,:,:)
    !< Jacobian determinant (np,np,np,ne)

    real(RDP), contiguous, intent(in)  :: Ji(:,:,:,:,:,:)
    !< inverse Jacobian matrix  (np,np,np,ne,3,3)

    real(RDP), contiguous, intent(in)  :: n(:,:,:,:,:)
    !< unit normal vector at element boundary faces (np,np,6,ne,3)

    real(RDP), contiguous, intent(in)  :: zeta(:,:,:,:)
    !< bulk viscosity or modulus ζ

    real(RDP), contiguous, intent(in)  :: eta(:,:,:,:)
    !< shear viscosity or modulus η

    real(RDP), contiguous, intent(in)  :: v(:,:,:,:,:)
    !< velocity or displacement (np,np,np,ne,3)

    real(RDP), contiguous, intent(out) :: fv(:,:,:,:,:)
    !< element volume contribution to stress operator (np,np,np,ne,3)

    real(RDP), contiguous, optional, intent(out) :: vb(:,:,:,:,:)
    !< element boundary values of velocity or displacement (np,np,6,ne+ng,3),

    real(RDP), contiguous, optional, intent(out) :: sb(:,:,:,:,:)
    !< stress normal to element boundary fluxes, sb = n⁻⋅τ⁻ (np,np,6,ne+ng,3)

    select case(size(Ms))

#ifdef __LIBXSMM__
    case(3:)
      call TPO_Stress_DLVI_XSMM(Ms, Ds, Jd, Ji, n, zeta, eta, v, fv, vb, sb)
#endif

    case default
      call TPO_Stress_DLVI_Gen(Ms, Ds, Jd, Ji, n, zeta, eta, v, fv, vb, sb)
    end select

  end subroutine TPO_Stress_DLVI_RDP

  !=============================================================================

end module TPO__Stress__3D_DLVI
