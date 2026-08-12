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

!> summary: 3D curvilinear elliptic operator with variable isotropic diffusivity
!> author:  Joerg Stiller
!> date:    2024/08/12
!===============================================================================

module TPO__Elliptic__3D_DLVI
  use Kind_Parameters, only: RDP
  use TPO__Elliptic__3D_DLVI__Hand
  use TPO__Elliptic__3D_DLVI__XSMM
  implicit none
  private

  public :: TPO_Elliptic_DLVI

  interface TPO_Elliptic_DLVI
    module procedure TPO_Elliptic_DLVI_RDP
  end interface

contains

  !-----------------------------------------------------------------------------
  !> 3D element elliptic operator for scalars
  !>
  !> The procedure evaluates the elliptic operator with constant lambda and
  !> variable nu for deformed 3D elements. Given the scalar variable u, the
  !> result is computed as
  !>   .......
  !> for each point (i,j,k) in every element e.
  !>
  !> Ms and Ds are the 1D standard mass and differention matrices for polynomial
  !> order P = np-1, respectively

  subroutine TPO_Elliptic_DLVI_RDP( Ms, Ds, Jd, G, lambda &
                                  , nu,  u,  v, Ji_n      &
                                  , nub, ub, qb           )

    real(RDP), intent(in)  :: Ms(:)        !< 1D standard mass matrix       (np)
    real(RDP), intent(in)  :: Ds(:,:)      !< 1D standard diff matrix    (np,np)
    real(RDP), intent(in)  :: Jd(:,:,:,:)  !< Jacobian determinant (np,np,np,ne)
    real(RDP), intent(in)  :: G(:,:,:,:,:) !< Laplacian metrics
    real(RDP), intent(in)  :: lambda       !< Helmholtz parameter
    real(RDP), intent(in)  :: nu(:,:,:,:)  !< diffusivity          (np,np,np,ne)
    real(RDP), intent(in)  :: u(:,:,:,:)   !< operand              (np,np,np,ne)
    real(RDP), intent(out) :: v(:,:,:,:)   !< result               (np,np,np,ne)

    real(RDP), optional, intent(in)    :: Ji_n(:,:,:,:,:) !< J⁻¹⋅n @ elem faces
    real(RDP), optional, intent(inout) :: nub(:,:,:,:) !< ν at element boundary
    real(RDP), optional, intent(inout) :: ub(:,:,:,:) !< element boundary values
    real(RDP), optional, intent(inout) :: qb(:,:,:,:) !< element boundary fluxes

    contiguous :: Ms, Ds, Jd, G, u, v, Ji_n, nub, ub, qb

    select case(size(Ms))

#ifdef __LIBXSMM__
    case(5:)
      call TPO_Elliptic_DLVI_XSMM(Ms, Ds, Jd, G, lambda, nu, u, v, &
                                  Ji_n, nub, ub, qb)
#endif

    case default
      call TPO_Elliptic_DLVI_Hand(Ms, Ds, Jd, G, lambda, nu, u, v, &
                                  Ji_n, nub, ub, qb)
    end select

  end subroutine TPO_Elliptic_DLVI_RDP

  !=============================================================================

end module TPO__Elliptic__3D_DLVI
