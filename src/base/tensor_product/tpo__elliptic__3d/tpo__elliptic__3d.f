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

!> summary:  3D element elliptic operators
!> author:   Joerg Stiller
!> date:     2020/05/29
!===============================================================================

module TPO__Elliptic__3D
  use Kind_Parameters, only: RDP
  use Standard_Element_Operators__1D
  use Mesh__3D
  use Spectral_Element_Mesh__3D

  use TPO__Elliptic__3D_RLCI
  use TPO__Elliptic__3D_RLVI
  use TPO__Elliptic__3D_DLCI
  use TPO__Elliptic__3D_DLVI

  private

  public :: TPO_Elliptic

  interface TPO_Elliptic
    module procedure TPO_Elliptic_CI_RDP
    module procedure TPO_Elliptic_VI_RDP
  end interface

contains

  !-----------------------------------------------------------------------------
  !> Automatic elliptic TPO with constant isotropic diffusivity
  !> and double precision

  subroutine TPO_Elliptic_CI_RDP(eop, sem, lambda, nu, u, v, ub, qb)
    class(StandardElementOperators_1D), intent(in) :: eop
    class(SpectralElementMesh_3D), intent(in) :: sem
    real(RDP),           intent(in)    :: lambda      !< Helmholtz parameter λ
    real(RDP),           intent(in)    :: nu          !< diffusivity ν
    real(RDP),           intent(in)    :: u(:,:,:,:)  !< 3D scalar field u
    real(RDP),           intent(out)   :: v(:,:,:,:)  !< result
    real(RDP), optional, intent(inout) :: ub(:,:,:,:) !< element boundary values
    real(RDP), optional, intent(inout) :: qb(:,:,:,:) !< element boundary fluxes

    contiguous :: u, v, ub, qb

    if (sem % mesh % regular) then
      call TPO_Elliptic_RLCI( Ms     = eop % w              &
                            , Ls     = eop % L              &
                            , dx     = sem % mesh % dx      &
                            , lambda = lambda               &
                            , nu     = nu                   &
                            , u      = u                    &
                            , v      = v                    &
                            , Ds     = eop % D              &
                            , ub     = ub                   &
                            , qb     = qb                   )
    else
      call TPO_Elliptic_DLCI( Ms     = eop % w              &
                            , Ds     = eop % D              &
                            , Jd     = sem % metrics % Jd   &
                            , G      = sem % metrics % G    &
                            , lambda = lambda               &
                            , nu     = nu                   &
                            , u      = u                    &
                            , v      = v                    &
                            , Ji_n   = sem % metrics % Ji_n &
                            , ub     = ub                   &
                            , qb     = qb                   )
    end if

  end subroutine TPO_Elliptic_CI_RDP

  !-----------------------------------------------------------------------------
  !> Automatic elliptic TPO with variable isotropic diffusivity
  !> and double precision

  subroutine TPO_Elliptic_VI_RDP(eop, sem, lambda, nu, u, v, nub, ub, qb)
    class(StandardElementOperators_1D), intent(in) :: eop
    class(SpectralElementMesh_3D), intent(in) :: sem
    real(RDP),           intent(in)    :: lambda       !< Helmholtz parameter λ
    real(RDP),           intent(in)    :: nu(:,:,:,:)  !< diffusivity ν
    real(RDP),           intent(in)    :: u(:,:,:,:)   !< 3D scalar field u
    real(RDP),           intent(out)   :: v(:,:,:,:)   !< result
    real(RDP), optional, intent(inout) :: nub(:,:,:,:) !< ν at element boundary
    real(RDP), optional, intent(inout) :: ub(:,:,:,:)  !< element boundary values
    real(RDP), optional, intent(inout) :: qb(:,:,:,:)  !< element boundary fluxes

    contiguous :: nu, u, v, nub, ub, qb

    if (sem % mesh % regular) then
      call TPO_Elliptic_RLVI( Ms     = eop % w              &
                            , Ds     = eop % D              &
                            , dx     = sem % mesh % dx      &
                            , lambda = lambda               &
                            , nu     = nu                   &
                            , u      = u                    &
                            , v      = v                    &
                            , nub    = nub                  &
                            , ub     = ub                   &
                            , qb     = qb                   )
    else
      call TPO_Elliptic_DLVI( Ms     = eop % w              &
                            , Ds     = eop % D              &
                            , Jd     = sem % metrics % Jd   &
                            , G      = sem % metrics % G    &
                            , lambda = lambda               &
                            , nu     = nu                   &
                            , u      = u                    &
                            , v      = v                    &
                            , Ji_n   = sem % metrics % Ji_n &
                            , nub    = nub                  &
                            , ub     = ub                   &
                            , qb     = qb                   )
    end if

  end subroutine TPO_Elliptic_VI_RDP

  !=============================================================================

end module TPO__Elliptic__3D
