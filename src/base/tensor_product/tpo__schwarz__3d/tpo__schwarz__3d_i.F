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

!> summary:  Evaluation of 3D constant isotropic Schwarz TPO
!> author:   Joerg Stiller
!> date:     2020/06/02
!===============================================================================

module TPO__Schwarz__3D_I
  use Kind_Parameters, only: RSP, RDP
  use TPO__Schwarz__3D_I__Hand
  use TPO__Schwarz__3D_I__XSMM
  implicit none
  private

  public :: TPO_Schwarz

  interface TPO_Schwarz
    module procedure TPO_Schwarz_I_RSP
    module procedure TPO_Schwarz_I_RDP
  end interface

contains

  !-----------------------------------------------------------------------------
  !> 3D constant isotropic operator for real(RSP)

  subroutine TPO_Schwarz_I_RSP(S, V, W, g, cfg, lambda, nu, f, u)
    real(RSP), intent(in)  :: S(:,:,:)       !< 1D eigenvectors   (np,np,nc)
    real(RSP), intent(in)  :: V(:,:)         !< 1D eigenvalues    (np,nc)
    real(RSP), intent(in)  :: W(:,:)         !< 1D weights        (np,nc)
    real(RSP), intent(in)  :: g(:,:)         !< metrics           (4,nd)
    integer,   intent(in)  :: cfg(:,:)       !< subdomain configs (3,nd)
    real(RSP), intent(in)  :: lambda         !< λ
    real(RSP), intent(in)  :: nu(:)          !< ν + νˢ
    real(RSP), intent(in)  :: f(:,:,:,:)     !< RHS               (np,np,np,nd)
    real(RSP), intent(out) :: u(:,:,:,:)     !< solution          (np,np,np,nd)

    select case(size(S,1))

#ifdef __LIBXSMM__
    case(8:)
      call TPO_Schwarz_I_XSMM(S, V, W, g, cfg, lambda, nu, f, u)
#endif

    case default
      call TPO_Schwarz_I_Hand(S, V, W, g, cfg, lambda, nu, f, u)

    end select

  end subroutine TPO_Schwarz_I_RSP

  !-----------------------------------------------------------------------------
  !> 3D constant isotropic operator for real(RDP)

  subroutine TPO_Schwarz_I_RDP(S, V, W, g, cfg, lambda, nu, f, u)
    real(RDP), intent(in)  :: S(:,:,:)       !< 1D eigenvectors
    real(RDP), intent(in)  :: V(:,:)         !< 1D eigenvalues
    real(RDP), intent(in)  :: W(:,:)         !< 1D weights
    real(RDP), intent(in)  :: g(:,:)         !< metrics
    integer,   intent(in)  :: cfg(:,:)       !< subdomain configs
    real(RDP), intent(in)  :: lambda         !< λ
    real(RDP), intent(in)  :: nu(:)          !< ν + νˢ
    real(RDP), intent(in)  :: f(:,:,:,:)     !< RHS
    real(RDP), intent(out) :: u(:,:,:,:)     !< solution

    select case(size(S,1))

#ifdef __LIBXSMM__
    case(8:)
      call TPO_Schwarz_I_XSMM(S, V, W, g, cfg, lambda, nu, f, u)
#endif

    case default
      call TPO_Schwarz_I_Hand(S, V, W, g, cfg, lambda, nu, f, u)

    end select

  end subroutine TPO_Schwarz_I_RDP

  !=============================================================================

end module TPO__Schwarz__3D_I
