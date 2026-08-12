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

!> summary:  Evaluation of  v = [(S₃xS₂xS₁) Λ (S₃xS₂xS₁)ᵀ] u
!>           with 3D mesh variables
!> author:   Joerg Stiller
!> date:     2020/06/02
!===============================================================================

module TPO__Spectral__3D_CA
  use Kind_Parameters, only: RDP
  use TPO__Spectral__3D_CA__Gen
  use TPO__Spectral__3D_CA__XSMM
  implicit none
  private

  public :: TPO_Spectral

  interface TPO_Spectral
    module procedure TPO_Spectral_CA_RDP
  end interface

contains

  !-----------------------------------------------------------------------------
  !> 3D (S₃ x S₂ x S₁) Λ (S₃ x S₂ x S₁)ᵀ operator

  subroutine TPO_Spectral_CA_RDP(S1, S2, S3, Lambda, u, v)
    real(RDP), intent(in)  :: S1(:,:)       !< eigenvectors for dir 1: S₁(n1,n1)
    real(RDP), intent(in)  :: S2(:,:)       !< eigenvectors for dir 2: S₂(n2,n2)
    real(RDP), intent(in)  :: S3(:,:)       !< eigenvectors for dir 3: S₃(n3,n3)
    real(RDP), intent(in)  :: Lambda(:,:,:) !< 3D eigenvalues: Λ(n1,n2,n3)
    real(RDP), intent(in)  :: u(:,:,:,:)    !< operand: u(n1.n2,n3,ne)
    real(RDP), intent(out) :: v(:,:,:,:)    !< result:  v(n1,n2,n3,ne)

    select case(min( size(S1,1), size(S2,1), size(S3,1) ))

#ifdef __LIBXSMM__
    case(6:)
      call TPO_Spectral_CA_XSMM(S1, S2, S3, Lambda, u, v)
#endif

    case default
      call TPO_Spectral_CA_Gen(S1, S2, S3, Lambda, u, v)

    end select

  end subroutine TPO_Spectral_CA_RDP

  !=============================================================================

end module TPO__Spectral__3D_CA
