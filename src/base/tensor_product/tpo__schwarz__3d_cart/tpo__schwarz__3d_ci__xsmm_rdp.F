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

!> summary:  3d generic isotropic Schwarz operator using LIBXSMM
!> author:   Joerg Stiller, Erik Pfister
!> date:     2020/06/02
!===============================================================================

module TPO__Schwarz__3D_CI__XSMM_RDP

#ifdef __LIBXSMM__

  use LIBXSMM
  use Kind_Parameters, only: RDP
  implicit none
  private

  public :: TPO_Schwarz_CI_XSMM

  interface TPO_Schwarz_CI_XSMM
    module procedure TPO_Schwarz_CI_XSMM_RDP
  end interface

contains

  !-----------------------------------------------------------------------------

  subroutine TPO_Schwarz_CI_XSMM_RDP(S, W, cfg, D_inv, f, u)
    real(RDP), intent(in)  :: S(:,:,:)       !< 1D eigenvectors   (np,np,nc)
    real(RDP), intent(in)  :: W(:,:)         !< 1D weights        (np,nc)
    integer,   intent(in)  :: cfg(:,:)       !< subdomain configs (3,nd)
    real(RDP), intent(in)  :: D_inv(:,:,:,:) !< inv. 3D eigenvals (np,np,np,nd)
    real(RDP), intent(in)  :: f(:,:,:,:)     !< RHS               (np,np,np,nd)
    real(RDP), intent(out) :: u(:,:,:,:)     !< solution          (np,np,np,nd)

    real(RDP), parameter :: alpha  = 1
    real(RDP), parameter :: beta   = 0

    !-----------------------------------------------------------------------------
    ! local variables

    type(LIBXSMM_DMMFunction), save :: xmm_1, xmm_2, xmm_3

    real(RDP) :: S_t  ( size(S,2), size(S,1), size(S,3) )
    real(RDP) :: WS   ( size(S,1), size(S,2), size(S,3) )
    real(RDP) :: WS_t ( size(S,1), size(S,2), size(S,3) )
    real(RDP) :: y    ( size(f,1), size(f,2), size(f,3) )
    real(RDP) :: z    ( size(f,1), size(f,2), size(f,3) )

    integer :: nc, nd, np
    integer :: i, j, k, l, p
    integer :: c, c1, c2, c3

    !-----------------------------------------------------------------------------
    ! initialization

    np = size(S,1)
    nc = size(S,3)
    nd = size(f,4)

    ! WS(:,:,c) = (WS)(:,:,c)
    do c = 1, nc
    do j = 1, np
    do i = 1, np
      WS(i,j,c) = S(i,j,c) * W(i,c)
    end do
    end do
    end do

    ! WS_t(:,:,c) = (WS)ᵀ(:,:,c)
    do c = 1, nc
    do j = 1, np
    do i = 1, np
      WS_t(j,i,c) = W(i,c) * S(i,j,c)
    end do
    end do
    end do

    do c = 1, nc
    do j = 1, np
    do i = 1, np
      S_t(j,i,c) = S(i,j,c)
    end do
    end do
    end do

    ! XSMM dispatch
    !$omp master
    call LIBXSMM_Dispatch(xmm_1, np   , np**2, np, alpha=alpha, beta=beta)
    call LIBXSMM_Dispatch(xmm_2, np   , np   , np, alpha=alpha, beta=beta)
    call LIBXSMM_Dispatch(xmm_3, np**2, np   , np, alpha=alpha, beta=beta)

    if (.not. ( libxsmm_available(xmm_1) .and.  &
                libxsmm_available(xmm_2) .and.  &
                libxsmm_available(xmm_3)        ) ) then

      stop "TPO_Div_R_XSMM: LIBXSMM_Dispatch failed"

    end if
    !$omp end master
    !$omp barrier

    !---------------------------------------------------------------------------
    ! evaluation

    !$omp do private(l)
    subdomains: do l = 1, nd

      ! element-boundary configurations
      c1 = cfg(1,l)
      c2 = cfg(2,l)
      c3 = cfg(3,l)

      ! y = Sᵀ x I x I f .......................................................
      call LIBXSMM_DMMCall(xmm_3, f(:,:,:,l), S(:,:,c3) , y(:,:,:))

      ! z = I x Sᵀ x I y .......................................................
      do k = 1, np
        call LIBXSMM_DMMCall(xmm_2, y(:,:,k), S(:,:,c2), z(:,:,k))
      end do

      ! y =I x I x Sᵀ z ........................................................
      call LIBXSMM_DMMCall(xmm_1, S_t(:,:,c1), z(:,:,:) , y(:,:,:))

      ! y = D⁻¹ y
      y = D_inv(:,:,:,l) * y

      ! z = I x I x WS y .......................................................
      call LIBXSMM_DMMCall(xmm_1, WS(:,:,c1), y(:,:,:) , z(:,:,:))

      ! y = I x WS x I z .......................................................
      do k = 1, np
        call LIBXSMM_DMMCall(xmm_2, z(:,:,k), WS_t(:,:,c2), y(:,:,k))
      end do

      ! u = WS x I x I y .......................................................
      call LIBXSMM_DMMCall(xmm_3, y(:,:,:), WS_t(:,:,c3) , u(:,:,:,l))

    end do subdomains

  end subroutine TPO_Schwarz_CI_XSMM_RDP

  !=============================================================================

#endif

end module TPO__Schwarz__3D_CI__XSMM_RDP
