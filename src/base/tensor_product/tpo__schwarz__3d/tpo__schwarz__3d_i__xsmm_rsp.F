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
!> author:   Joerg Stiller
!> date:     2021/06/10
!===============================================================================

module TPO__Schwarz__3D_I__XSMM_RSP

#ifdef __LIBXSMM__

  use LIBXSMM
  use Kind_Parameters, only: RSP
  implicit none
  private

  public :: TPO_Schwarz_I_XSMM

  interface TPO_Schwarz_I_XSMM
    module procedure TPO_Schwarz_I_XSMM_RSP
  end interface

contains

  !-----------------------------------------------------------------------------

  subroutine TPO_Schwarz_I_XSMM_RSP(S, V, W, g, cfg, lambda, nu, f, u)
    real(RSP), intent(in)  :: S(:,:,:)       !< 1D eigenvectors
    real(RSP), intent(in)  :: V(:,:)         !< 1D eigenvalues
    real(RSP), intent(in)  :: W(:,:)         !< 1D weights
    real(RSP), intent(in)  :: g(:,:)         !< metrics
    integer,   intent(in)  :: cfg(:,:)       !< subdomain configs
    real(RSP), intent(in)  :: lambda         !< λ
    real(RSP), intent(in)  :: nu(:)          !< ν + νˢ
    real(RSP), intent(in)  :: f(:,:,:,:)     !< RHS
    real(RSP), intent(out) :: u(:,:,:,:)     !< solution

    real(RSP), parameter :: alpha = 1
    real(RSP), parameter :: beta  = 0
    real(RSP), parameter :: eps   = epsilon(lambda)

    !-----------------------------------------------------------------------------
    ! local variables

    type(LIBXSMM_SMMFunction), save :: xmm_1, xmm_2, xmm_3

    real(RSP) :: S_t  ( size(S,2), size(S,1), size(S,3) )
    real(RSP) :: WS   ( size(S,1), size(S,2), size(S,3) )
    real(RSP) :: WS_t ( size(S,1), size(S,2), size(S,3) )
    real(RSP) :: y    ( size(f,1), size(f,2), size(f,3) )
    real(RSP) :: z    ( size(f,1), size(f,2), size(f,3) )
    real(RSP) :: a1, a2, a3, a4, d

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

      ! coefficients
      a1 = g(1,l) * nu(l)
      a2 = g(2,l) * nu(l)
      a3 = g(3,l) * nu(l)
      a4 = g(4,l) * lambda

      ! y = Sᵀ x I x I f .......................................................
      call LIBXSMM_SMMCall(xmm_3, f(:,:,:,l), S(:,:,c3) , y(:,:,:))

      ! z = I x Sᵀ x I y .......................................................
      do k = 1, np
        call LIBXSMM_SMMCall(xmm_2, y(:,:,k), S(:,:,c2), z(:,:,k))
      end do

      ! y =I x I x Sᵀ z ........................................................
      call LIBXSMM_SMMCall(xmm_1, S_t(:,:,c1), z(:,:,:) , y(:,:,:))

      ! y = D⁻¹ y
      do k = 1, np
      do j = 1, np
      do i = 1, np
        d = a1 * V(i,c1) + a2 * V(j,c2) + a3 * V(k,c3) + a4
        y(i,j,k) = y(i,j,k) / sign(max(abs(d),eps), d)
      end do
      end do
      end do

      ! z = I x I x WS y .......................................................
      call LIBXSMM_SMMCall(xmm_1, WS(:,:,c1), y(:,:,:) , z(:,:,:))

      ! y = I x WS x I z .......................................................
      do k = 1, np
        call LIBXSMM_SMMCall(xmm_2, z(:,:,k), WS_t(:,:,c2), y(:,:,k))
      end do

      ! u = WS x I x I y .......................................................
      call LIBXSMM_SMMCall(xmm_3, y(:,:,:), WS_t(:,:,c3) , u(:,:,:,l))

    end do subdomains

  end subroutine TPO_Schwarz_I_XSMM_RSP

  !=============================================================================

#endif

end module TPO__Schwarz__3D_I__XSMM_RSP
