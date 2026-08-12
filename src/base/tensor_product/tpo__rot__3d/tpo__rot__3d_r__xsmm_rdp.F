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

!> summary:  3D rotation element operator using LIBXSSM (R)
!> author:   Erik Pfister
!> date:     2019/11/06
!===============================================================================

module TPO__Rot__3D_R__XSMM_RDP

#ifdef __LIBXSMM__

  use Kind_Parameters, only: RDP
  use LIBXSMM
  implicit none
  private

  public :: TPO_Rot_R_XSMM

  interface TPO_Rot_R_XSMM
    module procedure TPO_Rot_R_XSMM_RDP
  end interface

contains

  !-----------------------------------------------------------------------------
  !> 3D rotation element operator based on LIBXSMM (R)

  subroutine TPO_Rot_R_XSMM_RDP(np, ne, Ds, dx, u, v)
    integer,   intent(in)  :: np               !< num points per direction
    integer,   intent(in)  :: ne               !< num elements
    real(RDP), intent(in)  :: Ds(np, np)       !< standard diff matrix
    real(RDP), intent(in)  :: dx(3)            !< element extensions
    real(RDP), intent(in)  :: u(np,np,np,ne,3) !< 3D vector field
    real(RDP), intent(out) :: v(np,np,np,ne,3) !< element-wise rotation of u

    type(LIBXSMM_DMMFunction), save :: xmm_1, xmm_2, xmm_3

    real(RDP) :: Ds_1(np,np)
    real(RDP) :: Ds_2t(np,np), Ds_3t(np, np)
    real(RDP), parameter :: alpha  = 1
    real(RDP), parameter :: beta   = 1

    real(RDP) :: g(3)
    integer   :: e, k

    !-----------------------------------------------------------------------------

    ! metric coefficients
    g = 2 / dx

    Ds_1  = g(1) * Ds
    Ds_2t = g(2) * transpose(Ds)
    Ds_3t = g(3) * transpose(Ds)

    ! XSMM dispatch
    !$omp master
    call LIBXSMM_Dispatch(xmm_1, np   , np**2, np, alpha=alpha, beta=beta)
    call LIBXSMM_Dispatch(xmm_2, np   , np   , np, alpha=alpha, beta=beta)
    call LIBXSMM_Dispatch(xmm_3, np**2, np   , np, alpha=alpha, beta=beta)

    if (.not. ( libxsmm_available(xmm_1) .and.  &
                libxsmm_available(xmm_2) .and.  &
                libxsmm_available(xmm_3)        ) ) then

      stop "TPO_Rot_R_XSMM: LIBXSMM_Dispatch failed"

    end if
    !$omp end master
    !$omp barrier

    !---------------------------------------------------------------------------
    ! evaluation

    !$omp do
    do e = 1, ne

      v(:,:,:,e,1) = 0
      v(:,:,:,e,2) = 0
      v(:,:,:,e,3) = 0

      ! v2 = du1/dx3
      call LIBXSMM_DMMCall(xmm_3, u(:,:,:,e,1), Ds_3t,  v(:,:,:,e,2))

      ! v3 =-du1/dx2
      do k = 1, np
        call LIBXSMM_DMMCall(xmm_2, u(:,:,k,e,1), -Ds_2t, v(:,:,k,e,3))
      end do

      ! v3+= du2/dx1
      call LIBXSMM_DMMCall(xmm_1, Ds_1, u(:,:,:,e,2), v(:,:,:,e,3))

      ! v1 =-du1/dx3
      call LIBXSMM_DMMCall(xmm_3, u(:,:,:,e,2), -Ds_3t,  v(:,:,:,e,1))

      ! v1+= du1/dx2
      do k = 1, np
        call LIBXSMM_DMMCall(xmm_2, u(:,:,k,e,3), Ds_2t, v(:,:,k,e,1))
      end do

      ! v2+=-du3/dx1
      call LIBXSMM_DMMCall(xmm_1, -Ds_1, u(:,:,:,e,3), v(:,:,:,e,2))

    end do

  end subroutine TPO_Rot_R_XSMM_RDP

  !=============================================================================

#endif

end module TPO__Rot__3D_R__XSMM_RDP
