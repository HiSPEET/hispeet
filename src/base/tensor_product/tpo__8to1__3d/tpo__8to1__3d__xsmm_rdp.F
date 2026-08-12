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

!> summary:   3d regular 8:1 projection operator using LIBXSMM for real(RDP)
!> author:    Joerg Stiller
!> date:      2024/08/07
!===============================================================================

module TPO__8To1__3D__XSMM_RDP

#ifdef __LIBXSMM__

  use Kind_Parameters, only: RDP
  use LIBXSMM
  implicit none
  private

  public :: TPO_8To1_XSMM

  interface TPO_8To1_XSMM
    module procedure TPO_8To1_XSMM_RDP
  end interface

contains

  !-----------------------------------------------------------------------------
  !> XSMM AxAxA operator

  subroutine TPO_8To1_XSMM_RDP(A, u, v)
    real(RDP), contiguous, intent(in)  :: A(:,:,:)   !< 1D projection operator
    real(RDP), contiguous, intent(in)  :: u(:,:,:,:) !< operand
    real(RDP), contiguous, intent(out) :: v(:,:,:,:) !< result

    real(RDP), parameter :: ZERO = 0
    real(RDP), parameter :: ONE  = 1

    !---------------------------------------------------------------------------
    ! local variables

    type(LIBXSMM_DMMFunction), save :: xmm_set_1, xmm_add_1
    type(LIBXSMM_DMMFunction), save :: xmm_set_2, xmm_add_2
    type(LIBXSMM_DMMFunction), save :: xmm_set_3, xmm_add_3

    real(RDP) :: At ( size(A,2), size(A,1), 2 )
    real(RDP) :: w  ( size(u,1), size(u,1), size(u,1), 2, 2, 2 )
    real(RDP) :: z3 ( size(u,1), size(u,1), size(v,1), 2, 2 )
    real(RDP) :: z2 ( size(u,1), size(v,1), size(v,1), 2 )

    integer   :: np, nq
    integer   :: e, k, l, m, o

    !---------------------------------------------------------------------------
    ! initialization

    np = size(u,1) ! = size(A,2)
    nq = size(v,1) ! = size(A,1)

    At(:,:,1) = transpose(A(:,:,1))
    At(:,:,2) = transpose(A(:,:,2))

    !-----------------------------------------------------------------------------
    ! XSMM dispatch
    !$omp master
    call LIBXSMM_Dispatch(xmm_set_3, np**2, nq   , np, alpha=ONE, beta=ZERO)
    call LIBXSMM_Dispatch(xmm_add_3, np**2, nq   , np, alpha=ONE, beta=ONE )
    call LIBXSMM_Dispatch(xmm_set_2, np   , nq   , np, alpha=ONE, beta=ZERO)
    call LIBXSMM_Dispatch(xmm_add_2, np   , nq   , np, alpha=ONE, beta=ONE )
    call LIBXSMM_Dispatch(xmm_set_1, nq   , nq**2, np, alpha=ONE, beta=ZERO)
    call LIBXSMM_Dispatch(xmm_add_1, nq   , nq**2, np, alpha=ONE, beta=ONE )
    !$omp end master
    !$omp barrier

    !---------------------------------------------------------------------------
    ! evaluation

    !$omp do
    do e = 1, size(v,4)

      ! preliminaries ..........................................................

      ! offset of fine element index
      o = 8 * (e-1)

      ! extract 2x2x2 cluster of fine element variable
      w(:,:,:,1,1,1) = u(:,:,:,o+1)
      w(:,:,:,2,1,1) = u(:,:,:,o+2)
      w(:,:,:,1,2,1) = u(:,:,:,o+3)
      w(:,:,:,2,2,1) = u(:,:,:,o+4)
      w(:,:,:,1,1,2) = u(:,:,:,o+5)
      w(:,:,:,2,1,2) = u(:,:,:,o+6)
      w(:,:,:,1,2,2) = u(:,:,:,o+7)
      w(:,:,:,2,2,2) = u(:,:,:,o+8)

      ! z3 = AxIxI w ...........................................................

      do m = 1, 2
      do l = 1, 2
        call LIBXSMM_DMMCall(xmm_set_3, w(:,:,:,l,m,1), At(:,:,1), z3(:,:,:,l,m))
        call LIBXSMM_DMMCall(xmm_add_3, w(:,:,:,l,m,2), At(:,:,2), z3(:,:,:,l,m))
      end do
      end do

      ! z2 = IxAxI z3 ..........................................................

      do l = 1, 2
      do k = 1, nq
        call LIBXSMM_DMMCall(xmm_set_2, z3(:,:,k,l,1), At(:,:,1), z2(:,:,k,l))
        call LIBXSMM_DMMCall(xmm_add_2, z3(:,:,k,l,2), At(:,:,2), z2(:,:,k,l))
      end do
      end do

      ! v^e = IxIxA z2 .........................................................

      call LIBXSMM_DMMCall(xmm_set_1, A(:,:,1), z2(:,:,:,1), v(:,:,:,e))
      call LIBXSMM_DMMCall(xmm_add_1, A(:,:,2), z2(:,:,:,2), v(:,:,:,e))

    end do

  end subroutine TPO_8To1_XSMM_RDP

  !=============================================================================

#endif

end module TPO__8To1__3D__XSMM_RDP
