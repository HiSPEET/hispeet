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

!> summary:  3D divergence element operator using LIBXSMM for real(RDP) (D)
!> author:   Jerome Michel, Joerg Stiller
!> date:     2021/08/17
!===============================================================================

module TPO__Div__3D_D__XSMM_RDP

#ifdef __LIBXSMM__

  use Kind_Parameters, only: RDP
  use LIBXSMM
  implicit none
  private

  public :: TPO_Div_D_XSMM

  interface TPO_Div_D_XSMM
    module procedure TPO_Div_D_XSMM_RDP
  end interface

contains

  !-----------------------------------------------------------------------------
  !> 3D divergence element operator based on LIBXSMM (R)

  subroutine TPO_Div_D_XSMM_RDP(np, ne, Ms, Ds, Jd, Ji, a, n, u, up, v)
    integer,   intent(in)  :: np                  !< num points per direction
    integer,   intent(in)  :: ne                  !< num elements
    real(RDP), intent(in)  :: Ms(np)              !< 1D standard mass matrix
    real(RDP), intent(in)  :: Ds(np,np)           !< 1D standard diff matrix
    real(RDP), intent(in)  :: Jd(np,np,np,ne)     !< element Jacobian determinant
    real(RDP), intent(in)  :: Ji(np,np,np,ne,3,3) !< inverse Jacobi matrix
    real(RDP), intent(in)  :: a(np,np,6,ne)       !< face area coefficients
    real(RDP), intent(in)  :: n(np,np,6,ne,3)     !< face unit normal vectors
    real(RDP), intent(in)  :: u(np,np,np,ne,3)    !< 3D vector field
    real(RDP), intent(in)  :: up(np,np,6,ne,3)    !< exterior traces u⁺ at faces
    real(RDP), intent(out) :: v(np,np,np,ne)      !< divergence of u

    optional :: Ms, Jd, a, n, up

    !---------------------------------------------------------------------------
    ! Local variables

    type(LIBXSMM_DMMFunction), save :: xmm_1, xmm_2, xmm_3

    real(RDP), parameter :: alpha  = 1
    real(RDP), parameter :: beta   = 0

    real(RDP) :: Ds_t(np,np)

    real(RDP), dimension(np,np,np) :: r, s, t, w

    real(RDP) :: c
    integer   :: e, f, i, j, k
    logical   :: fluxes

    !-----------------------------------------------------------------------------
    ! Initialization

    fluxes = present(Ms) .and. &
             present(Jd) .and. &
             present(a)  .and. &
             present(n)  .and. &
             present(up)

    Ds_t = transpose(Ds)

    ! XSMM dispatch
    !$omp master
    call LIBXSMM_Dispatch(xmm_1, np   , np**2, np, alpha=alpha, beta=beta)
    call LIBXSMM_Dispatch(xmm_2, np   , np   , np, alpha=alpha, beta=beta)
    call LIBXSMM_Dispatch(xmm_3, np**2, np   , np, alpha=alpha, beta=beta)
    !$omp end master
    !$omp barrier

    !---------------------------------------------------------------------------
    ! Evaluation

    !$omp do
    do e = 1, ne

      ! direction 1 ............................................................

      ! r = ∂u₁/∂ξ
      call LIBXSMM_DMMCall(xmm_1, Ds, u(:,:,:,e,1), r)

      ! s = ∂u₁/∂η
      do k = 1, np
        call LIBXSMM_DMMCall(xmm_2, u(:,:,k,e,1), Ds_t, s(:,:,k))
      end do

      ! t = ∂u₁/∂ζ
      call LIBXSMM_DMMCall(xmm_3, u(:,:,:,e,1), Ds_t, t)

      w = r * Ji(:,:,:,e,1,1) &
        + s * Ji(:,:,:,e,2,1) &
        + t * Ji(:,:,:,e,3,1)

      ! direction 2 ............................................................

      ! r = ∂u₂/∂ξ
      call LIBXSMM_DMMCall(xmm_1, Ds, u(:,:,:,e,2), r)

      ! s = ∂u₂/∂η
      do k = 1, np
        call LIBXSMM_DMMCall(xmm_2, u(:,:,k,e,2), Ds_t, s(:,:,k))
      end do

      ! t = ∂u₂/∂ζ
      call LIBXSMM_DMMCall(xmm_3, u(:,:,:,e,2), Ds_t, t)

      w = w + r * Ji(:,:,:,e,1,2) &
            + s * Ji(:,:,:,e,2,2) &
            + t * Ji(:,:,:,e,3,2)

      ! direction 3 ............................................................

      ! r = ∂u₃/∂ξ
      call LIBXSMM_DMMCall(xmm_1, Ds, u(:,:,:,e,3), r)

      ! s = ∂u₂/∂η
      do k = 1, np
        call LIBXSMM_DMMCall(xmm_2, u(:,:,k,e,3), Ds_t, s(:,:,k))
      end do

      ! t = ∂u₂/∂ζ
      call LIBXSMM_DMMCall(xmm_3, u(:,:,:,e,3), Ds_t, t)

      w = w + r * Ji(:,:,:,e,1,3) &
            + s * Ji(:,:,:,e,2,3) &
            + t * Ji(:,:,:,e,3,3)

      v(:,:,:,e) = w

      ! fluxes .................................................................

      if (fluxes) then

        ! faces 1 and 2
        do f = 1, 2
          i = 1 + (np - 1) * (f - 1)
          do k = 1, np
          do j = 1, np
            c = a(j,k,f,e) / (2 * Ms(i ) * Jd(i,j,k,e))
            v(i,j,k,e) = v(i,j,k,e)                                         &
                       + c * n(j,k,f,e,1) * (up(j,k,f,e,1) - u(i,j,k,e,1))  &
                       + c * n(j,k,f,e,2) * (up(j,k,f,e,2) - u(i,j,k,e,2))  &
                       + c * n(j,k,f,e,3) * (up(j,k,f,e,3) - u(i,j,k,e,3))
          end do
          end do
        end do

        ! faces 3 and 4
        do f = 3, 4
          j = 1 + (np - 1) * (f - 3)
          do k = 1, np
          do i = 1, np
            c = a(i,k,f,e) / (2 * Ms(j) * Jd(i,j,k,e))
            v(i,j,k,e) = v(i,j,k,e)                                         &
                       + c * n(i,k,f,e,1) * (up(i,k,f,e,1) - u(i,j,k,e,1))  &
                       + c * n(i,k,f,e,2) * (up(i,k,f,e,2) - u(i,j,k,e,2))  &
                       + c * n(i,k,f,e,3) * (up(i,k,f,e,3) - u(i,j,k,e,3))
          end do
          end do
        end do

        ! faces 5 and 6
        do f = 5, 6
          k = 1 + (np - 1) * (f - 5)
          do j = 1, np
          do i = 1, np
            c = a(i,j,f,e) / (2 * Ms(k) * Jd(i,j,k,e))
            v(i,j,k,e) = v(i,j,k,e)                                         &
                       + c * n(i,j,f,e,1) * (up(i,j,f,e,1) - u(i,j,k,e,1))  &
                       + c * n(i,j,f,e,2) * (up(i,j,f,e,2) - u(i,j,k,e,2))  &
                       + c * n(i,j,f,e,3) * (up(i,j,f,e,3) - u(i,j,k,e,3))
          end do
          end do
        end do

      end if

    end do

  end subroutine TPO_Div_D_XSMM_RDP

  !=============================================================================

#endif

end module TPO__Div__3D_D__XSMM_RDP
