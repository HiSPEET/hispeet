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

!> summary:  3d generic element elliptic operator (RLVI) using LIBXSMM
!> author:   Erik Pfister, Joerg Stiller
!> date:     2020/06/08
!===============================================================================

module TPO__Elliptic__3D_RLVI__XSMM_RDP

#ifdef __LIBXSMM__

  use Kind_Parameters, only: RDP
  use LIBXSMM
  implicit none
  private

  public :: TPO_Elliptic_RLVI_XSMM

  interface TPO_Elliptic_RLVI_XSMM
    module procedure TPO_Elliptic_RLVI_XSMM_RDP
  end interface

contains

  !-----------------------------------------------------------------------------
  !> Generic element elliptic operator (RLVI)

  subroutine TPO_Elliptic_RLVI_XSMM_RDP( Ms, Ds, dx, lambda &
                                       , nu,  u,  v         &
                                       , nub, ub, qb        )

    real(RDP), intent(in)  :: Ms(:)       !< 1D standard mass matrix        (np)
    real(RDP), intent(in)  :: Ds(:,:)     !< 1D standard diff matrix     (np,np)
    real(RDP), intent(in)  :: lambda      !< Helmholtz parameter
    real(RDP), intent(in)  :: dx(3)       !< element extensions
    real(RDP), intent(in)  :: nu(:,:,:,:) !< diffusivity           (np,np,np,ne)
    real(RDP), intent(in)  :: u(:,:,:,:)  !< operand               (np,np,np,ne)
    real(RDP), intent(out) :: v(:,:,:,:)  !< result                (np,np,np,ne)

    real(RDP), optional, intent(inout) :: nub(:,:,:,:) !< ν at element boundary
    real(RDP), optional, intent(inout) :: ub (:,:,:,:) !< u at element boundary
    real(RDP), optional, intent(inout) :: qb (:,:,:,:) !< element boundary flux

    contiguous :: Ms, Ds, nu, u, v, nub, ub, qb

    !---------------------------------------------------------------------------
    ! local variables

    type(LIBXSMM_DMMFunction), save :: xmm_set_1, xmm_set_2, xmm_set_3
    type(LIBXSMM_DMMFunction), save :: xmm_add_1, xmm_add_2, xmm_add_3

    real(RDP), parameter :: ZERO = 0
    real(RDP), parameter :: ONE  = 1

    real(RDP), dimension(size(Ms), size(Ms), size(Ms)) :: M, M_nu, z
    real(RDP), dimension(size(Ms), size(Ms)) :: Ds_t

    real(RDP) :: g(3), tmp
    integer   :: np, ne
    integer   :: e, i, j, k
    logical   :: get_traces

    !---------------------------------------------------------------------------
    ! initialization

    np = size(Ms)
    ne = size(u,4)

    get_traces = present(nub) .and. present(ub) .and. present(qb)

    ! element mass matrix
    tmp = product(dx) / 8
    do k = 1, np
    do j = 1, np
    do i = 1, np
      M(i,j,k) = tmp * Ms(k) * Ms(j) * Ms(i)
    end do
    end do
    end do

    ! metric coefficients
    g = 4 / dx**2

    ! initialize operators
    Ds_t = transpose(Ds)

    ! XSMM dispatch
    !$omp master
    call LIBXSMM_Dispatch(xmm_set_1, np   , np**2, np, alpha=ONE, beta=ZERO)
    call LIBXSMM_Dispatch(xmm_add_1, np   , np**2, np, alpha=ONE, beta=ONE )
    call LIBXSMM_Dispatch(xmm_set_2, np   , np   , np, alpha=ONE, beta=ZERO)
    call LIBXSMM_Dispatch(xmm_add_2, np   , np   , np, alpha=ONE, beta=ONE )
    call LIBXSMM_Dispatch(xmm_set_3, np**2, np   , np, alpha=ONE, beta=ZERO)
    call LIBXSMM_Dispatch(xmm_add_3, np**2, np   , np, alpha=ONE, beta=ONE )
    !$omp end master
    !$omp barrier

    ! coefficients
    g = 4 / dx**2

    !---------------------------------------------------------------------------
    ! evaluation

    !$omp do private(e)
    do e = 1, ne

      ! vᵉ = λ Mᵉ uᵉ, M_nu = Mᵉνᵉ ..............................................

      do k = 1, np
      do j = 1, np
      do i = 1, np
        v(i,j,k,e) = lambda * M(i,j,k) * u(i,j,k,e)
        M_nu(i,j,k) = M(i,j,k) * nu(i,j,k,e)
      end do
      end do
      end do

      ! direction 1 ............................................................

      ! z = [I x I x Dˢ] uᵉ
      call LIBXSMM_DMMCall(xmm_set_1, Ds, u(:,:,:,e), z)

      if (get_traces) then
        ! nub = ν;  ub = u;  qb = n⋅ν∇u   @ Γ₁ + Γ₂
        do k = 1, np
        do j = 1, np
          nub(j,k,1,e) =  nu( 1,j,k,e)
          nub(j,k,2,e) =  nu(np,j,k,e)
          ub (j,k,1,e) =  u ( 1,j,k,e)
          ub (j,k,2,e) =  u (np,j,k,e)
          qb (j,k,1,e) = -2/dx(1) * nu( 1,j,k,e) * z( 1,j,k)
          qb (j,k,2,e) =  2/dx(1) * nu(np,j,k,e) * z(np,j,k)
        end do
        end do
      end if

      ! z = g₁ [νᵉMᵉ] z
      do k = 1, np
      do j = 1, np
      do i = 1, np
        z(i,j,k) = g(1) * M_nu(i,j,k) * z(i,j,k)
      end do
      end do
      end do

      ! vᵉ += [I x I x (Dˢ)ᵀ] z
      call LIBXSMM_DMMCall(xmm_add_1, Ds_t, z, v(:,:,:,e))

      ! direction 2 ............................................................

      do k = 1, np

        ! z = [I x Dˢ x I] uᵉ
        call LIBXSMM_DMMCall(xmm_set_2, u(:,:,k,e), Ds_t, z(:,:,k))

        if (get_traces) then
         ! nub = ν;  ub = u;  qb = n⋅ν∇u  @ Γ₃ + Γ₄
         do i = 1, np
            nub(i,k,3,e) =  nu(i, 1,k,e)
            nub(i,k,4,e) =  nu(i,np,k,e)
            ub (i,k,3,e) =  u (i, 1,k,e)
            ub (i,k,4,e) =  u (i,np,k,e)
            qb (i,k,3,e) = -2/dx(2) * nu(i, 1,k,e) * z(i, 1,k)
            qb (i,k,4,e) =  2/dx(2) * nu(i,np,k,e) * z(i,np,k)
          end do
        end if

        ! z = g₂ [νᵉMᵉ] z
        do j = 1, np
        do i = 1, np
          z(i,j,k) = g(2) * M_nu(i,j,k) * z(i,j,k)
        end do
        end do

        ! vᵉ += [I x (Dˢ)ᵀ x I] z
        call LIBXSMM_DMMCall(xmm_add_2, z(:,:,k), Ds, v(:,:,k,e))
      end do

      ! direction 3 ............................................................

      ! z = [Dˢ x I x I] uᵉ
      call LIBXSMM_DMMCall(xmm_set_3, u(:,:,:,e), Ds_t, z)

      if (get_traces) then
        ! nub = ν;  ub = u;  qb = n⋅ν∇u   @ Γ₅ + Γ₆
        do j = 1, np
        do i = 1, np
          nub(i,j,5,e) =  nu(i,j, 1,e)
          nub(i,j,6,e) =  nu(i,j,np,e)
          ub (i,j,5,e) =  u (i,j, 1,e)
          ub (i,j,6,e) =  u (i,j,np,e)
          qb (i,j,5,e) = -2/dx(3) * nu(i,j, 1,e) * z(i,j, 1)
          qb (i,j,6,e) =  2/dx(3) * nu(i,j,np,e) * z(i,j,np)
        end do
        end do
      end if

      ! z = g₃ [νᵉMᵉ] z
      do k = 1, np
      do j = 1, np
      do i = 1, np
        z(i,j,k) = g(3) * M_nu(i,j,k) * z(i,j,k)
      end do
      end do
      end do

      ! vᵉ += [(Dˢ)ᵀ x I x I] z
      call LIBXSMM_DMMCall(xmm_add_3, z, Ds, v(:,:,:,e))

    end do
    !$omp end do

  end subroutine TPO_Elliptic_RLVI_XSMM_RDP

  !=============================================================================

#endif

end module TPO__Elliptic__3D_RLVI__XSMM_RDP
