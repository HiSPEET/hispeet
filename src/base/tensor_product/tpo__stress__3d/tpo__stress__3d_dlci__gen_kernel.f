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

!> summary:  3d generic curvilinear element stress operator
!> author:   Joerg Stiller
!> date:     2021/12/23
!===============================================================================

subroutine TPO_Stress_DLCI_Gen_RWP(Ms, Ds, Jd, Ji, n, zeta, eta, v, fv, vb, sb)

  real(RWP), contiguous, intent(in) :: Ms(:)
  !< 1D diagonal standard mass matrix (np)

  real(RWP), contiguous, intent(in) :: Ds(:,:)
  !< 1D standard diff matrix (np,np)

  real(RWP), contiguous, intent(in) :: Jd(:,:,:,:)
  !< Jacobian determinant (np,np,np,ne)

  real(RWP), contiguous, intent(in) :: Ji(:,:,:,:,:,:)
  !< inverse Jacobian matrix  (np,np,np,ne,3,3)

  real(RWP), contiguous, intent(in) :: n(:,:,:,:,:)
  !< unit normal vector at element boundary faces (np,np,6,ne,3)

  real(RWP), intent(in) :: zeta
  !< bulk viscosity or modulus ζ

  real(RWP), intent(in) :: eta
  !< shear viscosity or modulus η

  real(RWP), contiguous, intent(in) :: v(:,:,:,:,:)
  !< velocity or displacement (np,np,np,ne,3)

  real(RWP), contiguous, intent(out) :: fv(:,:,:,:,:)
  !< element volume contribution to stress operator (np,np,np,ne,3)

  real(RWP), contiguous, optional, intent(out) :: vb(:,:,:,:,:)
  !< element boundary values of velocity or displacement (np,np,6,ne+ng,3),

  real(RWP), contiguous, optional, intent(out) :: sb(:,:,:,:,:)
  !< stress normal to element boundary fluxes, sb = n⁻⋅τ⁻ (np,np,6,ne+ng,3)

  !-----------------------------------------------------------------------------
  ! local variables

  real(RWP), parameter :: HALF      = 0.5_RWP
  real(RWP), parameter :: TWO_THIRD = 2.0_RWP / 3.0_RWP

  real(RWP), dimension(size(Ms), size(Ms), size(Ms))       :: M
  real(RWP), dimension(size(Ms), size(Ms), size(Ms), 3)    :: w, z
  real(RWP), dimension(size(Ms), size(Ms), size(Ms), 3, 3) :: tau

  real(RWP) :: chi, dv(3), tmp
  integer   :: c, e, f, i, j, k, p, np, ne
  logical   :: get_traces

  !-----------------------------------------------------------------------------
  ! initialization

  np = size(Ms)
  ne = size(v,4)

  chi = zeta - TWO_THIRD * eta

  get_traces = present(vb) .and. present(sb)

  !-----------------------------------------------------------------------------
  ! evaluation

  !$omp do private(e)
  do e = 1, ne

    ! element mass matrix ......................................................

    do k = 1, np
    do j = 1, np
    do i = 1, np
      M(i,j,k) = Ms(k) * Ms(j) * Ms(i) * Jd(i,j,k,e)
    end do
    end do
    end do

    ! viscous stress tensor ....................................................

    tau = 0

    do c = 1, 3

      ! w₁ = [I x I x Dˢ] v(:,:,:,e,c) = ∂vc/∂ξ₁
      do k = 1, np
      do j = 1, np
      do i = 1, np
        tmp = 0
        do p = 1, np
          tmp = tmp + Ds(i,p) * v(p,j,k,e,c)
        end do
        w(i,j,k,1) = tmp
      end do
      end do
      end do

      ! w₂ = [I x Dˢ x I] v(:,:,:,e,c) = ∂vc/∂ξ₂
      do k = 1, np
      do j = 1, np
      do i = 1, np
        tmp = 0
        do p = 1, np
          tmp = tmp + Ds(j,p) * v(i,p,k,e,c)
        end do
        w(i,j,k,2) = tmp
      end do
      end do
      end do

      ! w₃ = [Dˢ x I x I] v(:,:,:,e,c) = ∂vc/∂ξ₃
      do k = 1, np
      do j = 1, np
      do i = 1, np
        tmp = 0
        do p = 1, np
          tmp = tmp + Ds(k,p) * v(i,j,p,e,c)
        end do
        w(i,j,k,3) = tmp
      end do
      end do
      end do

      ! tau = η(∇v + ∇vᵀ) + χ∇⋅v
      do k = 1, np
      do j = 1, np
      do i = 1, np

        dv(1) = w(i,j,k,1) * Ji(i,j,k,e,1,1) &
              + w(i,j,k,2) * Ji(i,j,k,e,2,1) &
              + w(i,j,k,3) * Ji(i,j,k,e,3,1)

        dv(2) = w(i,j,k,1) * Ji(i,j,k,e,1,2) &
              + w(i,j,k,2) * Ji(i,j,k,e,2,2) &
              + w(i,j,k,3) * Ji(i,j,k,e,3,2)

        dv(3) = w(i,j,k,1) * Ji(i,j,k,e,1,3) &
              + w(i,j,k,2) * Ji(i,j,k,e,2,3) &
              + w(i,j,k,3) * Ji(i,j,k,e,3,3)

        tau(i,j,k,1,c) = tau(i,j,k,1,c) + eta * dv(1)
        tau(i,j,k,2,c) = tau(i,j,k,2,c) + eta * dv(2)
        tau(i,j,k,3,c) = tau(i,j,k,3,c) + eta * dv(3)

        tau(i,j,k,c,1) = tau(i,j,k,c,1) + eta * dv(1)
        tau(i,j,k,c,2) = tau(i,j,k,c,2) + eta * dv(2)
        tau(i,j,k,c,3) = tau(i,j,k,c,3) + eta * dv(3)

        tau(i,j,k,1,1) = tau(i,j,k,1,1) + chi * dv(c)
        tau(i,j,k,2,2) = tau(i,j,k,2,2) + chi * dv(c)
        tau(i,j,k,3,3) = tau(i,j,k,3,3) + chi * dv(c)

      end do
      end do
      end do

    end do

    ! volume integral ..........................................................

    ! w = M ΣᵣJ⁻¹(1,r)⋅τ(r,1:3), exploiting symmetry of τ
    do k = 1, np
    do j = 1, np
    do i = 1, np

      w(i,j,k,1) = M(i,j,k) * ( Ji(i,j,k,e,1,1) * tau(i,j,k,1,1) &
                              + Ji(i,j,k,e,1,2) * tau(i,j,k,2,1) &
                              + Ji(i,j,k,e,1,3) * tau(i,j,k,3,1) )

      w(i,j,k,2) = M(i,j,k) * ( Ji(i,j,k,e,1,1) * tau(i,j,k,2,1) & ! = τ₁₂
                              + Ji(i,j,k,e,1,2) * tau(i,j,k,2,2) & ! τ₂₂
                              + Ji(i,j,k,e,1,3) * tau(i,j,k,3,2) )

      w(i,j,k,3) = M(i,j,k) * ( Ji(i,j,k,e,1,1) * tau(i,j,k,3,1) & ! = τ₁₃
                              + Ji(i,j,k,e,1,2) * tau(i,j,k,3,2) & ! = τ₂₃
                              + Ji(i,j,k,e,1,3) * tau(i,j,k,3,3) )
    end do
    end do
    end do

    ! z = [I x I x (Dˢ)ᵀ] w
    do k = 1, np
    do j = 1, np
    do i = 1, np
      z(i,j,k,1) = 0
      z(i,j,k,2) = 0
      z(i,j,k,3) = 0
      do p = 1, np
        z(i,j,k,1) = z(i,j,k,1) + Ds(p,i) * w(p,j,k,1)
        z(i,j,k,2) = z(i,j,k,2) + Ds(p,i) * w(p,j,k,2)
        z(i,j,k,3) = z(i,j,k,3) + Ds(p,i) * w(p,j,k,3)
      end do
    end do
    end do
    end do

    ! w = M ΣᵣJ⁻¹(2,r)⋅τ(r,1:3), exploiting symmetry of τ
    do k = 1, np
    do j = 1, np
    do i = 1, np

      w(i,j,k,1) = M(i,j,k) * ( Ji(i,j,k,e,2,1) * tau(i,j,k,1,1) &
                              + Ji(i,j,k,e,2,2) * tau(i,j,k,2,1) &
                              + Ji(i,j,k,e,2,3) * tau(i,j,k,3,1) )

      w(i,j,k,2) = M(i,j,k) * ( Ji(i,j,k,e,2,1) * tau(i,j,k,2,1) & ! = τ₁₂
                              + Ji(i,j,k,e,2,2) * tau(i,j,k,2,2) & ! τ₂₂
                              + Ji(i,j,k,e,2,3) * tau(i,j,k,3,2) )

      w(i,j,k,3) = M(i,j,k) * ( Ji(i,j,k,e,2,1) * tau(i,j,k,3,1) & ! = τ₁₃
                              + Ji(i,j,k,e,2,2) * tau(i,j,k,3,2) & ! = τ₂₃
                              + Ji(i,j,k,e,2,3) * tau(i,j,k,3,3) )
    end do
    end do
    end do

    ! z += [I x (Dˢ)ᵀ x I] w
    do k = 1, np
    do j = 1, np
    do i = 1, np
      do p = 1, np
        z(i,j,k,1) = z(i,j,k,1) + Ds(p,j) * w(i,p,k,1)
        z(i,j,k,2) = z(i,j,k,2) + Ds(p,j) * w(i,p,k,2)
        z(i,j,k,3) = z(i,j,k,3) + Ds(p,j) * w(i,p,k,3)
      end do
    end do
    end do
    end do

    ! w = M ΣᵣJ⁻¹(3,r)⋅τ(r,1:3), exploiting symmetry of τ
    do k = 1, np
    do j = 1, np
    do i = 1, np

      w(i,j,k,1) = M(i,j,k) * ( Ji(i,j,k,e,3,1) * tau(i,j,k,1,1) &
                              + Ji(i,j,k,e,3,2) * tau(i,j,k,2,1) &
                              + Ji(i,j,k,e,3,3) * tau(i,j,k,3,1) )

      w(i,j,k,2) = M(i,j,k) * ( Ji(i,j,k,e,3,1) * tau(i,j,k,2,1) & ! = τ₁₂
                              + Ji(i,j,k,e,3,2) * tau(i,j,k,2,2) & ! τ₂₂
                              + Ji(i,j,k,e,3,3) * tau(i,j,k,3,2) )

      w(i,j,k,3) = M(i,j,k) * ( Ji(i,j,k,e,3,1) * tau(i,j,k,3,1) & ! = τ₁₃
                              + Ji(i,j,k,e,3,2) * tau(i,j,k,3,2) & ! = τ₂₃
                              + Ji(i,j,k,e,3,3) * tau(i,j,k,3,3) )
    end do
    end do
    end do

    ! z += [(Dˢ)ᵀ x I x I] w
    do k = 1, np
    do j = 1, np
    do i = 1, np
      do p = 1, np
        z(i,j,k,1) = z(i,j,k,1) + Ds(p,k) * w(i,j,p,1)
        z(i,j,k,2) = z(i,j,k,2) + Ds(p,k) * w(i,j,p,2)
        z(i,j,k,3) = z(i,j,k,3) + Ds(p,k) * w(i,j,p,3)
      end do
    end do
    end do
    end do

    ! fv = -z
    do k = 1, np
    do j = 1, np
    do i = 1, np
      fv(i,j,k,e,1) = -z(i,j,k,1)
      fv(i,j,k,e,2) = -z(i,j,k,2)
      fv(i,j,k,e,3) = -z(i,j,k,3)
    end do
    end do
    end do

    ! boundary values ..........................................................

    if (get_traces) then

      ! vb = v,  sb = n⋅τ  @ Γ₁ ∪ Γ₂, exploiting symmetry of τ

      i = 1
      do f = 1, 2
        do k = 1, np
        do j = 1, np

          vb(j,k,f,e,1) = v(i,j,k,e,1)
          vb(j,k,f,e,2) = v(i,j,k,e,2)
          vb(j,k,f,e,3) = v(i,j,k,e,3)

          sb(j,k,f,e,1) = n(j,k,f,e,1) * tau(i,j,k,1,1) & ! n₁ τ₁₁
                        + n(j,k,f,e,2) * tau(i,j,k,2,1) & ! n₂ τ₂₁
                        + n(j,k,f,e,3) * tau(i,j,k,3,1)   ! n₃ τ₃₁

          sb(j,k,f,e,2) = n(j,k,f,e,1) * tau(i,j,k,2,1) & ! n₁ τ₁₂
                        + n(j,k,f,e,2) * tau(i,j,k,2,2) & ! n₂ τ₂₂
                        + n(j,k,f,e,3) * tau(i,j,k,3,2)   ! n₃ τ₃₂

          sb(j,k,f,e,3) = n(j,k,f,e,1) * tau(i,j,k,3,1) & ! n₁ τ₁₃
                        + n(j,k,f,e,2) * tau(i,j,k,3,2) & ! n₂ τ₂₃
                        + n(j,k,f,e,3) * tau(i,j,k,3,3)   ! n₃ τ₃₃

        end do
        end do
        i = np
      end do

      ! vb = v,  sb = n⋅τ  @ Γ₃ ∪ Γ₄

      j = 1
      do f = 3, 4
        do k = 1, np
        do i = 1, np

          vb(i,k,f,e,1) = v(i,j,k,e,1)
          vb(i,k,f,e,2) = v(i,j,k,e,2)
          vb(i,k,f,e,3) = v(i,j,k,e,3)

          sb(i,k,f,e,1) = n(i,k,f,e,1) * tau(i,j,k,1,1) & ! n₁ τ₁₁
                        + n(i,k,f,e,2) * tau(i,j,k,2,1) & ! n₂ τ₂₁
                        + n(i,k,f,e,3) * tau(i,j,k,3,1)   ! n₃ τ₃₁

          sb(i,k,f,e,2) = n(i,k,f,e,1) * tau(i,j,k,2,1) & ! n₁ τ₁₂
                        + n(i,k,f,e,2) * tau(i,j,k,2,2) & ! n₂ τ₂₂
                        + n(i,k,f,e,3) * tau(i,j,k,3,2)   ! n₃ τ₃₂

          sb(i,k,f,e,3) = n(i,k,f,e,1) * tau(i,j,k,3,1) & ! n₁ τ₁₃
                        + n(i,k,f,e,2) * tau(i,j,k,3,2) & ! n₂ τ₂₃
                        + n(i,k,f,e,3) * tau(i,j,k,3,3)   ! n₃ τ₃₃
        end do
        end do
        j = np
      end do

      ! vb = v,  sb = n⋅τ  @ Γ₅ ∪ Γ₆

      k = 1
      do f = 5, 6
        do j = 1, np
        do i = 1, np

          vb(i,j,f,e,1) = v(i,j,k,e,1)
          vb(i,j,f,e,2) = v(i,j,k,e,2)
          vb(i,j,f,e,3) = v(i,j,k,e,3)

          sb(i,j,f,e,1) = n(i,j,f,e,1) * tau(i,j,k,1,1) & ! n₁ τ₁₁
                        + n(i,j,f,e,2) * tau(i,j,k,2,1) & ! n₂ τ₂₁
                        + n(i,j,f,e,3) * tau(i,j,k,3,1)   ! n₃ τ₃₁

          sb(i,j,f,e,2) = n(i,j,f,e,1) * tau(i,j,k,2,1) & ! n₁ τ₁₂
                        + n(i,j,f,e,2) * tau(i,j,k,2,2) & ! n₂ τ₂₂
                        + n(i,j,f,e,3) * tau(i,j,k,3,2)   ! n₃ τ₃₂

          sb(i,j,f,e,3) = n(i,j,f,e,1) * tau(i,j,k,3,1) & ! n₁ τ₁₃
                        + n(i,j,f,e,2) * tau(i,j,k,3,2) & ! n₂ τ₂₃
                        + n(i,j,f,e,3) * tau(i,j,k,3,3)   ! n₃ τ₃₃

        end do
        end do
        k = np
      end do

    end if

  end do
  !$omp end do

end subroutine TPO_Stress_DLCI_Gen_RWP
