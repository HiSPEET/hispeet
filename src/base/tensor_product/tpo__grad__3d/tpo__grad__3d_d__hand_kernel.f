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

!-------------------------------------------------------------------------------
!> Parametrized 3d gradient kernel using hand-crafted suboperators (D)

subroutine PROC(TPO_Grad_D_Hand__,_NP_)(ne, Ms, Ds, Jd, Ji, a, n, u, up, v)
  integer,   intent(in)  :: ne                        !< num elements
  real(RWP), intent(in)  :: Ms(_NP_)                  !< 1D standard mass matrix
  real(RWP), intent(in)  :: Ds(_NP_,_NP_)             !< standard diff matrix
  real(RWP), intent(in)  :: Jd(_NP_,_NP_,_NP_,ne)     !< element Jacobian determinant
  real(RWP), intent(in)  :: Ji(_NP_,_NP_,_NP_,ne,3,3) !< element extensions
  real(RWP), intent(in)  :: a(_NP_,_NP_,6,ne)         !< face area coefficients
  real(RWP), intent(in)  :: n(_NP_,_NP_,6,ne,3)       !< face unit normal vectors
  real(RWP), intent(in)  :: u(_NP_,_NP_,_NP_,ne)      !< 3D scalar field
  real(RWP), intent(in)  :: up(_NP_,_NP_,6,ne)        !< exterior traces u⁺ at faces
  real(RWP), intent(out) :: v(_NP_,_NP_,_NP_,ne,3)    !< element-wise gradient of u

  optional :: Ms, Jd, a, n, up

  !-----------------------------------------------------------------------------
  ! local variables

  real(RWP), parameter :: ZERO = 0, HALF = 0.5, ONE = 1
  real(RWP) :: r(_NP_,_NP_,_NP_)
  real(RWP) :: s(_NP_,_NP_,_NP_)
  real(RWP) :: t(_NP_,_NP_,_NP_)
  real(RWP) :: Ds_t(_NP_,_NP_)
  real(RWP) :: c

  integer :: e, f, i, j, k
  logical :: fluxes

  !-----------------------------------------------------------------------------
  ! initialization

  fluxes = present(Ms) .and. &
           present(Jd) .and. &
           present(a)  .and. &
           present(n)  .and. &
           present(up)

  Ds_t = transpose(Ds)

  r = 0  ! avoid trouble with NaNs and memory warnings
  s = 0  ! avoid trouble with NaNs and memory warnings
  t = 0  ! avoid trouble with NaNs and memory warnings

  !-----------------------------------------------------------------------------
  ! evaluation

  !$omp do
  do e = 1, ne

    ! [r,s,t]ᵀ = ∇ˢu ...........................................................

    call PROC(IxIxQt__,_NP_)(Ds_t, ONE, ZERO, u(:,:,:,e), r)
    call PROC(IxQtxI__,_NP_)(Ds_t, ONE, ZERO, u(:,:,:,e), s)
    call PROC(QtxIxI__,_NP_)(Ds_t, ONE, ZERO, u(:,:,:,e), t)

    ! v = J⁻ᵀ ⋅ [r,s,t]^T ......................................................

    v(:,:,:,e,1) = r * Ji(:,:,:,e,1,1) + &
                   s * Ji(:,:,:,e,2,1) + &
                   t * Ji(:,:,:,e,3,1)

    v(:,:,:,e,2) = r * Ji(:,:,:,e,1,2) + &
                   s * Ji(:,:,:,e,2,2) + &
                   t * Ji(:,:,:,e,3,2)

    v(:,:,:,e,3) = r * Ji(:,:,:,e,1,3) + &
                   s * Ji(:,:,:,e,2,3) + &
                   t * Ji(:,:,:,e,3,3)

    if (fluxes) then

      ! faces 1 and 2
      do f = 1, 2
        i = 1 + (_NP_ - 1) * (f - 1)
        do k = 1, _NP_
        do j = 1, _NP_
          c = a(j,k,f,e) / (2 * Ms(i ) * Jd(i,j,k,e)) * (up(j,k,f,e) - u(i,j,k,e))
          v(i,j,k,e,1) = v(i,j,k,e,1) + c * n(j,k,f,e,1)
          v(i,j,k,e,2) = v(i,j,k,e,2) + c * n(j,k,f,e,2)
          v(i,j,k,e,3) = v(i,j,k,e,3) + c * n(j,k,f,e,3)
        end do
        end do
      end do

      ! faces 3 and 4
      do f = 3, 4
        j = 1 + (_NP_ - 1) * (f - 3)
        do k = 1, _NP_
        do i = 1, _NP_
          c = a(i,k,f,e) / (2 * Ms(j) * Jd(i,j,k,e)) * (up(i,k,f,e) - u(i,j,k,e))
          v(i,j,k,e,1) = v(i,j,k,e,1) + c * n(i,k,f,e,1)
          v(i,j,k,e,2) = v(i,j,k,e,2) + c * n(i,k,f,e,2)
          v(i,j,k,e,3) = v(i,j,k,e,3) + c * n(i,k,f,e,3)
        end do
        end do
      end do

      ! faces 5 and 6
      do f = 5, 6
        k = 1 + (_NP_ - 1) * (f - 5)
        do j = 1, _NP_
        do i = 1, _NP_
          c = a(i,j,f,e) / (2 * Ms(k) * Jd(i,j,k,e)) * (up(i,j,f,e) - u(i,j,k,e))
          v(i,j,k,e,1) = v(i,j,k,e,1) + c * n(i,j,f,e,1)
          v(i,j,k,e,2) = v(i,j,k,e,2) + c * n(i,j,f,e,2)
          v(i,j,k,e,3) = v(i,j,k,e,3) + c * n(i,j,f,e,3)
        end do
        end do
      end do

    end if

  end do

  !-----------------------------------------------------------------------------

end subroutine PROC(TPO_Grad_D_Hand__,_NP_)
