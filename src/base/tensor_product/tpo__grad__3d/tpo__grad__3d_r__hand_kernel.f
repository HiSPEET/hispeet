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
!> Parametrized 3d gradient kernel using hand-crafted suboperators (RCLI)

subroutine PROC(TPO_Grad_R_Hand__,_NP_)(ne, Ms, Ds, dx, u, up, v)
  integer,   intent(in)  :: ne                     !< num elements
  real(RWP), intent(in)  :: Ms(_NP_)               !< 1D standard mass matrix
  real(RWP), intent(in)  :: Ds(_NP_,_NP_)          !< standard diff matrix
  real(RWP), intent(in)  :: dx(3)                  !< element extensions
  real(RWP), intent(in)  :: u(_NP_,_NP_,_NP_,ne)   !< 3D scalar field
  real(RWP), intent(in)  :: up(_NP_,_NP_,6,ne)     !< exterior traces u⁺ at faces
  real(RWP), intent(out) :: v(_NP_,_NP_,_NP_,ne,3) !< element-wise gradient of u

  optional :: Ms, up

  !-----------------------------------------------------------------------------
  ! local variables

  real(RWP), parameter :: ZERO = 0
  real(RWP) :: Ds_t(_NP_,_NP_)
  real(RWP) :: g(3), tmp
  integer   :: e, f, i, j, k, n
  logical   :: fluxes

  !-----------------------------------------------------------------------------
  ! initialization

  fluxes = present(Ms) .and. present(up)

  Ds_t = transpose(Ds)

  ! metric coefficients
  g  = 2 / dx

  !---------------------------------------------------------------------------
  ! evaluation

  !$omp do
  do e = 1, ne

    ! v1 = du/dx1
    v(:,:,:,e,1) = 0  ! avoid trouble with NaNs
    call PROC(IxIxQt__,_NP_)(Ds_t, g(1), ZERO, u(:,:,:,e), v(:,:,:,e,1))
    ! v2 = du/dx2
    v(:,:,:,e,2) = 0  ! avoid trouble with NaNs
    call PROC(IxQtxI__,_NP_)(Ds_t, g(2), ZERO, u(:,:,:,e), v(:,:,:,e,2))
    ! v3 = du/dx3
    v(:,:,:,e,3) = 0  ! avoid trouble with NaNs
    call PROC(QtxIxI__,_NP_)(Ds_t, g(3), ZERO, u(:,:,:,e), v(:,:,:,e,3))

    if (fluxes) then

      ! faces 1 and 2
      i =  1
      n = -1
      do f = 1, 2
        tmp = n / (dx(1) * Ms(i))
        do k = 1, _NP_
        do j = 1, _NP_
          v(i,j,k,e,1) = v(i,j,k,e,1) + tmp * (up(j,k,f,e) - u(i,j,k,e))
        end do
        end do
        i = _NP_
        n = 1
      end do

      ! faces 3 and 4
      j =  1
      n = -1
      do f = 3, 4
        tmp = n / (dx(2) * Ms(j))
        do k = 1, _NP_
        do i = 1, _NP_
          v(i,j,k,e,2) = v(i,j,k,e,2) + tmp * (up(i,k,f,e) - u(i,j,k,e))
        end do
        end do
        j = _NP_
        n = 1
      end do

      ! faces 5 and 6
      k =  1
      n = -1
      do f = 5, 6
        tmp = n / (dx(3) * Ms(k))
        do j = 1, _NP_
        do i = 1, _NP_
          v(i,j,k,e,3) = v(i,j,k,e,3) + tmp * (up(i,j,f,e) - u(i,j,k,e))
        end do
        end do
        k = _NP_
        n = 1
      end do

    end if

  end do

  !-----------------------------------------------------------------------------

end subroutine PROC(TPO_Grad_R_Hand__,_NP_)
