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

!> summary:  Generic gradient of a vector field: 3D Cartesian curvilinear
!> author:   Jerome Michel, Joerg Stiller
!> date:     2021/08/02
!===============================================================================

subroutine TPO_Grad_D_Gen_RWP(np, ne, Ms, Ds, Jd, Ji, a, n, u, up, v)
  integer,   intent(in)  :: np                  !< num points per direction
  integer,   intent(in)  :: ne                  !< num elements
  real(RWP), intent(in)  :: Ms(np)              !< 1D standard mass matrix
  real(RWP), intent(in)  :: Ds(np,np)           !< 1D standard diff matrix
  real(RWP), intent(in)  :: Jd(np,np,np,ne)     !< element Jacobian determinant
  real(RWP), intent(in)  :: Ji(np,np,np,ne,3,3) !< inverse Jacobi matrix
  real(RWP), intent(in)  :: a(np,np,6,ne)       !< face area coefficients
  real(RWP), intent(in)  :: n(np,np,6,ne,3)     !< face unit normal vectors
  real(RWP), intent(in)  :: u(np,np,np,ne)      !< 3D scalar field u
  real(RWP), intent(in)  :: up(np,np,6,ne)      !< exterior traces u⁺ at faces
  real(RWP), intent(out) :: v(np,np,np,ne,3)    !< element-wise gradient of u

  optional :: Ms, Jd, a, n, up

  !-----------------------------------------------------------------------------
  ! local variables

  real(RWP) :: c, r, s, t
  integer   :: e, f, i, j, k, p
  logical   :: fluxes

  !-----------------------------------------------------------------------------
  ! initialization

  fluxes = present(Ms) .and. &
           present(Jd) .and. &
           present(a)  .and. &
           present(n)  .and. &
           present(up)

  !-----------------------------------------------------------------------------
  ! evaluation

  !$omp do private(e)
  do e = 1, ne

    do k = 1, np
    do j = 1, np
    do i = 1, np

      ! [r,s,t]ᵀ = ∇ˢu .........................................................

      r = 0
      s = 0
      t = 0

      do p = 1, np
        r = r + Ds(i,p) * u(p,j,k,e)
        s = s + Ds(j,p) * u(i,p,k,e)
        t = t + Ds(k,p) * u(i,j,p,e)
      end do

      ! v = J⁻ᵀ ⋅ [r,s,t]^T ....................................................

      v(i,j,k,e,1) = r * Ji(i,j,k,e,1,1) &
                   + s * Ji(i,j,k,e,2,1) &
                   + t * Ji(i,j,k,e,3,1)

      v(i,j,k,e,2) = r * Ji(i,j,k,e,1,2) &
                   + s * Ji(i,j,k,e,2,2) &
                   + t * Ji(i,j,k,e,3,2)

      v(i,j,k,e,3) = r * Ji(i,j,k,e,1,3) &
                   + s * Ji(i,j,k,e,2,3) &
                   + t * Ji(i,j,k,e,3,3)

    end do
    end do
    end do

    if (fluxes) then

      ! faces 1 and 2
      do f = 1, 2
        i = 1 + (np - 1) * (f - 1)
        do k = 1, np
        do j = 1, np
          c = a(j,k,f,e) / (2 * Ms(i ) * Jd(i,j,k,e)) * (up(j,k,f,e) - u(i,j,k,e))
          v(i,j,k,e,1) = v(i,j,k,e,1) + c * n(j,k,f,e,1)
          v(i,j,k,e,2) = v(i,j,k,e,2) + c * n(j,k,f,e,2)
          v(i,j,k,e,3) = v(i,j,k,e,3) + c * n(j,k,f,e,3)
        end do
        end do
      end do

      ! faces 3 and 4
      do f = 3, 4
        j = 1 + (np - 1) * (f - 3)
        do k = 1, np
        do i = 1, np
          c = a(i,k,f,e) / (2 * Ms(j) * Jd(i,j,k,e)) * (up(i,k,f,e) - u(i,j,k,e))
          v(i,j,k,e,1) = v(i,j,k,e,1) + c * n(i,k,f,e,1)
          v(i,j,k,e,2) = v(i,j,k,e,2) + c * n(i,k,f,e,2)
          v(i,j,k,e,3) = v(i,j,k,e,3) + c * n(i,k,f,e,3)
        end do
        end do
      end do

      ! faces 5 and 6
      do f = 5, 6
        k = 1 + (np - 1) * (f - 5)
        do j = 1, np
        do i = 1, np
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

end subroutine TPO_Grad_D_Gen_RWP
