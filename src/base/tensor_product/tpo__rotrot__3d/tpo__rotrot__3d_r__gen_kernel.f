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

!> summary:  Generic rotrot of a vector field: 3D Cartesian equidistant
!> author:   Joerg Stiller, Erik Pfister
!> date:     2017/01/19
!===============================================================================

subroutine TPO_RotRot_R_Gen_RWP(np, ne, Ds, dx, u, v)
  integer,   intent(in)  :: np               !< number of points per direction
  integer,   intent(in)  :: ne               !< number of elements
  real(RWP), intent(in)  :: Ds(np,np)        !< 1D standard diff matrix
  real(RWP), intent(in)  :: dx(3)            !< element extensions
  real(RWP), intent(in)  :: u(np,np,np,ne,3) !< 3D vector field
  real(RWP), intent(out) :: v(np,np,np,ne,3) !< v = rot(rot(u))

  !-----------------------------------------------------------------------------
  ! local variables

  real(RWP) :: A(np,np), DA(np,np)
  real(RWP) :: div_u(np,np,np)
  real(RWP) :: g(3), tmp, tmp1, tmp2

  integer :: e, i, j, k, p
  integer :: vec_len

  !-----------------------------------------------------------------------------
  ! initialization

  ! OpenACC vector length
  if (np < 8) then
    vec_len = 128
  else
    vec_len = 256
  end if

  ! Initialization
  A  = transpose(Ds)

  do j = 1, np
  do i = 1, np
    DA(j,i) = sum(Ds(i,:) * Ds(:,j))
  end do
  end do

  ! metric coefficients
  g = 2 / dx

  !-----------------------------------------------------------------------------
  ! evaluation

  !$acc data present(u,v) copyin(A,g) async
  !$acc parallel async &
  !$acc & device_type(nvidia) num_workers(1024/vec_len) vector_length(vec_len)
  !$acc loop gang worker

  !$omp do private(e)
  do e = 1, ne

    ! div_u = du1/dx1 ..........................................................

    !$acc loop collapse(3) vector
    do k = 1, np
    do j = 1, np
    do i = 1, np

      tmp = 0
      do p = 1, np
        tmp = tmp + A(p,i) * u(p,j,k,e,1)
      end do
      div_u(i,j,k) = g(1) * tmp

    end do
    end do
    end do

    ! div_u+= du2/dx2 ..........................................................

    !$acc loop collapse(3) vector
    do k = 1, np
    do j = 1, np
    do i = 1, np

      tmp = 0
      do p = 1, np
        tmp = tmp + A(p,j) * u(i,p,k,e,2)
      end do
      div_u(i,j,k) = div_u(i,j,k) + g(2) * tmp

    end do
    end do
    end do

    ! div_u+= du3/dx3 ..........................................................

    !$acc loop collapse(3) vector
    do k = 1, np
    do j = 1, np
    do i = 1, np

      tmp = 0
      do p = 1, np
        tmp = tmp + A(p,k) * u(i,j,p,e,3)
      end do
      div_u(i,j,k) = div_u(i,j,k) + g(3) * tmp

    end do
    end do
    end do

    ! v1 = d(div_u)/dx1 ........................................................

    !$acc loop collapse(3) vector
    do k = 1, np
    do j = 1, np
    do i = 1, np

      tmp = 0
      do p = 1, np
        tmp = tmp + A(p,i) * div_u(p,j,k)
      end do
      v(i,j,k,e,1) = g(1) * tmp

    end do
    end do
    end do

    ! v1+=-d( du1/dx1 )/dx1 - d( du1/dx2 )/dx2 .................................

    !$acc loop collapse(3) vector
    do k = 1, np
    do j = 1, np
    do i = 1, np

      tmp1 = 0
      tmp2 = 0
      do p = 1, np
        tmp1 = tmp1 + DA(p,i) * u(p,j,k,e,1)
        tmp2 = tmp2 + DA(p,j) * u(i,p,k,e,1)
      end do
      v(i,j,k,e,1) = v(i,j,k,e,1) - (g(1)*g(1)*tmp1 + g(2)*g(2)*tmp2)

    end do
    end do
    end do

    ! v1+=-d( du1/dx3 )/dx3 ....................................................

    !$acc loop collapse(3) vector
    do k = 1, np
    do j = 1, np
    do i = 1, np

      tmp = 0
      do p = 1, np
        tmp = tmp + DA(p,k) * u(i,j,p,e,1)
      end do
      v(i,j,k,e,1) = v(i,j,k,e,1) - g(3)*g(3)*tmp

    end do
    end do
    end do

    ! v2 = d(div_u)/dx2 ........................................................

    !$acc loop collapse(3) vector
    do k = 1, np
    do j = 1, np
    do i = 1, np

      tmp = 0
      do p = 1, np
        tmp = tmp + A(p,j) * div_u(i,p,k)
      end do
      v(i,j,k,e,2) = g(2) * tmp

    end do
    end do
    end do

    ! v2+=-d( du2/dx1 )/dx1 - d( du2/dx2 )/dx2 ..................................

    !$acc loop collapse(3) vector
    do k = 1, np
    do j = 1, np
    do i = 1, np

      tmp1 = 0
      tmp2 = 0
      do p = 1, np
        tmp1 = tmp1 + DA(p,i) * u(p,j,k,e,2)
        tmp2 = tmp2 + DA(p,j) * u(i,p,k,e,2)
      end do
      v(i,j,k,e,2) = v(i,j,k,e,2) - (g(1)*g(1)*tmp1 + g(2)*g(2)*tmp2)

    end do
    end do
    end do

    ! v2+=-d( du2/dx3 )/dx3 .....................................................

    !$acc loop collapse(3) vector
    do k = 1, np
    do j = 1, np
    do i = 1, np

      tmp = 0
      do p = 1, np
        tmp = tmp + DA(p,k) * u(i,j,p,e,2)
      end do
      v(i,j,k,e,2) = v(i,j,k,e,2) - g(3)*g(3)*tmp

    end do
    end do
    end do

    ! v3 = d(div_u)/dx3 ........................................................

    !$acc loop collapse(3) vector
    do k = 1, np
    do j = 1, np
    do i = 1, np

      tmp = 0
      do p = 1, np
        tmp = tmp + A(p,k) * div_u(i,j,p)
      end do
      v(i,j,k,e,3) = g(3) * tmp

    end do
    end do
    end do

    ! v3+=-d( du3/dx1 )/dx1 - d( du3/dx2 )/dx2 ..................................

    !$acc loop collapse(3) vector
    do k = 1, np
    do j = 1, np
    do i = 1, np

      tmp1 = 0
      tmp2 = 0
      do p = 1, np
        tmp1 = tmp1 + DA(p,i) * u(p,j,k,e,3)
        tmp2 = tmp2 + DA(p,j) * u(i,p,k,e,3)
      end do
      v(i,j,k,e,3) = v(i,j,k,e,3) - (g(1)*g(1)*tmp1 + g(2)*g(2)*tmp2)

    end do
    end do
    end do

    ! v3+=-d( du3/dx3 )/dx3 ....................................................

    !$acc loop collapse(3) vector
    do k = 1, np
    do j = 1, np
    do i = 1, np

      tmp = 0
      do p = 1, np
        tmp = tmp + DA(p,k) * u(i,j,p,e,3)
      end do
      v(i,j,k,e,3) = v(i,j,k,e,3) - g(3)*g(3)*tmp

    end do
    end do
    end do

  end do
  !$omp end do

  !$acc end parallel
  !$acc end data

!===============================================================================

end subroutine TPO_RotRot_R_Gen_RWP
