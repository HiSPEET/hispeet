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

!> summary:  Generic divergence of a vector field: 3D Cartesian equidistant
!> author:   Joerg Stiller
!> date:     2017/01/19
!>
!>### Generic divergence of a vector field: 3D Cartesian equidistant
!===============================================================================

subroutine CART__TPO_Div__gen(np, ne, Ds, dx, u, v, transp)

  !-----------------------------------------------------------------------------
  ! modules

  use Kind_Parameters, only: RNP
  implicit none

  !-----------------------------------------------------------------------------
  ! arguments

  integer,   intent(in)  :: np               !< number of points per direction
  integer,   intent(in)  :: ne               !< number of elements
  real(RNP), intent(in)  :: Ds(np,np)        !< 1D standard diff matrix
  real(RNP), intent(in)  :: dx(3)            !< element extensions
  real(RNP), intent(in)  :: u(np,np,np,ne,3) !< 3D vector field
  real(RNP), intent(out) :: v(np,np,np,ne)   !< element-wise divergence of u
  logical, optional, intent(in) :: transp    !< if present, apply Ds^T

  !-----------------------------------------------------------------------------
  ! local variables

  real(RNP), allocatable :: A(:,:)
  real(RNP) :: g(3), tmp

  integer :: e, i, j, k, m
  integer :: vec_len

  !-----------------------------------------------------------------------------
  ! initialization

  ! OpenACC vector length
  if (np < 8) then
    vec_len = 128
  else
    vec_len = 256
  end if

  ! workspace
  allocate(A(np,np))

  ! adjust diff matrix for optimal memory access
  if (present(transp)) then
    do j = 1, np
    do i = 1, np
      A(i,j) = Ds(i,j)
    end do
    end do
  else
    do j = 1, np
    do i = 1, np
      A(j,i) = Ds(i,j)
    end do
    end do
  end if

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

    ! v = du1/dx1 ..............................................................

    !$acc loop collapse(3) vector
    do k = 1, np
    do j = 1, np
    do i = 1, np

      tmp = 0
      do m = 1, np
        tmp = tmp + A(m,i) * u(m,j,k,e,1)
      end do
      v(i,j,k,e) = g(1) * tmp

    end do
    end do
    end do

    ! v += du2/dx2 .............................................................

    !$acc loop collapse(3) vector
    do k = 1, np
    do j = 1, np
    do i = 1, np

      tmp = 0
      do m = 1, np
        tmp = tmp + A(m,j) * u(i,m,k,e,2)
      end do
      v(i,j,k,e) = v(i,j,k,e) + g(2) * tmp

    end do
    end do
    end do

    ! v += du3/dx3 .............................................................

    !$acc loop collapse(3) vector
    do k = 1, np
    do j = 1, np
    do i = 1, np

      tmp = 0
      do m = 1, np
        tmp = tmp + A(m,k) * u(i,j,m,e,3)
      end do
      v(i,j,k,e) = v(i,j,k,e) + g(3) * tmp

    end do
    end do
    end do

  end do
  !$omp end do

  !$acc end parallel
  !$acc end data

!===============================================================================

end subroutine CART__TPO_Div__gen
