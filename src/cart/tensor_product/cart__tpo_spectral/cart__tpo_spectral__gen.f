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

!> summary:  Spectral operator backward-forward marching, separate, generic
!> author:   Joerg Stiller
!> date:     2017/05/18
!>
!>### Spectral operator backward-forward marching, separate, generic
!===============================================================================

subroutine CART__TPO_Spectral__gen(na, nb, nc, ne, A, B, C, Lambda, u, v)

  !-----------------------------------------------------------------------------
  ! modules

  use Kind_Parameters, only: RNP
  implicit none

  !-----------------------------------------------------------------------------
  ! arguments

  integer,   intent(in)  :: na               !< dimension of A
  integer,   intent(in)  :: nb               !< dimension of B
  integer,   intent(in)  :: nc               !< dimension of C
  integer,   intent(in)  :: ne               !< number of elements
  real(RNP), intent(in)  :: A(na,na)         !< eigenvectors for direction 1
  real(RNP), intent(in)  :: B(nb,nb)         !< eigenvectors for direction 2
  real(RNP), intent(in)  :: C(nc,nc)         !< eigenvectors for direction 3
  real(RNP), intent(in)  :: Lambda(na,nb,nc) !< eigenvalues
  real(RNP), intent(in)  :: u(na,nb,nc,ne)   !< operand
  real(RNP), intent(out) :: v(na,nb,nc,ne)   !< result

  !-----------------------------------------------------------------------------
  ! local variables

  real(RNP), allocatable :: z1(:,:,:), z2(:,:,:), z3(:,:,:)
  real(RNP) :: tmp

  integer :: e, i, j, k, l
  integer :: vec_len

  !-----------------------------------------------------------------------------
  ! initialization

  ! OpenACC vector length
  if (min(na,nb,nc) < 8) then
    vec_len = 128
  else
    vec_len = 256
  end if

  allocate(z1(na,nb,nc), z2(na,nb,nc), z3(na,nb,nc))

  !-----------------------------------------------------------------------------
  ! evaluation

  !$acc data present(u,v) copyin(A,B,C,Lambda) async
  !$acc parallel async &
  !$acc & device_type(nvidia) num_workers(1024/vec_len) vector_length(vec_len)
  !$acc loop gang worker private(z1,z2,z3)

  !$omp do private(e)
  elements: do e = 1, ne

    ! apply C^t ................................................................

    !$acc loop collapse(3) vector
    do k = 1, nc
    do j = 1, nb
    do i = 1, na
      tmp = 0
      do l = 1, nc
        tmp = tmp + C(l,k) * u(i,j,l,e)
      end do
      z3(i,j,k) = tmp
    end do
    end do
    end do

    ! apply B^t ................................................................

    !$acc loop vector collapse(3)
    do k = 1, nc
    do j = 1, nb
    do i = 1, na
      tmp = 0
      do l = 1, nb
        tmp = tmp + B(l,j) * z3(i,l,k)
      end do
      z2(i,j,k) = tmp
    end do
    end do
    end do

    ! apply Lambda A^t ..............................................................

    !$acc loop vector collapse(3)
    do k = 1, nc
    do j = 1, nb
    do i = 1, na
      tmp = 0
      do l = 1, na
        tmp = tmp + A(l,i) * z2(l,j,k)
      end do
      z1(i,j,k) = tmp * Lambda(i,j,k)
    end do
    end do
    end do

    ! apply A ..................................................................

    !$acc loop vector collapse(3)
    do k = 1, nc
    do j = 1, nb
    do i = 1, na
      tmp = 0
      do l = 1, na
        tmp = tmp + A(i,l) * z1(l,j,k)
      end do
      z2(i,j,k) = tmp
    end do
    end do
    end do

    ! apply B ..................................................................

    !$acc loop vector collapse(3)
    do k = 1, nc
    do j = 1, nb
    do i = 1, na
      tmp = 0
      do l = 1, nb
        tmp = tmp + B(j,l) * z2(i,l,k)
      end do
      z3(i,j,k) = tmp
    end do
    end do
    end do

    ! apply C ..................................................................

    !$acc loop collapse(3) vector
    do k = 1, nc
    do j = 1, nb
    do i = 1, na
      tmp = 0
      do l = 1, nc
        tmp = tmp + C(k,l) * z3(i,j,l)
      end do
      v(i,j,k,e) = tmp
    end do
    end do
    end do

  end do elements
  !$omp end do

  !$acc end parallel
  !$acc end data

!===============================================================================

end subroutine  CART__TPO_Spectral__gen
