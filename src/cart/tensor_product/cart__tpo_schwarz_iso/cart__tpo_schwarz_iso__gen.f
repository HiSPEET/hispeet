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

!> summary:  Isotropic Schwarz operator: bf-marching, separate, generic
!> author:   Joerg Stiller
!> date:     2017/10/16
!>
!>### Isotropic Schwarz operator: bf-marching, separate, generic
!===============================================================================

subroutine CART__TPO_Schwarz_Iso__gen(ns, nc, nd, S, W, cfg, D_inv, f, u)

  !-----------------------------------------------------------------------------
  ! modules

  use Kind_Parameters, only: RNP
  implicit none

  !-----------------------------------------------------------------------------
  ! arguments

  integer,   intent(in)  :: ns                 !< dimension of S
  integer,   intent(in)  :: nc                 !< number of configurations
  integer,   intent(in)  :: nd                 !< number of subdomains
  real(RNP), intent(in)  :: S(ns,ns,nc)        !< 1D eigenvectors
  real(RNP), intent(in)  :: W(ns,nc)           !< 1D weights
  integer,   intent(in)  :: cfg(3,nd)          !< subdomain configurations
  real(RNP), intent(in)  :: D_inv(ns,ns,ns,nd) !< inverse 3D eigenvalues
  real(RNP), intent(in)  :: f(ns,ns,ns,nd)     !< RHS
  real(RNP), intent(out) :: u(ns,ns,ns,nd)     !< solution

  !-----------------------------------------------------------------------------
  ! local variables

  real(RNP), allocatable :: WS_t(:,:,:), y(:,:,:), z(:,:,:)
  real(RNP) :: tmp

  integer :: i, j, k, l, m
  integer :: c, c1, c2, c3
  integer :: vec_len

  !-----------------------------------------------------------------------------
  ! initialization

  ! OpenACC vector length
  if (ns < 8) then
    vec_len = 128
  else
    vec_len = 256
  end if

  allocate(WS_t(ns,ns,nc), y(ns,ns,ns), z(ns,ns,ns))

  do c = 1, nc
  do j = 1, ns
  do i = 1, ns
    WS_t(j,i,c) = W(i,c) * S(i,j,c)
  end do
  end do
  end do

  !-----------------------------------------------------------------------------
  ! evaluation

  !$acc data async present(cfg, D_inv, f, u) copyin(S, WS_t)
  !$acc parallel async &
  !$acc & device_type(nvidia) num_workers(1024/vec_len) vector_length(vec_len)
  !$acc loop gang worker private(y,z)

  !$omp do private(l)
  subdomains: do l = 1, nd

    ! element-boundary configurations
    c1 = cfg(1,l)
    c2 = cfg(2,l)
    c3 = cfg(3,l)

    ! y = S^t x I x I f ........................................................

    !$acc loop collapse(3) vector
    do k = 1, ns
    do j = 1, ns
    do i = 1, ns
      tmp = 0
      do m = 1, ns
        tmp = tmp + S(m,k,c3) * f(i,j,m,l)
      end do
      y(i,j,k) = tmp
    end do
    end do
    end do

    ! z = I x S^t x I y ........................................................

    !$acc loop collapse(3) vector
    do k = 1, ns
    do j = 1, ns
    do i = 1, ns
      tmp = 0
      do m = 1, ns
        tmp = tmp + S(m,j,c2) * y(i,m,k)
      end do
      z(i,j,k) = tmp
    end do
    end do
    end do

    ! y = D⁻¹ (I x I x S^t) z ..................................................

    !$acc loop collapse(3) vector
    do k = 1, ns
    do j = 1, ns
    do i = 1, ns
      tmp = 0
      do m = 1, ns
        tmp = tmp + S(m,i,c1) * z(m,j,k)
      end do
      y(i,j,k) = tmp * D_inv(i,j,k,l)
    end do
    end do
    end do

    ! z = I x I x WS y .........................................................

    !$acc loop collapse(3) vector
    do k = 1, ns
    do j = 1, ns
    do i = 1, ns
      tmp = 0
      do m = 1, ns
        tmp = tmp + WS_t(m,i,c1) * y(m,j,k)
      end do
      z(i,j,k) = tmp
    end do
    end do
    end do

    ! y = I x WS x I z .........................................................

    !$acc loop collapse(3) vector
    do k = 1, ns
    do j = 1, ns
    do i = 1, ns
      tmp = 0
      do m = 1, ns
        tmp = tmp + WS_t(m,j,c2) * z(i,m,k)
      end do
      y(i,j,k) = tmp
    end do
    end do
    end do

    ! u = WS x I x I y .........................................................

    !$acc loop collapse(3) vector
    do k = 1, ns
    do j = 1, ns
    do i = 1, ns
      tmp = 0
      do m = 1, ns
        tmp = tmp + WS_t(m,k,c3) * y(i,j,m)
      end do
      u(i,j,k,l) = tmp
    end do
    end do
    end do

  end do subdomains
  !$omp end do

  !$acc end parallel
  !$acc end data

!===============================================================================

end subroutine CART__TPO_Schwarz_Iso__gen
