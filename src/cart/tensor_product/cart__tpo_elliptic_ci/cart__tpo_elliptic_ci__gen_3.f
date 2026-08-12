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

!> summary:  element diffusion operator: 3D Cartesian equidistant, np = 3
!> author:   Joerg Stiller
!> date:     2017/11/14
!>
!>### element diffusion operator: 3D Cartesian equidistant, np = 3
!===============================================================================

subroutine CART__TPO_Diffusion_CI__gen_3(np, ne, Ms, Ls, lambda, nu, dx, u, v)

  !-----------------------------------------------------------------------------
  ! modules

  use Kind_Parameters,  only: RNP
  implicit none

  !-----------------------------------------------------------------------------
  ! arguments

  integer,   intent(in)  :: np          !< number of points per direction
  integer,   intent(in)  :: ne          !< number of elements
  real(RNP), intent(in)  :: Ms(3)       !< 1D standard mass matrix
  real(RNP), intent(in)  :: Ls(3,3)     !< 1D standard stiffness matrix
  real(RNP), intent(in)  :: lambda      !< Helmholtz parameter
  real(RNP), intent(in)  :: nu          !< diffusivity
  real(RNP), intent(in)  :: dx(3)       !< element extensions
  real(RNP), intent(in)  :: u(3,3,3,ne) !< operand
  real(RNP), intent(out) :: v(3,3,3,ne) !< result

  real(RNP), allocatable :: M(:,:,:), M_u(:,:,:), Lm(:,:)
  real(RNP) :: g(3), tmp

  integer :: e, i, j, k

  !-----------------------------------------------------------------------------
  ! initialization

  ! workspace
  allocate(M(np,np,np), M_u(np,np,np), Lm(np,np))

  ! element mass matrix
  tmp = product(dx) / 8
  do k = 1, 3
  do j = 1, 3
  do i = 1, 3
    M(i,j,k) = tmp * Ms(k) * Ms(j) * Ms(i)
  end do
  end do
  end do

  ! mass-weighted stiffness matrix: Lm = Ms^-1 Ls = (Ls Ms^-1)^T
  do j = 1, 3
  do i = 1, 3
    Lm(i,j) = Ls(i,j) / Ms(i)
  end do
  end do

  ! coefficients
  g = 4 * nu / dx**2

  !-----------------------------------------------------------------------------
  ! evaluation

  !$acc data present(u,v) copyin(g,M,Lm) async
  !$acc parallel async
  !$acc loop gang worker private(M_u)

  !$omp do private(e)
  do e = 1, ne

    call SubOp_0(lambda, M, u(:,:,:,e), M_u, v(:,:,:,e))
    call SubOp_1(g(1), Lm, M_u, v(:,:,:,e))
    call SubOp_2(g(2), Lm, M_u, v(:,:,:,e))
    call SubOp_3(g(3), Lm, M_u, v(:,:,:,e))

  end do
  !$omp end do

  !$acc end parallel
  !$acc end data

contains

!-------------------------------------------------------------------------------
!> Projection to mass matrix and application of linear operator

subroutine SubOp_0(lambda, M, u, M_u, v)
  !$acc routine vector
  real(RNP), intent(in)  :: lambda
  real(RNP), intent(in)  :: M(27)
  real(RNP), intent(in)  :: u(27)
  real(RNP), intent(out) :: M_u(27)
  real(RNP), intent(out) :: v(27)

  integer :: l

  !$acc loop vector
  !DIR$ SIMD
  do l = 1, 27
    M_u(l) = M(l) * u(l)
    v(l) = lambda * M_u(l)
  end do

end subroutine SubOp_0

!-------------------------------------------------------------------------------
!> Diffusion operator in direction 1

subroutine SubOp_1(g1, Lm, M_u, v)
  !$acc routine vector
  real(RNP), intent(in)    :: g1
  real(RNP), intent(in)    :: Lm(3,3)
  real(RNP), intent(in)    :: M_u(3,3,3)
  real(RNP), intent(inout) :: v(3,3,3)

  real(RNP) :: tmp1, tmp2, tmp3
  integer   :: j, k, p

  !$acc loop collapse(2) independent vector
  do k = 1, 3
  do j = 1, 3
    tmp1 = 0
    tmp2 = 0
    tmp3 = 0
    do p = 1, 3
      tmp1 = tmp1 + Lm(p,1) * M_u(p,j,k)
      tmp2 = tmp2 + Lm(p,2) * M_u(p,j,k)
      tmp3 = tmp3 + Lm(p,3) * M_u(p,j,k)
    end do
    v(1,j,k) = v(1,j,k) + g1 * tmp1
    v(2,j,k) = v(2,j,k) + g1 * tmp2
    v(3,j,k) = v(3,j,k) + g1 * tmp3
  end do
  end do

end subroutine SubOp_1

!-------------------------------------------------------------------------------
!> Diffusion operator in direction 2

subroutine SubOp_2(g2, Lm, M_u, v)
  !$acc routine vector
  real(RNP), intent(in)    :: g2
  real(RNP), intent(in)    :: Lm(3,3)
  real(RNP), intent(in)    :: M_u(3,3,3)
  real(RNP), intent(inout) :: v(3,3,3)

  real(RNP) :: tmp1, tmp2, tmp3
  integer   :: j, k, p

  !$acc loop collapse(2) independent vector
  do k = 1, 3
  do j = 1, 3
    tmp1 = 0
    tmp2 = 0
    tmp3 = 0
    do p = 1, 3
      tmp1 = tmp1 + Lm(p,j) * M_u(1,p,k)
      tmp2 = tmp2 + Lm(p,j) * M_u(2,p,k)
      tmp3 = tmp3 + Lm(p,j) * M_u(3,p,k)
    end do
    v(1,j,k) = v(1,j,k) + g2 * tmp1
    v(2,j,k) = v(2,j,k) + g2 * tmp2
    v(3,j,k) = v(3,j,k) + g2 * tmp3
  end do
  end do

end subroutine SubOp_2

!-------------------------------------------------------------------------------
!> Diffusion operator in direction 3

subroutine SubOp_3(g3, Lm, M_u, v)
  !$acc routine vector
  real(RNP), intent(in)    :: g3
  real(RNP), intent(in)    :: Lm(3,3)
  real(RNP), intent(in)    :: M_u(9,3)
  real(RNP), intent(inout) :: v(9,3)

  real(RNP) :: tmp
  integer   :: ij, k, p

  !$acc loop collapse(2) independent vector
  do k = 1, 3
  do ij = 1, 9
    tmp = 0
    do p = 1, 3
      tmp = tmp + Lm(p,k) * M_u(ij,p)
    end do
    v(ij,k) = v(ij,k) + g3 * tmp
  end do
  end do

end subroutine SubOp_3

!===============================================================================

end subroutine CART__TPO_Diffusion_CI__gen_3
