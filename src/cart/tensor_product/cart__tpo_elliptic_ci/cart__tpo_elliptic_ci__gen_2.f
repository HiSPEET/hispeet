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

!> summary:  element diffusion operator: 3D Cartesian equidistant, np = 2
!> author:   Joerg Stiller
!> date:     2017/11/14
!>
!>### element diffusion operator: 3D Cartesian equidistant, np = 2
!===============================================================================

subroutine CART__TPO_Diffusion_CI__gen_2(np, ne, Ms, Ls, lambda, nu, dx, u, v)

  !-----------------------------------------------------------------------------
  ! modules

  use Kind_Parameters,  only: RNP
  implicit none

  !-----------------------------------------------------------------------------
  ! arguments

  integer,   intent(in)  :: np          !< number of points per direction
  integer,   intent(in)  :: ne          !< number of elements
  real(RNP), intent(in)  :: Ms(2)       !< 1D standard mass matrix
  real(RNP), intent(in)  :: Ls(2,2)     !< 1D standard stiffness matrix
  real(RNP), intent(in)  :: lambda      !< Helmholtz parameter
  real(RNP), intent(in)  :: nu          !< diffusivity
  real(RNP), intent(in)  :: dx(3)       !< element extensions
  real(RNP), intent(in)  :: u(2,2,2,ne) !< operand
  real(RNP), intent(out) :: v(2,2,2,ne) !< result

  real(RNP), allocatable :: M(:,:,:), M_u(:,:,:), Lm(:,:)
  real(RNP) :: g(3), tmp

  integer :: e, i, j, k

  !-----------------------------------------------------------------------------
  ! initialization

  ! workspace
  allocate(M(np,np,np), M_u(np,np,np), Lm(np,np))

  ! element mass matrix
  tmp = product(dx) / 8
  do k = 1, 2
  do j = 1, 2
  do i = 1, 2
    M(i,j,k) = tmp * Ms(k) * Ms(j) * Ms(i)
  end do
  end do
  end do

  ! mass-weighted stiffness matrix: Lm = Ms^-1 Ls = (Ls Ms^-1)^T
  do j = 1, 2
  do i = 1, 2
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

    call Projection(M, u(:,:,:,e), M_u)
    call Diffusion(lambda, g, Lm, M_u, v(:,:,:,e))

  end do
  !$omp end do

  !$acc end parallel
  !$acc end data

contains

!-------------------------------------------------------------------------------
!> Projection to element mass matrix

subroutine Projection(M, u, M_u)
  !$acc routine vector
  real(RNP), intent(in)  :: M(8)
  real(RNP), intent(in)  :: u(8)
  real(RNP), intent(out) :: M_u(8)

  integer :: l

  !$acc loop vector
  do l = 1, 8
    M_u(l) = M(l) * u(l)
  end do

end subroutine Projection

!-------------------------------------------------------------------------------
!> Element operator, including linear part and diffusion in all directions

subroutine Diffusion(lambda, g, Lm, M_u, v)
  !$acc routine vector
  real(RNP), intent(in)  :: lambda
  real(RNP), intent(in)  :: g(3)
  real(RNP), intent(in)  :: Lm(2,2)
  real(RNP), intent(in)  :: M_u(2,2,2)
  real(RNP), intent(out) :: v(2,2,2)

  integer :: j, k

  !$acc loop collapse(2) independent vector
  do k = 1, 2
  do j = 1, 2

    v(1,j,k) = lambda * M_u(1,j,k)                                   &
             + g(1) * (Lm(1,1) * M_u(1,j,k) + Lm(2,1) * M_u(2,j,k))  &
             + g(2) * (Lm(1,j) * M_u(1,1,k) + Lm(2,j) * M_u(1,2,k))  &
             + g(3) * (Lm(1,k) * M_u(1,j,1) + Lm(2,k) * M_u(1,j,2))

    v(2,j,k) = lambda * M_u(2,j,k)                                   &
             + g(1) * (Lm(1,2) * M_u(1,j,k) + Lm(2,2) * M_u(2,j,k))  &
             + g(2) * (Lm(1,j) * M_u(2,1,k) + Lm(2,j) * M_u(2,2,k))  &
             + g(3) * (Lm(1,k) * M_u(2,j,1) + Lm(2,k) * M_u(2,j,2))
  end do
  end do

end subroutine Diffusion

!===============================================================================

end subroutine CART__TPO_Diffusion_CI__gen_2
