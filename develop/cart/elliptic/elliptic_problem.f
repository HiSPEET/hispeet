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

!> summary:  Interface to elliptic problems
!> author:   Joerg Stiller
!> date:     2019/01/27
!===============================================================================

module Elliptic_Problem
  use Kind_Parameters, only: RNP

  implicit none
  private

  public :: EllipticProblem

  !-----------------------------------------------------------------------------
  !> Abstract type for defining and handling an elliptic problem

  type, abstract :: EllipticProblem

    real(RNP) :: lambda = 0  !< Helmholtz parameter
    real(RNP) :: nu_0   = 1  !< diffusivity mean value ν₀
    real(RNP) :: nu_1   = 0  !< diffusivity fluctuation amplitude ν₁
    real(RNP) :: d_nu   = 0  !< diffusivity phase shift
    integer   :: k_nu   = 1  !< diffusivity wave number
    integer   :: k_u    = 1  !< solution wave number

  contains

    procedure :: SetProblem
    procedure :: GetSource

    procedure(GetExactSolution),       deferred :: GetExactSolution
    procedure(GetExactGradient),       deferred :: GetExactGradient
    procedure(GetExactLaplacian),      deferred :: GetExactLaplacian
    procedure(GetDiffusivity),         deferred :: GetDiffusivity
    procedure(GetDiffusivityGradient), deferred :: GetDiffusivityGradient

  end type EllipticProblem

  abstract interface

    !---------------------------------------------------------------------------
    !> Exact solution

    subroutine GetExactSolution(problem, x, u)
      import
      class(EllipticProblem), intent(in)  :: problem
      real(RNP),              intent(in)  :: x(:,:,:,:,:) !< mesh points
      real(RNP),              intent(out) :: u(:,:,:,:)   !< solution, u(x)
    end subroutine GetExactSolution

    !---------------------------------------------------------------------------
    !> Exact solution gradient

    subroutine GetExactGradient(problem, x, grad_u)
      import
      class(EllipticProblem), intent(in)  :: problem
      real(RNP),              intent(in)  :: x(:,:,:,:,:)      !< mesh points
      real(RNP),              intent(out) :: grad_u(:,:,:,:,:) !< ∇u(x)
    end subroutine GetExactGradient

    !---------------------------------------------------------------------------
    !> Exact laplacian

    subroutine GetExactLaplacian(problem, x, laplace_u)
      import
      class(EllipticProblem), intent(in)  :: problem
      real(RNP),              intent(in)  :: x(:,:,:,:,:)       !< mesh points
      real(RNP),              intent(out) :: laplace_u(:,:,:,:) !< ∇²u(x)
    end subroutine GetExactLaplacian

    !---------------------------------------------------------------------------
    !> Diffusivity

    subroutine GetDiffusivity(problem, x, nu)
      import
      class(EllipticProblem), intent(in)  :: problem
      real(RNP),              intent(in)  :: x(:,:,:,:,:)  !< mesh points
      real(RNP),              intent(out) :: nu(:,:,:,:)   !< ν(x)
    end subroutine GetDiffusivity

    !---------------------------------------------------------------------------
    !> Diffusivity gradient

    subroutine GetDiffusivityGradient(problem, x, grad_nu)
      import
      class(EllipticProblem), intent(in)  :: problem
      real(RNP),              intent(in)  :: x(:,:,:,:,:)       !< mesh points
      real(RNP),              intent(out) :: grad_nu(:,:,:,:,:) !< ∇ν(x)
    end subroutine GetDiffusivityGradient

  end interface

contains

!-------------------------------------------------------------------------------
!> Set problem

subroutine SetProblem(problem, lambda, nu_0, nu_1, d_nu, k_nu, k_u)
  class(EllipticProblem), intent(inout) :: problem
  real(RNP), optional :: lambda  !< Helmholtz parameter
  real(RNP), optional :: nu_0    !< diffusivity mean value ν₀
  real(RNP), optional :: nu_1    !< diffusivity fluctuation amplitude ν₁
  real(RNP), optional :: d_nu    !< diffusivity phase shift
  integer  , optional :: k_nu    !< diffusivity wave number
  integer  , optional :: k_u     !< solution wave number

  if (present(lambda))  problem % lambda = lambda
  if (present(nu_0  ))  problem % nu_0   = nu_0
  if (present(nu_1  ))  problem % nu_1   = nu_1
  if (present(d_nu  ))  problem % d_nu   = d_nu
  if (present(k_nu  ))  problem % k_nu   = k_nu
  if (present(k_u   ))  problem % k_u    = k_u

end subroutine SetProblem

!-------------------------------------------------------------------------------
!> Exact source

subroutine GetSource(problem, x, f)
  class(EllipticProblem), intent(in)  :: problem
  real(RNP),              intent(in)  :: x(:,:,:,:,:)  !< mesh points
  real(RNP),              intent(out) :: f(:,:,:,:)    !< source, f(x)

  real(RNP), dimension(:,:,:,:),   allocatable, save :: nu, u, laplace_u
  real(RNP), dimension(:,:,:,:,:), allocatable, save :: grad_nu, grad_u

  integer :: np, ne
  integer :: e, i, j, k

  allocate(nu        , mold=f)
  allocate(u         , mold=f)
  allocate(laplace_u , mold=f)
  allocate(grad_nu   , mold=x)
  allocate(grad_u    , mold=x)

  call problem % GetExactSolution       (x, u)
  call problem % GetExactGradient       (x, grad_u)
  call problem % GetExactLaplacian      (x, laplace_u)
  call problem % GetDiffusivity         (x, nu)
  call problem % GetDiffusivityGradient (x, grad_nu)

  np = size(x,1)
  ne = size(x,4)

  associate(lambda => problem % lambda)

    !$omp do
    do e = 1, ne
      do k = 1, np
      do j = 1, np
      do i = 1, np

        f(i,j,k,e) = lambda * u(i,j,k,e)                        &
                    - ( nu(i,j,k,e) * laplace_u(i,j,k,e)        &
                      + grad_nu(i,j,k,e,1) * grad_u(i,j,k,e,1)  &
                      + grad_nu(i,j,k,e,2) * grad_u(i,j,k,e,2)  &
                      + grad_nu(i,j,k,e,3) * grad_u(i,j,k,e,3) )
      end do
      end do
      end do
    end do

  end associate

  deallocate(nu, u, laplace_u, grad_nu, grad_u)

end subroutine GetSource

!===============================================================================

end module Elliptic_Problem
