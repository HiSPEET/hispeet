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

!> summary:  Defines a simple periodic 3D test problem
!> author:   Joerg Stiller
!> date:     2019/01/27
!===============================================================================

module Elliptic_Problem__Simple_3D
  use Kind_Parameters, only: RNP
  use Elliptic_Problem

  implicit none
  private

  public :: EllipticProblem_Simple3D

  !-----------------------------------------------------------------------------
  !> Type defining a 3D test problem

  type, extends(EllipticProblem) :: EllipticProblem_Simple3D
  contains

    procedure :: GetExactSolution
    procedure :: GetExactGradient
    procedure :: GetExactLaplacian
    procedure :: GetDiffusivity
    procedure :: GetDiffusivityGradient

  end type EllipticProblem_Simple3D

contains

!===============================================================================
! GetExactSolution

!---------------------------------------------------------------------------
!> Exact solution

subroutine GetExactSolution(problem, x, u)
  class(EllipticProblem_Simple3D), intent(in)  :: problem
  real(RNP),              intent(in)  :: x(:,:,:,:,:) !< mesh points
  real(RNP),              intent(out) :: u(:,:,:,:)   !< solution, u(x)

  integer :: n

  n = size(x(:,:,:,:,1))
  call GetExactSolution_X(problem % k_u, n, x, u)

end subroutine GetExactSolution

!-------------------------------------------------------------------------------
!> Exact solution -- explicit

subroutine GetExactSolution_X(k, n, x, u)
  integer,   intent(in)  :: k      !< wave number
  integer,   intent(in)  :: n      !< number of points
  real(RNP), intent(in)  :: x(n,3) !< mesh points
  real(RNP), intent(out) :: u(n)   !< solution, u(x)

  real(RNP) :: x1, x2, x3
  integer   :: i

  do i = 1, n

    x1 = x(i,1)
    x2 = x(i,2)
    x3 = x(i,3)

    u(i) = sin(k * x1) * sin(k * x2) * sin(k * x3)

  end do

end subroutine GetExactSolution_X

!===============================================================================
! GetExactGradient

!---------------------------------------------------------------------------
!> Exact solution gradient

subroutine GetExactGradient(problem, x, grad_u)
  class(EllipticProblem_Simple3D), intent(in)  :: problem
  real(RNP), intent(in)  :: x(:,:,:,:,:)      !< mesh points
  real(RNP), intent(out) :: grad_u(:,:,:,:,:) !< ∇u(x)

  integer :: n

  n = size(x(:,:,:,:,1))
  call GetExactGradient_X(problem % k_u, n, x, grad_u)

end subroutine GetExactGradient

!-------------------------------------------------------------------------------
!> Exact gradient -- explicit

subroutine GetExactGradient_X(k, n, x, grad_u)
  integer,   intent(in)  :: k           !< wave number
  integer,   intent(in)  :: n           !< number of points
  real(RNP), intent(in)  :: x(n,3)      !< mesh points
  real(RNP), intent(out) :: grad_u(n,3) !< gradient, grad u(x)

  real(RNP) :: x1, x2, x3
  integer   :: i

  do i = 1, n

    x1 = x(i,1)
    x2 = x(i,2)
    x3 = x(i,3)

    grad_u(i,1) = k * cos(k * x1) * sin(k * x2) * sin(k * x3)
    grad_u(i,2) = k * sin(k * x1) * cos(k * x2) * sin(k * x3)
    grad_u(i,3) = k * sin(k * x1) * sin(k * x2) * cos(k * x3)

  end do

end subroutine GetExactGradient_X

!===============================================================================
! GetExactLaplacian

!---------------------------------------------------------------------------
!> Exact laplacian

subroutine GetExactLaplacian(problem, x, laplace_u)
  class(EllipticProblem_Simple3D), intent(in)  :: problem
  real(RNP), intent(in)  :: x(:,:,:,:,:)       !< mesh points
  real(RNP), intent(out) :: laplace_u(:,:,:,:) !< ∇²u(x)

  integer :: n

  n = size(x(:,:,:,:,1))
  call GetExactLaplacian_X(problem % k_u, n, x, laplace_u)

end subroutine GetExactLaplacian

!-------------------------------------------------------------------------------
!> Exact laplacian -- explicit

subroutine GetExactLaplacian_X(k, n, x, laplace_u)
  integer,   intent(in)  :: k             !< wave number
  integer,   intent(in)  :: n             !< number of points
  real(RNP), intent(in)  :: x(n,3)        !< mesh points
  real(RNP), intent(out) :: laplace_u(n)  !< laplacian, laplace u(x)

  real(RNP) :: x1, x2, x3
  integer   :: i

  do i = 1, n

    x1 = x(i,1)
    x2 = x(i,2)
    x3 = x(i,3)

    laplace_u(i) = -3*k*k * sin(k * x1) * sin(k * x2) * sin(k * x3)

  end do

end subroutine GetExactLaplacian_X

!===============================================================================
! GetDiffusivity

!---------------------------------------------------------------------------
!> Diffusivity

subroutine GetDiffusivity(problem, x, nu)
  class(EllipticProblem_Simple3D), intent(in)  :: problem
  real(RNP), intent(in)  :: x(:,:,:,:,:)  !< mesh points
  real(RNP), intent(out) :: nu(:,:,:,:)   !< ν(x)

  integer :: n

  n = size(x(:,:,:,:,1))

  call GetDiffusivity_X( problem % nu_0, &
                         problem % nu_1, &
                         problem % k_nu, &
                         problem % d_nu, &
                         n, x, nu        )

end subroutine GetDiffusivity

!-------------------------------------------------------------------------------
!> Diffusivity -- explicit

subroutine GetDiffusivity_X(nu_0, nu_1, k, d, n, x, nu)
  real(RNP), intent(in)  :: nu_0    !< constant part ν₀
  real(RNP), intent(in)  :: nu_1    !< fluctuation amplitude ν₁
  integer,   intent(in)  :: k       !< fluctuation wave number
  real(RNP), intent(in)  :: d       !< fluctuation phase shift
  integer,   intent(in)  :: n       !< number of points
  real(RNP), intent(in)  :: x(n,3)  !< mesh points
  real(RNP), intent(out) :: nu(n)   !< diffusivity ν(x)

  real(RNP) :: x1, x2, x3
  integer   :: i

  do i = 1, n

    x1 = x(i,1)
    x2 = x(i,2)
    x3 = x(i,3)

    nu(i) = nu_0                      &
          + nu_1 * sin(k * (x1 - d))  &
                 * sin(k * (x2 - d))  &
                 * sin(k * (x3 - d))
  end do

end subroutine GetDiffusivity_X

!===============================================================================
! GetDiffusivityGradient

!---------------------------------------------------------------------------
!> Diffusivity gradient

subroutine GetDiffusivityGradient(problem, x, grad_nu)
  class(EllipticProblem_Simple3D), intent(in)  :: problem
  real(RNP), intent(in)  :: x(:,:,:,:,:)       !< mesh points
  real(RNP), intent(out) :: grad_nu(:,:,:,:,:) !< ∇ν(x)

  integer :: n

  n = size(x(:,:,:,:,1))

  call GetDiffusivityGradient_X( problem % nu_1, &
                                 problem % k_nu, &
                                 problem % d_nu, &
                                 n, x, grad_nu   )

end subroutine GetDiffusivityGradient

!-------------------------------------------------------------------------------
!> Diffusivity gradient -- explicit

subroutine GetDiffusivityGradient_X(nu_1, k, d, n, x, grad_nu)
  real(RNP), intent(in)  :: nu_1          !< fluctuation amplitude ν₁
  integer,   intent(in)  :: k             !< fluctuation wave number
  real(RNP), intent(in)  :: d             !< fluctuation phase shift
  integer,   intent(in)  :: n             !< number of points
  real(RNP), intent(in)  :: x(n,3)        !< mesh points
  real(RNP), intent(out) :: grad_nu(n,3)  !< diffusivity gradient ∇ν(x)

  real(RNP) :: s1, s2, s3
  real(RNP) :: c1, c2, c3
  integer   :: i

  do i = 1, n

    s1 = sin(k * (x(i,1) - d))
    c1 = cos(k * (x(i,1) - d))

    s2 = sin(k * (x(i,2) - d))
    c2 = cos(k * (x(i,2) - d))

    s3 = sin(k * (x(i,3) - d))
    c3 = cos(k * (x(i,3) - d))

    grad_nu(i,1) = k * nu_1 * c1 * s2 * s3
    grad_nu(i,2) = k * nu_1 * s1 * c2 * s3
    grad_nu(i,3) = k * nu_1 * s1 * s2 * c3

  end do

end subroutine GetDiffusivityGradient_X

!===============================================================================

end module Elliptic_Problem__Simple_3D
