!> summary:  Defines a simple periodic 2D test problem
!> author:   Joerg Stiller
!> date:     2019/01/27
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Elliptic_Problem__Simple_2D
  use Kind_Parameters, only: RNP
  use Elliptic_Problem

  implicit none
  private

  public :: EllipticProblem_Simple2D

  !-----------------------------------------------------------------------------
  !> Type defining a 2D test problem

  type, extends(EllipticProblem) :: EllipticProblem_Simple2D
  contains

    procedure :: GetExactSolution
    procedure :: GetExactGradient
    procedure :: GetExactLaplacian
    procedure :: GetDiffusivity
    procedure :: GetDiffusivityGradient

  end type EllipticProblem_Simple2D

contains

!===============================================================================
! GetExactSolution

!---------------------------------------------------------------------------
!> Exact solution

subroutine GetExactSolution(problem, x, u)
  class(EllipticProblem_Simple2D), intent(in)  :: problem
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

  real(RNP) :: x1, x2
  integer   :: i

  do i = 1, n

    x1 = x(i,1)
    x2 = x(i,2)

    u(i) = sin(k * x1) * sin(k * x2)

  end do

end subroutine GetExactSolution_X

!===============================================================================
! GetExactGradient

!---------------------------------------------------------------------------
!> Exact solution gradient

subroutine GetExactGradient(problem, x, grad_u)
  class(EllipticProblem_Simple2D), intent(in)  :: problem
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

  real(RNP) :: x1, x2
  integer   :: i

  do i = 1, n

    x1 = x(i,1)
    x2 = x(i,2)

    grad_u(i,1) = -k * cos(k * x1) * sin(k * x2)
    grad_u(i,2) = -k * sin(k * x1) * cos(k * x2)
    grad_u(i,3) =  0

  end do

end subroutine GetExactGradient_X

!===============================================================================
! GetExactLaplacian

!---------------------------------------------------------------------------
!> Exact laplacian

subroutine GetExactLaplacian(problem, x, laplace_u)
  class(EllipticProblem_Simple2D), intent(in)  :: problem
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

  real(RNP) :: x1, x2
  integer   :: i

  do i = 1, n

    x1 = x(i,1)
    x2 = x(i,2)

    laplace_u(i) = -2*k*k * sin(k * x1) * sin(k * x2)

  end do

end subroutine GetExactLaplacian_X

!===============================================================================
! GetDiffusivity

!---------------------------------------------------------------------------
!> Diffusivity

subroutine GetDiffusivity(problem, x, nu)
  class(EllipticProblem_Simple2D), intent(in)  :: problem
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

  real(RNP) :: x1, x2
  integer   :: i

  do i = 1, n

    x1 = x(i,1)
    x2 = x(i,2)

    nu(i) = nu_0  +  nu_1 * sin(k * (x1 - d)) * sin(k * (x2 - d))

  end do

end subroutine GetDiffusivity_X

!===============================================================================
! GetDiffusivityGradient

!---------------------------------------------------------------------------
!> Diffusivity gradient

subroutine GetDiffusivityGradient(problem, x, grad_nu)
  class(EllipticProblem_Simple2D), intent(in)  :: problem
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

  real(RNP) :: s1, s2
  real(RNP) :: c1, c2
  integer   :: i

  do i = 1, n

    s1 = sin(k * (x(i,1) - d))
    c1 = cos(k * (x(i,1) - d))

    s2 = sin(k * (x(i,2) - d))
    c2 = cos(k * (x(i,2) - d))

    grad_nu(i,1) = -k * nu_1 * c1 * s2
    grad_nu(i,2) = -k * nu_1 * s1 * c2
    grad_nu(i,3) =  0

  end do

end subroutine GetDiffusivityGradient_X

!===============================================================================

end module Elliptic_Problem__Simple_2D
