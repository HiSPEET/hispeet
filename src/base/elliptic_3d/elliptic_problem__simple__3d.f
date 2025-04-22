!> summary:  Defines a simple periodic test problem with 1D/2D/3D solution
!> author:   Joerg Stiller
!> date:     2019/01/27, 2024/09/17
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> Test problems with simple analytical solution
!>
!>       λ - ∇⋅ν∇u = f
!>
!> with manufactured solution
!>
!>   - 1D:  u = sin(k_u x)
!>   - 2D:  u = sin(k_u x) sin(k_u y)
!>   - 3D:  u = sin(k_u x) sin(k_u y) sin(k_u z)
!>
!> for constant λ and possibly variable diffusivity
!>
!>   - 1D:  ν = ν₀ + ν₁ sin(k_nu (x - d))
!>   - 2D:  ν = ν₀ + ν₁ sin(k_nu (x - d)) sin(k_nu (y - d))
!>   - 3D:  ν = ν₀ + ν₁ sin(k_nu (x - d)) sin(k_nu (y - d)) sin(k_nu (z - d))
!>
!===============================================================================

module Elliptic_Problem__Simple__3D
  use Kind_Parameters, only: RNP
  use Elliptic_Problem__3D

  implicit none
  private

  public :: EllipticProblem_Simple_3D

  !-----------------------------------------------------------------------------
  !> Type defining a simple elliptic test problem

  type, extends(EllipticProblem_3D) :: EllipticProblem_Simple_3D

    integer :: k_u = 1  !< solution wave number
    integer :: dim = 3  !< solution dimensionality {1,2,3}

  contains

    procedure :: GetExactSolution
    procedure :: GetExactGradient
    procedure :: GetExactLaplacian
    procedure :: GetDiffusivity
    procedure :: GetDiffusivityGradient

  end type EllipticProblem_Simple_3D

  ! constructor
  interface EllipticProblem_Simple_3D
    module procedure New_Problem
  end interface

contains

  !=============================================================================
  ! Create new object of type EllipticProblem_Simple_3D

  function New_Problem(lambda, nu_0, nu_1, d_nu, k_nu, k_u, dim) result(this)
    type(EllipticProblem_Simple_3D) :: this
    real(RNP), optional :: lambda !< Helmholtz parameter
    real(RNP), optional :: nu_0   !< diffusivity mean value ν₀
    real(RNP), optional :: nu_1   !< diffusivity fluctuation amplitude ν₁
    real(RNP), optional :: d_nu   !< diffusivity fluctuation phase shift
    integer  , optional :: k_nu   !< diffusivity fluctuation wave number
    integer  , optional :: k_u    !< solution wave number
    integer  , optional :: dim    !< solution dimensionality

    if (present(lambda))  this % lambda = lambda
    if (present(nu_0  ))  this % nu_0   = nu_0
    if (present(nu_1  ))  this % nu_1   = nu_1
    if (present(d_nu  ))  this % d_nu   = d_nu
    if (present(k_nu  ))  this % k_nu   = k_nu
    if (present(k_u   ))  this % k_u    = k_u
    if (present(dim   ))  this % dim    = dim

  end function New_Problem

  !=============================================================================
  ! GetExactSolution

  !-----------------------------------------------------------------------------
  !> Exact solution

  subroutine GetExactSolution(this, x, u)
    class(EllipticProblem_Simple_3D), intent(in) :: this
    real(RNP), intent(in)  :: x(:,:,:,:,:) !< mesh points
    real(RNP), intent(out) :: u(:,:,:,:)   !< solution, u(x)

    integer :: n

    n = size(x(:,:,:,:,1))

    select case(this % dim)
    case(1)
      call GetExactSolution1D(this % k_u, n, x, u)
    case(2)
      call GetExactSolution2D(this % k_u, n, x, u)
    case default
      call GetExactSolution3D(this % k_u, n, x, u)
    end select

  end subroutine GetExactSolution

  !-----------------------------------------------------------------------------
  !> Exact 1D solution

  subroutine GetExactSolution1D(k, n, x, u)
    integer,   intent(in)  :: k      !< wave number
    integer,   intent(in)  :: n      !< number of points
    real(RNP), intent(in)  :: x(n,3) !< mesh points
    real(RNP), intent(out) :: u(n)   !< solution, u(x)

    integer   :: i

    do i = 1, n

      u(i) = sin(k * x(i,1))

    end do

  end subroutine GetExactSolution1D

  !-----------------------------------------------------------------------------
  !> Exact 2D solution

  subroutine GetExactSolution2D(k, n, x, u)
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

  end subroutine GetExactSolution2D

  !-----------------------------------------------------------------------------
  !> Exact 3D solution

  subroutine GetExactSolution3D(k, n, x, u)
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

  end subroutine GetExactSolution3D

  !=============================================================================
  ! GetExactGradient

  !-----------------------------------------------------------------------------
  !> Exact solution gradient

  subroutine GetExactGradient(this, x, grad_u)
    class(EllipticProblem_Simple_3D), intent(in) :: this
    real(RNP), intent(in)  :: x(:,:,:,:,:)      !< mesh points
    real(RNP), intent(out) :: grad_u(:,:,:,:,:) !< ∇u(x)

    integer :: n

    n = size(x(:,:,:,:,1))

    select case(this % dim)
    case(1)
      call GetExactGradient1D(this % k_u, n, x, grad_u)
    case(2)
      call GetExactGradient2D(this % k_u, n, x, grad_u)
    case default
      call GetExactGradient3D(this % k_u, n, x, grad_u)
    end select

  end subroutine GetExactGradient

  !-----------------------------------------------------------------------------
  !> Exact gradient of 1D solution

  subroutine GetExactGradient1D(k, n, x, grad_u)
    integer,   intent(in)  :: k           !< wave number
    integer,   intent(in)  :: n           !< number of points
    real(RNP), intent(in)  :: x(n,3)      !< mesh points
    real(RNP), intent(out) :: grad_u(n,3) !< gradient, grad u(x)

    integer   :: i

    do i = 1, n

      grad_u(i,1) =  k * cos(k * x(i,1))
      grad_u(i,2) =  0
      grad_u(i,3) =  0

    end do

  end subroutine GetExactGradient1D

  !-----------------------------------------------------------------------------
  !> Exact gradient of 2D solution

  subroutine GetExactGradient2D(k, n, x, grad_u)
    integer,   intent(in)  :: k           !< wave number
    integer,   intent(in)  :: n           !< number of points
    real(RNP), intent(in)  :: x(n,3)      !< mesh points
    real(RNP), intent(out) :: grad_u(n,3) !< gradient, grad u(x)

    real(RNP) :: x1, x2
    integer   :: i

    do i = 1, n

      x1 = x(i,1)
      x2 = x(i,2)

      grad_u(i,1) = k * cos(k * x1) * sin(k * x2)
      grad_u(i,2) = k * sin(k * x1) * cos(k * x2)
      grad_u(i,3) = 0

    end do

  end subroutine GetExactGradient2D

  !-----------------------------------------------------------------------------
  !> Exact gradient of 3D solution

  subroutine GetExactGradient3D(k, n, x, grad_u)
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

  end subroutine GetExactGradient3D

  !=============================================================================
  ! GetExactLaplacian

  !-----------------------------------------------------------------------------
  !> Exact laplacian

  subroutine GetExactLaplacian(this, x, laplace_u)
    class(EllipticProblem_Simple_3D), intent(in) :: this
    real(RNP), intent(in)  :: x(:,:,:,:,:)       !< mesh points
    real(RNP), intent(out) :: laplace_u(:,:,:,:) !< ∇²u(x)

    integer :: n

    n = size(x(:,:,:,:,1))

    select case(this % dim)
    case(1)
      call GetExactLaplacian1D(this % k_u, n, x, laplace_u)
    case(2)
      call GetExactLaplacian2D(this % k_u, n, x, laplace_u)
    case default
      call GetExactLaplacian3D(this % k_u, n, x, laplace_u)
    end select

  end subroutine GetExactLaplacian

  !-----------------------------------------------------------------------------
  !> Exact laplacian of 1D solution

  subroutine GetExactLaplacian1D(k, n, x, laplace_u)
    integer,   intent(in)  :: k             !< wave number
    integer,   intent(in)  :: n             !< number of points
    real(RNP), intent(in)  :: x(n,3)        !< mesh points
    real(RNP), intent(out) :: laplace_u(n)  !< laplacian, laplace u(x)

    integer   :: i

    do i = 1, n

      laplace_u(i) = -k*k * sin(k * x(i,1))

    end do

  end subroutine GetExactLaplacian1D

  !-----------------------------------------------------------------------------
  !> Exact laplacian of 2D solution

  subroutine GetExactLaplacian2D(k, n, x, laplace_u)
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

  end subroutine GetExactLaplacian2D

  !-----------------------------------------------------------------------------
  !> Exact laplacian of 3D solution

  subroutine GetExactLaplacian3D(k, n, x, laplace_u)
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

  end subroutine GetExactLaplacian3D

  !=============================================================================
  ! GetDiffusivity

  !-----------------------------------------------------------------------------
  !> Diffusivity

  subroutine GetDiffusivity(this, x, nu)
    class(EllipticProblem_Simple_3D), intent(in) :: this
    real(RNP), intent(in)  :: x(:,:,:,:,:)  !< mesh points
    real(RNP), intent(out) :: nu(:,:,:,:)   !< ν(x)

    integer :: n

    n = size(x(:,:,:,:,1))

    select case(this % dim)
    case(1)
      call GetDiffusivity1D( this % nu_0, &
                             this % nu_1, &
                             this % k_nu, &
                             this % d_nu, &
                             n, x, nu     )
    case(2)
      call GetDiffusivity2D( this % nu_0, &
                             this % nu_1, &
                             this % k_nu, &
                             this % d_nu, &
                             n, x, nu     )
    case default
      call GetDiffusivity3D( this % nu_0, &
                             this % nu_1, &
                             this % k_nu, &
                             this % d_nu, &
                             n, x, nu     )
    end select

  end subroutine GetDiffusivity

  !-----------------------------------------------------------------------------
  !> 1D diffusivity

  subroutine GetDiffusivity1D(nu_0, nu_1, k, d, n, x, nu)
    real(RNP), intent(in)  :: nu_0    !< constant part ν₀
    real(RNP), intent(in)  :: nu_1    !< fluctuation amplitude ν₁
    integer,   intent(in)  :: k       !< fluctuation wave number
    real(RNP), intent(in)  :: d       !< fluctuation phase shift
    integer,   intent(in)  :: n       !< number of points
    real(RNP), intent(in)  :: x(n,3)  !< mesh points
    real(RNP), intent(out) :: nu(n)   !< diffusivity ν(x)

    integer :: i

    do i = 1, n

      nu(i) = nu_0  +  nu_1 * sin(k * (x(i,1) - d))

    end do

  end subroutine GetDiffusivity1D

  !-----------------------------------------------------------------------------
  !> 2D diffusivity

  subroutine GetDiffusivity2D(nu_0, nu_1, k, d, n, x, nu)
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

  end subroutine GetDiffusivity2D

  !-----------------------------------------------------------------------------
  !> 3D diffusivity

  subroutine GetDiffusivity3D(nu_0, nu_1, k, d, n, x, nu)
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

  end subroutine GetDiffusivity3D

  !=============================================================================
  ! GetDiffusivityGradient

  !-----------------------------------------------------------------------------
  !> Diffusivity gradient

  subroutine GetDiffusivityGradient(this, x, grad_nu)
    class(EllipticProblem_Simple_3D), intent(in) :: this
    real(RNP), intent(in)  :: x(:,:,:,:,:)       !< mesh points
    real(RNP), intent(out) :: grad_nu(:,:,:,:,:) !< ∇ν(x)

    integer :: n

    n = size(x(:,:,:,:,1))

    select case(this % dim)
    case(1)
      call GetDiffusivityGradient1D( this % nu_1,  &
                                     this % k_nu,  &
                                     this % d_nu,  &
                                     n, x, grad_nu )
    case(2)
      call GetDiffusivityGradient2D( this % nu_1,  &
                                     this % k_nu,  &
                                     this % d_nu,  &
                                     n, x, grad_nu )
    case default
      call GetDiffusivityGradient3D( this % nu_1,  &
                                     this % k_nu,  &
                                     this % d_nu,  &
                                     n, x, grad_nu )
    end select

  end subroutine GetDiffusivityGradient

  !-----------------------------------------------------------------------------
  !> 1D diffusivity gradient

  subroutine GetDiffusivityGradient1D(nu_1, k, d, n, x, grad_nu)
    real(RNP), intent(in)  :: nu_1          !< fluctuation amplitude ν₁
    integer,   intent(in)  :: k             !< fluctuation wave number
    real(RNP), intent(in)  :: d             !< fluctuation phase shift
    integer,   intent(in)  :: n             !< number of points
    real(RNP), intent(in)  :: x(n,3)        !< mesh points
    real(RNP), intent(out) :: grad_nu(n,3)  !< diffusivity gradient ∇ν(x)

    integer   :: i

    do i = 1, n

      grad_nu(i,1) =  k * nu_1 * cos(k * (x(i,1) - d))
      grad_nu(i,2) =  0
      grad_nu(i,3) =  0

    end do

  end subroutine GetDiffusivityGradient1D

  !-----------------------------------------------------------------------------
  !> 2D diffusivity gradient

  subroutine GetDiffusivityGradient2D(nu_1, k, d, n, x, grad_nu)
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

      grad_nu(i,1) = k * nu_1 * c1 * s2
      grad_nu(i,2) = k * nu_1 * s1 * c2
      grad_nu(i,3) = 0

    end do

  end subroutine GetDiffusivityGradient2D

  !-----------------------------------------------------------------------------
  !> 3D diffusivity gradient

  subroutine GetDiffusivityGradient3D(nu_1, k, d, n, x, grad_nu)
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

  end subroutine GetDiffusivityGradient3D

  !=============================================================================

end module Elliptic_Problem__Simple__3D
