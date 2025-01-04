!> summary:  Test problem based on the Taylor-Green vortex
!> author:   Joerg Stiller
!> date:     2022/10/31
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Elliptic_Problem__TGV_Pressure__3D
  use Kind_Parameters, only: RNP
  use Constants,       only: PI
  use Elliptic_Problem__3D

  implicit none
  private

  public :: EllipticProblem_TGV_Pressure_3D

  !-----------------------------------------------------------------------------
  !> Type defining a 2D test problem

  type, extends(EllipticProblem_3D) :: EllipticProblem_TGV_Pressure_3D

    real(RNP) :: t     = 0                     !< time
    real(RNP) :: vt(3) = [1.000, 1.000, 0.000] !< translation velocity
    real(RNP) :: xt(3) = [0.000, 0.125, 0.000] !< initial displacement

  contains

    procedure :: GetExactSolution
    procedure :: GetExactGradient
    procedure :: GetExactLaplacian
    procedure :: GetDiffusivity
    procedure :: GetDiffusivityGradient

  end type EllipticProblem_TGV_Pressure_3D

  ! constructor
  interface EllipticProblem_TGV_Pressure_3D
    module procedure New_Problem
  end interface

contains

  !=============================================================================
  ! Create new object of type EllipticProblem_TGV_Pressure_3D

  function New_Problem(lambda, nu_0, nu_1, d_nu, k_nu) result(this)
    type(EllipticProblem_TGV_Pressure_3D) :: this
    real(RNP), optional :: lambda  !< Helmholtz parameter
    real(RNP), optional :: nu_0    !< diffusivity mean value ν₀
    real(RNP), optional :: nu_1    !< diffusivity fluctuation amplitude ν₁
    real(RNP), optional :: d_nu    !< diffusivity fluctuation phase shift
    integer  , optional :: k_nu    !< diffusivity fluctuation wave number

    if (present(lambda))  this % lambda = lambda
    if (present(nu_0  ))  this % nu_0   = nu_0
    if (present(nu_1  ))  this % nu_1   = nu_1
    if (present(d_nu  ))  this % d_nu   = d_nu
    if (present(k_nu  ))  this % k_nu   = k_nu

  end function New_Problem

  !=============================================================================
  ! GetExactSolution

  !-----------------------------------------------------------------------------
  !> Exact solution

  subroutine GetExactSolution(this, x, u)
    class(EllipticProblem_TGV_Pressure_3D), intent(in) :: this
    real(RNP), intent(in)  :: x(:,:,:,:,:) !< mesh points
    real(RNP), intent(out) :: u(:,:,:,:)   !< solution, u(x)

    integer :: n

    n = size(x(:,:,:,:,1))
    call GetExactSolution_X(this, n, x, u)

  end subroutine GetExactSolution

  !-----------------------------------------------------------------------------
  !> Exact solution -- explicit

  subroutine GetExactSolution_X(this, n, x, u)
    class(EllipticProblem_TGV_Pressure_3D), intent(in) :: this
    integer,   intent(in)  :: n      !< number of points
    real(RNP), intent(in)  :: x(n,3) !< mesh points
    real(RNP), intent(out) :: u(n)   !< solution, u(x)

    real(RNP) :: a, c, phi1, phi2, t
    integer   :: i

    t = this % t
    a = real( -16 * PI**2 * this % nu_0, RNP)
    c = real( exp(a * t) / 4           , RNP)

    do i = 1, n
      phi1 = real(4 * PI * (x(i,1) - this % vt(1)*t - this % xt(1)), RNP)
      phi2 = real(4 * PI * (x(i,2) - this % vt(2)*t - this % xt(2)), RNP)
      u(i) = c * (cos(phi1) + cos(phi2))
    end do

  end subroutine GetExactSolution_X

  !=============================================================================
  ! GetExactGradient

  !-----------------------------------------------------------------------------
  !> Exact solution gradient

  subroutine GetExactGradient(this, x, grad_u)
    class(EllipticProblem_TGV_Pressure_3D), intent(in) :: this
    real(RNP), intent(in)  :: x(:,:,:,:,:)      !< mesh points
    real(RNP), intent(out) :: grad_u(:,:,:,:,:) !< ∇u(x)

    integer :: n

    n = size(x(:,:,:,:,1))
    call GetExactGradient_X(this, n, x, grad_u)

  end subroutine GetExactGradient

  !-----------------------------------------------------------------------------
  !> Exact gradient -- explicit

  subroutine GetExactGradient_X(this, n, x, grad_u)
    class(EllipticProblem_TGV_Pressure_3D), intent(in) :: this
    integer,   intent(in)  :: n           !< number of points
    real(RNP), intent(in)  :: x(n,3)      !< mesh points
    real(RNP), intent(out) :: grad_u(n,3) !< gradient, grad u(x)

    real(RNP) :: a, c, phi1, phi2, t
    integer   :: i

    t = this % t
    a = real( -16 * PI**2 * this % nu_0, RNP)
    c = real( -PI * exp(a * t)            , RNP)

    do i = 1, n

      phi1 = real(4 * PI * (x(i,1) - this % vt(1)*t - this % xt(1)), RNP)
      phi2 = real(4 * PI * (x(i,2) - this % vt(2)*t - this % xt(2)), RNP)

      grad_u(i,1) = c * sin(phi1)
      grad_u(i,2) = c * sin(phi2)
      grad_u(i,3) = 0

    end do

  end subroutine GetExactGradient_X

  !=============================================================================
  ! GetExactLaplacian

  !-----------------------------------------------------------------------------
  !> Exact laplacian

  subroutine GetExactLaplacian(this, x, laplace_u)
    class(EllipticProblem_TGV_Pressure_3D), intent(in) :: this
    real(RNP), intent(in)  :: x(:,:,:,:,:)       !< mesh points
    real(RNP), intent(out) :: laplace_u(:,:,:,:) !< ∇²u(x)

    integer :: n

    n = size(x(:,:,:,:,1))
    call GetExactLaplacian_X(this, n, x, laplace_u)

  end subroutine GetExactLaplacian

  !-----------------------------------------------------------------------------
  !> Exact laplacian -- explicit

  subroutine GetExactLaplacian_X(this, n, x, laplace_u)
    class(EllipticProblem_TGV_Pressure_3D), intent(in) :: this
    integer,   intent(in)  :: n             !< number of points
    real(RNP), intent(in)  :: x(n,3)        !< mesh points
    real(RNP), intent(out) :: laplace_u(n)  !< laplacian, laplace u(x)

    real(RNP) :: a, c, phi1, phi2, t
    integer   :: i

    t = this % t
    a = real( -16 * PI**2 * this % nu_0, RNP)
    c = real( - 4 * PI**2 * exp(a * t) , RNP)

    do i = 1, n

      phi1 = real(4 * PI * (x(i,1) - this % vt(1)*t - this % xt(1)), RNP)
      phi2 = real(4 * PI * (x(i,2) - this % vt(2)*t - this % xt(2)), RNP)

      laplace_u(i) = c * (cos(phi1) + cos(phi2))

    end do

  end subroutine GetExactLaplacian_X

  !=============================================================================
  ! GetDiffusivity

  !-----------------------------------------------------------------------------
  !> Diffusivity

  subroutine GetDiffusivity(this, x, nu)
    class(EllipticProblem_TGV_Pressure_3D), intent(in) :: this
    real(RNP), intent(in)  :: x(:,:,:,:,:)  !< mesh points
    real(RNP), intent(out) :: nu(:,:,:,:)   !< ν(x)

    integer :: n

    n = size(x(:,:,:,:,1))

    call GetDiffusivity_X( this % nu_0, &
                           this % nu_1, &
                           this % k_nu, &
                           this % d_nu, &
                           n, x, nu     )

  end subroutine GetDiffusivity

  !-----------------------------------------------------------------------------
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

  !=============================================================================
  ! GetDiffusivityGradient

  !-----------------------------------------------------------------------------
  !> Diffusivity gradient

  subroutine GetDiffusivityGradient(this, x, grad_nu)
    class(EllipticProblem_TGV_Pressure_3D), intent(in) :: this
    real(RNP), intent(in)  :: x(:,:,:,:,:)       !< mesh points
    real(RNP), intent(out) :: grad_nu(:,:,:,:,:) !< ∇ν(x)

    integer :: n

    n = size(x(:,:,:,:,1))

    call GetDiffusivityGradient_X( this % nu_1,  &
                                   this % k_nu,  &
                                   this % d_nu,  &
                                   n, x, grad_nu )

  end subroutine GetDiffusivityGradient

  !-----------------------------------------------------------------------------
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

      grad_nu(i,1) = k * nu_1 * c1 * s2
      grad_nu(i,2) = k * nu_1 * s1 * c2
      grad_nu(i,3) = 0

    end do

  end subroutine GetDiffusivityGradient_X

  !=============================================================================

end module Elliptic_Problem__TGV_Pressure__3D
