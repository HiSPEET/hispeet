!> summary:  3D elliptic test case featuring a spherical front
!> author:   Joerg Stiller
!> date:     2025/01/03
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> Implements the exact solution for the test case proposed in [1]
!>
!> 1. Cerveny J, Dobrev V, Kolev T: Nonconforming mesh refinement for high-order
!>    finite elements. SIAM J SCI COMPUT 41, C367–C392, 2019
!===============================================================================

module Elliptic_Problem__Sphere__3D
  use Kind_Parameters, only: RNP
  use Constants, only: ZERO, ONE
  use Array_Assignments
  use Elliptic_Problem__3D

  implicit none
  private

  public :: EllipticProblem_Sphere_3D

  !-----------------------------------------------------------------------------
  !> Type defining a simple elliptic test problem

  type, extends(EllipticProblem_3D) :: EllipticProblem_Sphere_3D

    real(RNP) :: x_c(3) = -0.05 !< sphere center
    real(RNP) :: r_0    =  0.7  !< sphere radius
    real(RNP) :: alpha  =  200  !< radial scaling factor

  contains

    procedure :: GetExactSolution
    procedure :: GetExactGradient
    procedure :: GetExactLaplacian
    procedure :: GetDiffusivity
    procedure :: GetDiffusivityGradient

  end type EllipticProblem_Sphere_3D

  ! constructor
  interface EllipticProblem_Sphere_3D
    module procedure New_Problem
  end interface

  !-----------------------------------------------------------------------------
  !> tolerance

  real(RNP), parameter :: eps = 1E-15

contains

  !=============================================================================
  ! Create new object of type EllipticProblem_Sphere_3D

  function New_Problem(lambda, x_c, r_0, alpha) result(this)
    type(EllipticProblem_Sphere_3D) :: this
    real(RNP), optional :: lambda !< Helmholtz parameter
    real(RNP), optional :: x_c(3) !< sphere center
    real(RNP), optional :: r_0    !< sphere radius
    real(RNP), optional :: alpha  !< radial scaling factor

    if (present(lambda))  this % lambda = lambda
    if (present(x_c   ))  this % x_c    = x_c
    if (present(r_0   ))  this % r_0    = r_0
    if (present(alpha ))  this % alpha  = alpha

  end function New_Problem

  !=============================================================================
  ! GetExactSolution

  !-----------------------------------------------------------------------------
  !> Exact solution

  subroutine GetExactSolution(this, x, u)
    class(EllipticProblem_Sphere_3D), intent(in) :: this
    real(RNP), intent(in)  :: x(:,:,:,:,:) !< mesh points
    real(RNP), intent(out) :: u(:,:,:,:)   !< solution, u(x)

    integer :: n

    n = size(x(:,:,:,:,1))
    call GetExactSolution_X(this%x_c, this%r_0, this%alpha, n, x, u)

  end subroutine GetExactSolution


  !-----------------------------------------------------------------------------
  !> Exact solution -- explicit

  subroutine GetExactSolution_X(x_c, r_0, alpha, n, x, u)
    real(RNP), intent(in)  :: x_c(3) !< sphere center
    real(RNP), intent(in)  :: r_0    !< sphere radius
    real(RNP), intent(in)  :: alpha  !< radial scaling factor
    integer,   intent(in)  :: n      !< number of points
    real(RNP), intent(in)  :: x(n,3) !< mesh points
    real(RNP), intent(out) :: u(n)   !< solution, u(x)

    real(RNP) :: r
    integer   :: i

    do i = 1, n

      r = sqrt( (x(i,1) - x_c(1))**2 &
              + (x(i,2) - x_c(2))**2 &
              + (x(i,3) - x_c(3))**2 )

      u(i) = atan( alpha * (r - r_0) )

    end do

  end subroutine GetExactSolution_X

  !=============================================================================
  ! GetExactGradient

  !-----------------------------------------------------------------------------
  !> Exact solution gradient

  subroutine GetExactGradient(this, x, grad_u)
    class(EllipticProblem_Sphere_3D), intent(in) :: this
    real(RNP), intent(in)  :: x(:,:,:,:,:)      !< mesh points
    real(RNP), intent(out) :: grad_u(:,:,:,:,:) !< ∇u(x)

    integer :: n

    n = size(x(:,:,:,:,1))
    call GetExactGradient_X(this%x_c, this%r_0, this%alpha, n, x, grad_u)

  end subroutine GetExactGradient

  !-----------------------------------------------------------------------------
  !> Exact gradient -- explicit

  subroutine GetExactGradient_X(x_c, r_0, alpha, n, x, grad_u)
    real(RNP), intent(in)  :: x_c(3)      !< sphere center
    real(RNP), intent(in)  :: r_0         !< sphere radius
    real(RNP), intent(in)  :: alpha       !< radial scaling factor
    integer,   intent(in)  :: n           !< number of points
    real(RNP), intent(in)  :: x(n,3)      !< mesh points
    real(RNP), intent(out) :: grad_u(n,3) !< gradient, grad u(x)

    real(RNP) :: d, r, u, u_r
    integer   :: i

    do i = 1, n

      r = sqrt( (x(i,1) - x_c(1))**2 &
              + (x(i,2) - x_c(2))**2 &
              + (x(i,3) - x_c(3))**2 )

      d   = r - r_0
      u   = atan( alpha * d )              ! solution
      u_r = alpha / ((alpha * d)**2 + 1)   ! 1st radial derivative

      grad_u(i,1) = u_r * (x(i,1) - x_c(1)) / max(r, eps)
      grad_u(i,2) = u_r * (x(i,2) - x_c(2)) / max(r, eps)
      grad_u(i,3) = u_r * (x(i,3) - x_c(3)) / max(r, eps)

    end do

  end subroutine GetExactGradient_X

  !=============================================================================
  ! GetExactLaplacian

  !-----------------------------------------------------------------------------
  !> Exact laplacian

  subroutine GetExactLaplacian(this, x, laplace_u)
    class(EllipticProblem_Sphere_3D), intent(in) :: this
    real(RNP), intent(in)  :: x(:,:,:,:,:)       !< mesh points
    real(RNP), intent(out) :: laplace_u(:,:,:,:) !< ∇²u(x)

    integer :: n

    n = size(x(:,:,:,:,1))
    call GetExactLaplacian_X(this%x_c, this%r_0, this%alpha, n, x, laplace_u)

  end subroutine GetExactLaplacian

  !-----------------------------------------------------------------------------
  !> Exact laplacian -- explicit

  subroutine GetExactLaplacian_X(x_c, r_0, alpha, n, x, laplace_u)
    real(RNP), intent(in)  :: x_c(3)       !< sphere center
    real(RNP), intent(in)  :: r_0          !< sphere radius
    real(RNP), intent(in)  :: alpha        !< radial scaling factor
    integer,   intent(in)  :: n            !< number of points
    real(RNP), intent(in)  :: x(n,3)       !< mesh points
    real(RNP), intent(out) :: laplace_u(n) !< laplacian, laplace u(x)

    real(RNP) :: d, r, u_r, u_rr
    integer   :: i

    do i = 1, n

      r = sqrt( (x(i,1) - x_c(1))**2 &
              + (x(i,2) - x_c(2))**2 &
              + (x(i,3) - x_c(3))**2 )

      d    = r - r_0
      u_r  = alpha / ((alpha * d)**2 + 1) ! 1st radial derivative
      u_rr = -2 * alpha * d * u_r**2      ! 2nd radial derivative

      laplace_u(i) = u_rr + 2 * u_r / max(r, eps)

    end do

  end subroutine GetExactLaplacian_X

  !=============================================================================
  ! GetDiffusivity

  !-----------------------------------------------------------------------------
  !> Diffusivity

  subroutine GetDiffusivity(this, x, nu)
    class(EllipticProblem_Sphere_3D), intent(in) :: this
    real(RNP), intent(in)  :: x(:,:,:,:,:)  !< mesh points
    real(RNP), intent(out) :: nu(:,:,:,:)   !< ν(x)

    call SetArray(nu, ONE)

    if (size(x) == 0) return

  end subroutine GetDiffusivity

  !=============================================================================
  ! GetDiffusivityGradient

  !-----------------------------------------------------------------------------
  !> Diffusivity gradient

  subroutine GetDiffusivityGradient(this, x, grad_nu)
    class(EllipticProblem_Sphere_3D), intent(in) :: this
    real(RNP), intent(in)  :: x(:,:,:,:,:)       !< mesh points
    real(RNP), intent(out) :: grad_nu(:,:,:,:,:) !< ∇ν(x)

    call SetArray(grad_nu, ZERO)

    if (this%nu_0 == 0 .or. size(x) == 0) return

  end subroutine GetDiffusivityGradient

  !=============================================================================

end module Elliptic_Problem__Sphere__3D
