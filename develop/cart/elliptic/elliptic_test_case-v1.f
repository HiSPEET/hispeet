!> summary:  Test case for elliptic equations, 2π-periodic
!> author:   Joerg Stiller
!> date:     2013/06/16, extended to variable diffusivity 2019/01/08
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Test case for elliptic equations, 2π-periodic
!===============================================================================

module Elliptic_Test_Case
  use Kind_Parameters, only: RNP
  use Constants,       only: PI
  implicit none
  private

  public :: GetExactSolution
  public :: GetExactGradient
  public :: GetExactLaplacian
  public :: GetDiffusivity
  public :: GetDiffusivityGradient
  public :: GetSource

contains

!-------------------------------------------------------------------------------
!> Exact solution

subroutine GetExactSolution(k, n, x, u)
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

    u(i) = cos(k * (x1 - 3*x2 + 2*x3))   &
         * sin(k * (1 + x1))             &
         * sin(k * (1 - x2))             &
         * sin(k * (2*x1 + x2))          &
         * sin(k * (3*x1 - 2*x2 + 2*x3))

  end do

end subroutine GetExactSolution

!-------------------------------------------------------------------------------
!> Exact gradient

subroutine GetExactGradient(k, n, x, grad_u)
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

    grad_u(i,1) = ( k * (2*cos(k*(1 + x1))                                &
                  - 5 * cos(k*(1 + 5*x1 + 2*x2))                          &
                  + 3 * cos(k - 3*k*x1 - 2*k*x2)                          &
                  + 5 * cos(k*(1 - 5*x1 + 4*x2 - 4*x3))                   &
                  - cos(k*(-1 + x1 - 6*x2 + 4*x3))                        &
                  + 3 * cos(k*(1 + 3*x1 - 6*x2 + 4*x3))                   &
                  - 7 * cos(k*(1 + 7*x1 - 4*x2 + 4*x3))) * sin(k - k*x2)  &
                  ) / 8

    grad_u(i,2) = k * sin(k*(1 + x1))                            &
                    * ( ( cos(k*(x1 - 3*x2 + 2*x3))              &
                        * ( -2 * cos(k*(1 + x1 - 4*x2 + 2*x3))   &
                          + cos(k*(-1 + x1 - 2*x2 + 2*x3))       &
                          + cos(k*(1 + 5*x1 - 2*x2 + 2*x3))      &
                          )                                      &
                        ) / 2                                    &
                      + 3 * sin(k*(2*x1 + x2))                   &
                          * sin(k - k*x2)                        &
                          * sin(k*(x1 - 3*x2 + 2*x3))            &
                          * sin(k*(3*x1 - 2*x2 + 2*x3))          &
                      )

    grad_u(i,3) = 2 * k * cos(k*(4*x1 - 5*x2 + 4*x3))  &
                        * sin(k*(1 + x1))              &
                        * sin(k*(2*x1 + x2))           &
                        * sin(k - k*x2)

  end do

end subroutine GetExactGradient

!-------------------------------------------------------------------------------
!> Exact laplacian

subroutine GetExactLaplacian(k, n, x, laplace_u)
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

    laplace_u(i) = ( -4 * cos(k - k*x2) * sin(k*(1 + x1))                    &
                        * (     sin(2*k*(2*x1 + x2))                         &
                          + 3 * sin(2*k*(x1 - 3*x2 + 2*x3))                  &
                          - 2 * sin(2*k*(3*x1 - 2*x2 + 2*x3))                &
                          )                                                  &
                   + sin(k - k*x2) * ( -2 * sin(k*(1 + x1))                  &
                                     + 15 * sin(k*(1 + 5*x1 + 2*x2))         &
                                     +  7 * sin(k - 3*k*x1 - 2*k*x2)         &
                                     + 29 * sin(k*(1 - 5*x1 + 4*x2 - 4*x3))  &
                                     + 27 * sin(k*(-1 + x1 - 6*x2 + 4*x3))   &
                                     - 31 * sin(k*(1 + 3*x1 - 6*x2 + 4*x3))  &
                                     + 41 * sin(k*(1 + 7*x1 - 4*x2 + 4*x3))  &
                                     )                                       &
                   ) * k**2 / 4

  end do

end subroutine GetExactLaplacian

!-------------------------------------------------------------------------------
!> Diffusivity

subroutine GetDiffusivity(nu_0, nu_1, k, d, n, x, nu)
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

end subroutine GetDiffusivity

!-------------------------------------------------------------------------------
!> Diffusivity

subroutine GetDiffusivityGradient(nu_1, k, d, n, x, grad_nu)
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

    grad_nu(i,1) = -k * nu_1 * c1 * s2 * s3
    grad_nu(i,2) = -k * nu_1 * s1 * c2 * s3
    grad_nu(i,3) = -k * nu_1 * s1 * s2 * c3

  end do

end subroutine GetDiffusivityGradient

!-------------------------------------------------------------------------------
!> Exact source

subroutine GetSource(lambda, nu_0, nu_1, k_nu, d_nu, k_u, n, x, f)
  real(RNP), intent(in)  :: lambda  !< Helmholtz parameter
  real(RNP), intent(in)  :: nu_0    !< constant part ν₀
  real(RNP), intent(in)  :: nu_1    !< fluctuation amplitude ν₁
  integer,   intent(in)  :: k_nu    !< fluctuation wave number
  real(RNP), intent(in)  :: d_nu    !< fluctuation phase shift
  integer,   intent(in)  :: k_u     !< solution wave number
  integer,   intent(in)  :: n       !< number of points
  real(RNP), intent(in)  :: x(n,3)  !< mesh points
  real(RNP), intent(out) :: f(n)    !< source f(x)

  real(RNP), allocatable :: u(:), grad_u(:,:), laplace_u(:)
  real(RNP), allocatable :: nu(:), grad_nu(:,:)
  integer :: i

  allocate(u(n), grad_u(n,3), laplace_u(n), nu(n), grad_nu(n,3))

  call GetExactSolution(k_u, n, x, u)
  call GetExactGradient(k_u, n, x, grad_u)
  call GetExactLaplacian(k_u, n, x, laplace_u)
  call GetDiffusivity(nu_0, nu_1, k_nu, d_nu, n, x, nu)
  call GetDiffusivityGradient(nu_1, k_nu, d_nu, n, x, grad_nu)

  do i = 1, n

    f(i) = lambda * u(i)                 &
         - ( nu(i) * laplace_u(i)        &
           + grad_nu(i,1) * grad_u(i,1)  &
           + grad_nu(i,2) * grad_u(i,2)  &
           + grad_nu(i,3) * grad_u(i,3) )

  end do

end subroutine GetSource

!===============================================================================

end module Elliptic_Test_Case
