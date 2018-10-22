!> summary:  Test case for elliptic equations, 2π-periodic
!> author:   Joerg Stiller
!> date:     2013/06/16, revised 2017/05/06
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
!> \brief   Exact laplacian
!> \author  Joerg Stiller

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

!===============================================================================

end module Elliptic_Test_Case
