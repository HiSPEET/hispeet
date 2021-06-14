!> summary:  Lagrange interpolation, including differentiation and integration
!> author:   Joerg Stiller
!> date:     2019/03/17
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Lagrange_Interpolation
  use Kind_parameters, only: RNP
  use Gauss_Jacobi
  implicit none
  private

  public :: LagrangePolynomial
  public :: LagrangePolynomialDerivative
  public :: LagrangeDiffMatrix
  public :: GaussLagrangeWeights

contains

!-------------------------------------------------------------------------------
!> Returns the k-th Lagrange polynomial to points xc(0:n) at position x

pure function LagrangePolynomial(k, xc, x) result(y)
  integer,   intent(in) ::  k       !< polynomial ID, 0 ≤ k ≤ n
  real(RNP), intent(in) ::  xc(0:)  !< collocation points, xc(i-1) < xc(i)
  real(RNP), intent(in) ::  x       !< evaluation position
  real(RNP)             ::  y

  integer :: i

  y = 1
  do i = 0, ubound(xc,1)
    if (i == k) cycle
    y = y * (x - xc(i)) / (xc(k) - xc(i))
  end do

end function LagrangePolynomial

!-------------------------------------------------------------------------------
!> Derivative of the k-th Lagrange polynomial to points xc(0:n) at xc(p)

pure function LagrangePolynomialDerivative(k, xc, p) result(dy)
  integer,   intent(in) ::  k       !< polynomial ID, 0 ≤ k ≤ n
  real(RNP), intent(in) ::  xc(0:)  !< collocation points, xc(i-1) < xc(i)
  integer,   intent(in) ::  p       !< point ID, 0 ≤ p ≤ n
  real(RNP)             ::  dy

  integer :: i

  if (k /= p) then
    dy = 1 / (xc(k) - xc(p))
    do i = 0, ubound(xc,1)
      if (i == k .or. i == p) cycle
      dy = dy * (xc(p) - xc(i)) / (xc(k) - xc(i))
    end do
  else
    dy = 0
    do i = 0, ubound(xc,1)
      if (i == k) cycle
      dy = dy + 1 / (xc(k) - xc(i))
    end do
  end if

end function LagrangePolynomialDerivative

!-------------------------------------------------------------------------------
!> Returns the Lagrange differentiation matrix
!>
!> Let xc(0:n) denote the collocation points l_j(x) the Lagrange polynomial
!> associated with xc(j). The colllocation differentiation matrix is then
!> defined as
!>
!>     D(i,j) = l_j'(xc(i))    for 0 <= i,j <= n

pure function LagrangeDiffMatrix(xc) result(D)
  real(RNP), intent(in) :: xc(0:)  !< collocation points, xc(i-1) < xc(i)
  real(RNP) :: D(0:ubound(xc,1),0:ubound(xc,1))

  integer :: i, j, n

  n = ubound(xc,1)
  do i = 0, n
  do j = 0, n
    D(i,j) = LagrangePolynomialDerivative(j, xc, i)
  end do
  end do

end function LagrangeDiffMatrix

!-------------------------------------------------------------------------------
!> Returns Gauss quadrature weights for integrating the Lagrange interpolant
!>
!> The collocation points `xcᵢ` define Lagrange polynomials `𝓁ᵢ(x)` which
!> interpolate any given set of coefficients (values) `ucᵢ`. Mapping the
!> interval `[xc(0), xc(n)]` to the reference domain [-1,1], the integral
!> over the Langrange interpolant can be evaluated using a Gauss quadrature
!> of degree `q = n/2`. This yields
!>
!>     ∫ ∑ᵢ ucᵢ 𝓁ᵢ(x) dx ≃ ∑ᵢ wcᵢ ucᵢ
!>
!> with weights
!>
!>     wcᵢ = L/2 sum(k=0:q) wᵏ 𝓁ᵢ(x(ξᵏ))
!>
!> where `ξᵏ` are the Gauss points in [-1,1], `wᵏ` the corresponding weights
!> and
!>
!>     x(ξ) = [ x(0) + x(n)  +  ξ (x(n) - x(0)) ] / 2

pure function GaussLagrangeWeights(xc) result(wc)
  real(RNP), intent(in) :: xc(0:)  !< collocation points, xc(i-1) < xc(i)
  real(RNP) :: wc(0:ubound(xc,1))

  real(RNP), allocatable :: xq(:), wq(:), c(:,:)
  real(RNP) :: a, b
  integer   :: n, q
  integer   :: i, k

  n = ubound(xc,1)
  q = n / 2

  ! Gauss points and weights in [-1,1]
  allocate(xq(0:q), source = GaussPoints(q))      ! xq = ξ
  allocate(wq(0:q), source = GaussWeights(xq))    ! wq = w

  ! evaluate Lagrange polynomials at Gauss points
  allocate(c(0:n,0:q))
  a = (xc(n) + xc(0)) / 2
  b = (xc(n) - xc(0)) / 2
  do k = 0, q
  do i = 0, n
    c(i,k) = b * LagrangePolynomial(i, xc, x = a + b*xq(k))
  end do
  end do

  ! scaled Gauss weights at collocation points
  wc = matmul(c, wq)

end function GaussLagrangeWeights

!===============================================================================

end module Lagrange_Interpolation
