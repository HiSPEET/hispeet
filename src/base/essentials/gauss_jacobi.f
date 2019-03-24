!> summary:  Implementation of Jacobi polynomials and Gauss quadratures
!> author:   Joerg Stiller
!> date:     2013/05/13
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Implementation of Jacobi polynomials and Gauss quadratures
!>
!> Provides procedures for evaluating
!>
!>   *  the values, derivatives and zeros of Jacobi polynomials,
!>   *  the Gauss-Legendre (GL) points and weights, and the related
!>      Lagrange polynomials and their derivatives,
!>   *  the Gauss-Lobatto-Legendre (GLL) points and weights, and the
!>      related Lagrange polynomials and their derivatives.
!>   *  the left-sided Gauss-Radau-Legendre (GRL) points and weights,
!>      and the related Lagrange polynomials and their derivatives.
!>
!> Implementation follows G.E. Karniadakis & S.J. Sherwin,
!> Spectral/hp Element Methods for CFD. Oxford Univ. Press, 2005.
!==============================================================================

module Gauss_Jacobi
  use Kind_parameters, only: RNP
  use Constants,       only: ZERO, HALF, ONE, FOUR, PI
  implicit none
  private

  public :: JacobiPolynomial
  public :: JacobiPolynomialDerivative
  public :: JacobiPolynomialZeros

  public :: GL_Points
  public :: GL_Weights
  public :: GL_Polynomial
  public :: GL_PolynomialDerivative
  public :: GL_DiffMatrix

  public :: GLL_Points
  public :: GLL_Weights
  public :: GLL_Polynomial
  public :: GLL_PolynomialDerivative
  public :: GLL_DiffMatrix

  public :: GRL_Points
  public :: GRL_Weights
  public :: GRL_Polynomial
  public :: GRL_PolynomialDerivative
  public :: GRL_DiffMatrix

  !> Mininmal admissible distance distance to collocation points.
  !> If the distance of a given point falls below `TOL`, the closest
  !> collocation point is used instead.
  real(RNP), parameter :: TOL = 1.0e-10_RNP

contains

!===============================================================================
! Jacobi polynomials

!-------------------------------------------------------------------------------
!> Returns the Jacobi polynomial P^{a,b}_{n} at position x in [-1,1]

pure function JacobiPolynomial(a, b, n, x) result(y)
  real(RNP), intent(in) :: a  !< first exponent, a = α
  real(RNP), intent(in) :: b  !< second exponent, b = β
  integer,   intent(in) :: n  !< order of the Jacobi polynomial
  real(RNP), intent(in) :: x  !< position in [-1,1]
  real(RNP)             :: y  !< P^{a,b}_{n}(x)

  real(RNP) :: P0, P1, a1, a2, a3, a4
  integer   :: i

  P0 = 1
  P1 = (a-b+(a+b+2)*x) / 2

  if (n < 2) then
    y = (1-n)*P0 + n*P1

  else if (a > -1 .and. b > -1) then
    do i = 1, n-1
      a1 = 2 * (i+1) * (i+a+b+1) * (2*i+a+b)
      a2 = (2*i+a+b+1) * (a*a-b*b)
      a3 = (2*i+a+b) * (2*i+a+b+1) * (2*i+a+b+2)
      a4 = 2 * (i+a)*(i+b) * (2*i+a+b+2)
      y  = ((a2+a3*x)*P1 - a4*P0)/a1
      P0 = P1
      P1 = y
    end do

  else ! a or b out of range, complain by returning very large value ;-)
    y = huge(y)
  end if

end function JacobiPolynomial

!-------------------------------------------------------------------------------
!> Returns the derivative of the Jacobi polynomial P^{a,b}_{n}(x)

pure function JacobiPolynomialDerivative(a, b, n, x) result(dy)
  real(RNP), intent(in) :: a   !< first exponent
  real(RNP), intent(in) :: b   !< second exponent
  integer,   intent(in) :: n   !< order of the Jacobi polynomial
  real(RNP), intent(in) :: x   !< position in [-1,1]
  real(RNP)             :: dy  !< first derivative of P^{a,b}_{n} at x

  dy = HALF * (a + b + n + 1) * JacobiPolynomial(a+1, b+1, n-1, x)

end function JacobiPolynomialDerivative

!-------------------------------------------------------------------------------
!> Returns the zeros of the Jacobi polynomial P^{a,b}_{n}(x)

pure function JacobiPolynomialZeros(a, b, n) result(x)
  real(RNP), intent(in)  :: a        !< first exponent
  real(RNP), intent(in)  :: b        !< second exponent
  integer,   intent(in)  :: n        !< order of the Jacobi polynomial
  real(RNP)              :: x(0:n-1) !< zeros in ascending order

  integer,   parameter :: MAXIT = 100          ! maximum number of iterations
  real(RNP), parameter :: EPS   = epsilon(ONE) ! accuracy threshold
  real(RNP) :: JP, JD, r, s, d
  integer   :: i, j

  do i= 0, n-1
    r = -real(cos(PI*real(2*i+1,RNP) / real(2*n,RNP)), RNP)
    if (i > 0) r = HALF * (r + x(i-1))
    do j = 1, MAXIT
      if (i > 0) then
        s = sum(1/(r - x(0:i-1)))
      else
        s = 1
      end if
      JP = JacobiPolynomial(a, b, n, r)
      JD = JacobiPolynomialDerivative(a, b, n, r)
      d = JP / (JD - s*JP)
      r = r - d
      if (abs(d) < EPS) exit
    end do
    x(i) = r
  end do

end function JacobiPolynomialZeros

!===============================================================================
! Gauss-Legendre quadrature and related Lagrangre polynomials

!-------------------------------------------------------------------------------
!> Gauss-Legendre quadrature points of degree q in [-1,1]

pure function GL_Points(q) result(x)
  integer, intent(in) :: q       !< degree of the quadrature polynomial
  real(RNP)           :: x(0:q)  !< quadrature points

  x = JacobiPolynomialZeros(ZERO, ZERO, n=q+1)

end function GL_Points

!-------------------------------------------------------------------------------
!> Quadrature weights related to the Gauss-Legendre points x(0:q)

pure function GL_Weights(x) result(w)
  real(RNP), intent(in) :: x(0:)            !< quadrature points
  real(RNP)             :: w(0:ubound(x,1)) !< quadrature weights

  integer :: i, q

  q = ubound(x,1)
  forall(i = 0:q)
    w(i) = 2/((1-x(i)**2) * JacobiPolynomialDerivative(ZERO,ZERO,q+1,x(i))**2)
  end forall

end function GL_Weights

!-------------------------------------------------------------------------------
!> Lagrange polynomial l_k to Gauss-Legendre points xi(0:n)

pure function GL_Polynomial(k, xi, x) result(y)
  integer,   intent(in) :: k       !< polynomial ID, 0 <= k <= q = ubound(x,1)
  real(RNP), intent(in) :: xi(0:)  !< Gauss-Legendre points
  real(RNP), intent(in) :: x       !< position in [-1, 1]
  real(RNP)             :: y       !< y = l_k(x)

  integer :: q

  q = ubound(xi,1)
  if (abs(x - xi(k)) < TOL) then
    y = 1
  else
    y = JacobiPolynomial(ZERO, ZERO, q+1, x)  &
      / (JacobiPolynomialDerivative(ZERO, ZERO, q+1, xi(k)) * (x-xi(k)))
  end if

end function GL_Polynomial

!-------------------------------------------------------------------------------
!> Derivative of the GL Lagrange polynomial l_k at the p-th GL point

pure function GL_PolynomialDerivative(k, xi, p) result(dy)
  integer,   intent(in) :: k       !< polynomial ID, 0 <= k <= q = ubound(x,1)
  real(RNP), intent(in) :: xi(0:)  !< Gauss-Legendre points
  integer,   intent(in) :: p       !< point ID, 0 <= p <= q
  real(RNP)             :: dy      !< dy = pi'_k(xi(p))

  integer :: q

  q = ubound(xi,1)
  if (k == p) then
    dy =  xi(p) / (1 - xi(p)**2)
  else
    dy =  JacobiPolynomialDerivative(ZERO, ZERO, q+1, xi(p)) &
       / (JacobiPolynomialDerivative(ZERO, ZERO, q+1, xi(k)) * (xi(p)-xi(k)))
  end if

end function GL_PolynomialDerivative

!-------------------------------------------------------------------------------
!> Returns the Gauss-Legendre differentiation matrix
!>
!> Let xi(0:q) denote the GL points of degree q and l_j(x) the Lagrange
!> polynomial associated with xi(j). The colllocation differentiation matrix
!> is then defined as
!>
!>     D(i,j) = l_j'(xi(i))    for 0 <= i,j <= q

pure function GL_DiffMatrix(xi) result(D)
  real(RNP), intent(in) :: xi(0:)                !< GL points
  real(RNP) :: D(0:ubound(xi,1),0:ubound(xi,1))  !< GL differentiation matrix

  integer :: i, j

  forall(i = 0:ubound(xi,1), j = 0:ubound(xi,1))
    D(i,j) = GL_PolynomialDerivative(j, xi, i)
  end forall

end function GL_DiffMatrix

!===============================================================================
! Gauss-Lobatto-Legendre quadrature and related Lagrangre polynomials

!-------------------------------------------------------------------------------
!> Gauss-Lobatto-Legendre quadrature points of degree q

pure function GLL_Points(q) result(x)
  integer, intent(in) :: q       !< degree of the quadrature polynomial
  real(RNP)           :: x(0:q)  !< quadrature points

  x(0)     = -1
  x(1:q-1) =  JacobiPolynomialZeros(ONE, ONE, n=q-1)
  x(q)     =  1

end function GLL_Points

!-------------------------------------------------------------------------------
!> Quadrature weights at Gauss-Lobatto-Legendre points x(0:q)

pure function GLL_Weights(x) result(w)
  real(RNP), intent(in) :: x(0:)             !< quadrature points
  real(RNP)             :: w(0:ubound(x,1))  !< quadrature weights

  integer :: i, q

  q = ubound(x,1)
  forall(i = 0:q)
    w(i) = 2 / (q*(q+1) * JacobiPolynomial(ZERO, ZERO, q, x(i))**2)
  end forall

end function GLL_Weights

!-------------------------------------------------------------------------------
!> Lagrange polynomial l_k to Gauss-Lobatto-Legendre points xi(0:n)

pure function GLL_Polynomial(k, xi, x) result(y)
  integer,   intent(in) :: k       !< polynomial ID, 0 <= k <= q = ubound(x,1)
  real(RNP), intent(in) :: xi(0:)  !< Gauss-Lobatto-Legendre points
  real(RNP), intent(in) :: x       !< position in [-1, 1]
  real(RNP)             :: y       !< y = l_k(x)

  integer :: q

  q = ubound(xi,1)
  if (abs(x - xi(k)) < TOL) then
    y = 1
  else
    y = (x-1)*(x+1) * JacobiPolynomialDerivative(ZERO, ZERO, q, x)  &
      / (q*(q+1) * JacobiPolynomial(ZERO, ZERO, q, xi(k)) * (x-xi(k)))
  end if

end function GLL_Polynomial

!-------------------------------------------------------------------------------
!> Derivative of the GLL Lagrange polynomial l_k at the p-th GLL point

pure function GLL_PolynomialDerivative(k, xi, p) result(dy)
  integer,   intent(in) :: k       !< polynomial ID, 0 <= k <= q = ubound(x,1)
  real(RNP), intent(in) :: xi(0:)  !< Gauss-Lobatto-Legendre points
  integer,   intent(in) :: p       !< point ID, 0 <= p <= q
  real(RNP)             :: dy      !< dy = pi'_k(xi(p))

  integer :: q

  q = ubound(xi,1)
  if      (k == 0 .and. p == 0) then;  dy = -q * (q+1) / FOUR
  else if (k == q .and. p == q) then;  dy =  q * (q+1) / FOUR
  else if (k == p)              then;  dy =  0
  else
    dy =  JacobiPolynomial(ZERO, ZERO, q, xi(p)) &
       / (JacobiPolynomial(ZERO, ZERO, q, xi(k)) * (xi(p)-xi(k)))
  end if

end function GLL_PolynomialDerivative

!-------------------------------------------------------------------------------
!> Returns the Gauss-Lobatto-Legendre differentiation matrix
!>
!> Let xi(0:q) denote the GLL points of degree q and l_j(x) the Lagrange
!> polynomial associated with xi(j). The colllocation differentiation matrix
!> is then defined as
!>
!>     D(i,j) = l_j'(xi(i))    for 0 <= i,j <= q

pure function GLL_DiffMatrix(xi) result(D)
  real(RNP), intent(in) :: xi(0:)                !< GLL points
  real(RNP) :: D(0:ubound(xi,1),0:ubound(xi,1))  !< GLL differentiation matrix

  integer :: i, j

  forall(i = 0:ubound(xi,1), j = 0:ubound(xi,1))
    D(i,j) = GLL_PolynomialDerivative(j, xi, i)
  end forall

end function GLL_DiffMatrix

!===============================================================================
! Gauss-Radau-Legendre quadrature and related Lagrangre polynomials

!-------------------------------------------------------------------------------
!> Left-sided Gauss-Radau-Legendre quadrature points of degree q

pure function GRL_Points(q) result(x)
  integer, intent(in) :: q       !< degree of the quadrature polynomial
  real(RNP)           :: x(0:q)  !< quadrature points

  x(0)   = -1
  x(1:q) =  JacobiPolynomialZeros(ZERO, ONE, n=q)

end function GRL_Points

!-------------------------------------------------------------------------------
!> Quadrature weights at Gauss-Radau-Legendre points x(0:q)

pure function GRL_Weights(x) result(w)
  real(RNP), intent(in) :: x(0:)             !< quadrature points
  real(RNP)             :: w(0:ubound(x,1))  !< quadrature weights

  integer :: i, q

  q = ubound(x,1)
  forall(i = 0:q)
    w(i) = (1 - x(i)) / ((q+1) * JacobiPolynomial(ZERO, ZERO, q, x(i)))**2
  end forall

end function GRL_Weights

!-------------------------------------------------------------------------------
!> Lagrange polynomial l_k to Gauss-Radau-Legendre points xi(0:n)

pure function GRL_Polynomial(k, xi, x) result(y)
  integer,   intent(in) :: k       !< polynomial ID, 0 <= k <= q = ubound(x,1)
  real(RNP), intent(in) :: xi(0:)  !< Gauss-Radau-Legendre points
  real(RNP), intent(in) :: x       !< position in [-1, 1]
  real(RNP)             :: y       !< y = l_k(x)

  integer :: q

  q = ubound(xi,1)
  if (abs(x - xi(k)) < TOL) then
    y = 1
  else
    y = (x+1) * JacobiPolynomial(ZERO, ONE, q, x)                          &
      / ( ( (xi(k) + 1) * JacobiPolynomialDerivative(ZERO, ONE, q, xi(k))  &
          + JacobiPolynomial(ZERO, ONE, q, xi(k))                          &
          ) * (x - xi(k))                                                  &
        )
  end if

end function GRL_Polynomial

!-------------------------------------------------------------------------------
!> Derivative of the GRL Lagrange polynomial l_k at the p-th GL point

pure function GRL_PolynomialDerivative(k, xi, p) result(dy)
  integer,   intent(in) :: k       !< polynomial ID, 0 <= k <= q = ubound(x,1)
  real(RNP), intent(in) :: xi(0:)  !< Gauss-Radau-Legendre points
  integer,   intent(in) :: p       !< point ID, 0 <= p <= q
  real(RNP)             :: dy      !< dy = pi'_k(xi(p))

  integer :: q

  q = ubound(xi,1)
  if      (k == 0 .and. p == 0) then;  dy = -q * (q+2) / FOUR
  else if (k == p)              then;  dy =  1 / (2 * (1 - xi(p)))
  else
    dy =  JacobiPolynomial(ZERO, ZERO, q, xi(p)) * (1 - xi(k)) &
       / (JacobiPolynomial(ZERO, ZERO, q, xi(k)) * (1 - xi(p)) * (xi(p)-xi(k)))
  end if

end function GRL_PolynomialDerivative

!-------------------------------------------------------------------------------
!> Returns the Gauss-Radau-Legendre differentiation matrix
!>
!> Let xi(0:q) denote the GLL points of degree q and l_j(x) the Lagrange
!> polynomial associated with xi(j). The colllocation differentiation matrix
!> is then defined as
!>
!>     D(i,j) = l_j'(xi(i))    for 0 <= i,j <= q

pure function GRL_DiffMatrix(xi) result(D)
  real(RNP), intent(in) :: xi(0:)                !< GRL points
  real(RNP) :: D(0:ubound(xi,1),0:ubound(xi,1))  !< GRL differentiation matrix

  integer :: i, j

  forall(i = 0:ubound(xi,1), j = 0:ubound(xi,1))
    D(i,j) = GRL_PolynomialDerivative(j, xi, i)
  end forall

end function GRL_DiffMatrix

!===============================================================================

end module Gauss_Jacobi