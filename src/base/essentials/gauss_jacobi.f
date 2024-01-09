!> summary:  Implementation of Jacobi polynomials and Gauss quadratures
!> author:   Joerg Stiller
!> date:     2013/05/13
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Implementation of Jacobi polynomials and Gauss quadratures
!>
!> Provides procedures for evaluating
!>
!>   -  the values, derivatives and zeros of Jacobi polynomials,
!>   -  the Gauss-Legendre (G) points and weights, and the related
!>      Lagrange polynomials and their derivatives,
!>   -  the Gauss-Lobatto-Legendre (L) points and weights, and the related
!>      Lagrange polynomials and their derivatives.
!>   -  the left/right-sided Gauss-Radau-Legendre (RL/RR) points and weights,
!>      and the related Lagrange polynomials and their derivatives.
!>
!> Implementation follows G.E. Karniadakis & S.J. Sherwin,
!> Spectral/hp Element Methods for CFD. Oxford Univ. Press, 2005.
!==============================================================================

module Gauss_Jacobi
  use Kind_parameters, only: RNP
  use Constants,       only: ZERO, HALF, ONE, TWO, FOUR, PI
  implicit none
  private

  public :: JacobiPolynomial
  public :: JacobiPolynomialDerivative
  public :: JacobiPolynomialZeros

  public :: GaussPoints
  public :: GaussWeights
  public :: GaussPolynomial
  public :: GaussPolynomialDerivative
  public :: GaussDiffMatrix

  public :: LobattoPoints
  public :: LobattoWeights
  public :: LobattoPolynomial
  public :: LobattoPolynomialDerivative
  public :: LobattoDiffMatrix

  public :: RadauPoints
  public :: RadauWeights
  public :: RadauPolynomial
  public :: RadauPolynomialDerivative
  public :: RadauDiffMatrix

  !> Mininmal admissible distance distance to collocation points.
  !> If the distance of a given point falls below `TOL`, the closest
  !> collocation point is used instead.
  real(RNP), parameter :: TOL = 1.0e-10_RNP

contains

  !=============================================================================
  ! Jacobi polynomials

  !-----------------------------------------------------------------------------
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

  !-----------------------------------------------------------------------------
  !> Returns the derivative of the Jacobi polynomial P^{a,b}_{n}(x)

  pure function JacobiPolynomialDerivative(a, b, n, x) result(dy)
    real(RNP), intent(in) :: a   !< first exponent
    real(RNP), intent(in) :: b   !< second exponent
    integer,   intent(in) :: n   !< order of the Jacobi polynomial
    real(RNP), intent(in) :: x   !< position in [-1,1]
    real(RNP)             :: dy  !< first derivative of P^{a,b}_{n} at x

    dy = HALF * (a + b + n + 1) * JacobiPolynomial(a+1, b+1, n-1, x)

  end function JacobiPolynomialDerivative

  !-----------------------------------------------------------------------------
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

  !=============================================================================
  ! Gauss-Legendre quadrature and related Lagrangre polynomials

  !-----------------------------------------------------------------------------
  !> Gauss-Legendre quadrature points of degree q in [-1,1]

  pure function GaussPoints(q) result(x)
    integer, intent(in) :: q       !< degree of the quadrature polynomial
    real(RNP)           :: x(0:q)  !< quadrature points

    x = JacobiPolynomialZeros(ZERO, ZERO, n=q+1)

  end function GaussPoints

  !-----------------------------------------------------------------------------
  !> Quadrature weights related to the Gauss-Legendre points x(0:q)

  pure function GaussWeights(x) result(w)
    real(RNP), intent(in) :: x(0:)            !< quadrature points
    real(RNP)             :: w(0:ubound(x,1)) !< quadrature weights

    integer :: i, q

    q = ubound(x,1)
    forall(i = 0:q)
      w(i) = 2/((1-x(i)**2) * JacobiPolynomialDerivative(ZERO,ZERO,q+1,x(i))**2)
    end forall

  end function GaussWeights

  !-----------------------------------------------------------------------------
  !> Lagrange polynomial l_k to Gauss-Legendre points xi(0:n)

  pure function GaussPolynomial(k, xi, x) result(y)
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

  end function GaussPolynomial

  !-----------------------------------------------------------------------------
  !> Derivative of the Gauss Lagrange polynomial l_k at the p-th Gauss point

  pure function GaussPolynomialDerivative(k, xi, p) result(dy)
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

  end function GaussPolynomialDerivative

  !-----------------------------------------------------------------------------
  !> Returns the Gauss-Legendre differentiation matrix
  !>
  !> Let xi(0:q) denote the Gauss points of degree q and l_j(x) the Lagrange
  !> polynomial associated with xi(j). The colllocation differentiation matrix
  !> is then defined as
  !>
  !>     D(i,j) = l_j'(xi(i))    for 0 <= i,j <= q

  pure function GaussDiffMatrix(xi) result(D)
    real(RNP), intent(in) :: xi(0:)                !< Gauss points
    real(RNP) :: D(0:ubound(xi,1),0:ubound(xi,1))  !< Gauss diff matrix

    integer :: i, j

    forall(i = 0:ubound(xi,1), j = 0:ubound(xi,1))
      D(i,j) = GaussPolynomialDerivative(j, xi, i)
    end forall

  end function GaussDiffMatrix

  !=============================================================================
  ! Gauss-Lobatto-Legendre quadrature and related Lagrangre polynomials

  !-----------------------------------------------------------------------------
  !> Gauss-Lobatto-Legendre quadrature points of degree q

  pure function LobattoPoints(q) result(x)
    integer, intent(in) :: q       !< degree of the quadrature polynomial
    real(RNP)           :: x(0:q)  !< quadrature points

    x(0)     = -1
    x(1:q-1) =  JacobiPolynomialZeros(ONE, ONE, n=q-1)
    x(q)     =  1

  end function LobattoPoints

  !-----------------------------------------------------------------------------
  !> Quadrature weights at Gauss-Lobatto-Legendre points x(0:q)

  pure function LobattoWeights(x) result(w)
    real(RNP), intent(in) :: x(0:)             !< quadrature points
    real(RNP)             :: w(0:ubound(x,1))  !< quadrature weights

    integer :: i, q

    q = ubound(x,1)
    forall(i = 0:q)
      w(i) = 2 / (q*(q+1) * JacobiPolynomial(ZERO, ZERO, q, x(i))**2)
    end forall

  end function LobattoWeights

  !-----------------------------------------------------------------------------
  !> Lagrange polynomial l_k to Gauss-Lobatto-Legendre points xi(0:n)

  pure function LobattoPolynomial(k, xi, x) result(y)
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

  end function LobattoPolynomial

  !-----------------------------------------------------------------------------
  !> Derivative of the Lobatto Lagrange polynomial l_k at the p-th Lobatto point

  pure function LobattoPolynomialDerivative(k, xi, p) result(dy)
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

  end function LobattoPolynomialDerivative

  !-----------------------------------------------------------------------------
  !> Returns the Gauss-Lobatto-Legendre differentiation matrix
  !>
  !> Let xi(0:q) denote the Lobatto points of degree q and l_j(x) the Lagrange
  !> polynomial associated with xi(j). The colllocation differentiation matrix
  !> is then defined as
  !>
  !>     D(i,j) = l_j'(xi(i))    for 0 <= i,j <= q

  pure function LobattoDiffMatrix(xi) result(D)
    real(RNP), intent(in) :: xi(0:)                !< Lobatto points
    real(RNP) :: D(0:ubound(xi,1),0:ubound(xi,1))  !< Lobatto diff matrix

    integer :: i, j

    forall(i = 0:ubound(xi,1), j = 0:ubound(xi,1))
      D(i,j) = LobattoPolynomialDerivative(j, xi, i)
    end forall

  end function LobattoDiffMatrix

  !=============================================================================
  ! Gauss-Radau-Legendre quadrature and related Lagrangre polynomials

  !-----------------------------------------------------------------------------
  !> Gauss-Radau-Legendre quadrature points of degree q

  pure function RadauPoints(q, right) result(x)
    integer,           intent(in) :: q      !< degree of quadrature polynomial
    logical, optional, intent(in) :: right  !< switch to Radau right [F]
    real(RNP) :: x(0:q)  !< quadrature points

    logical :: left

    if (present(right)) then
      left = .not. right
    else
      left = .true.
    end if

    if (left) then
      x(0)     = -1
      x(1:q)   =  JacobiPolynomialZeros(ZERO, ONE, n=q)
    else
      x(0:q-1) =  JacobiPolynomialZeros(ONE, ZERO, n=q)
      x(q)     =  1
    end if

  end function RadauPoints

  !-----------------------------------------------------------------------------
  !> Quadrature weights at Gauss-Radau-Legendre points x(0:q)

  pure function RadauWeights(x) result(w)
    real(RNP), intent(in) :: x(0:)             !< quadrature points
    real(RNP)             :: w(0:ubound(x,1))  !< quadrature weights

    integer :: i, q

    q = ubound(x,1)

    if (abs(ONE + x(0)) < TOL) then
      ! Radau left
      w(0) = TWO / (q+1)**2
      do i = 1, q
        w(i) = (1 - x(i)) / ((q+1) * JacobiPolynomial(ZERO, ZERO, q, x(i)))**2
      end do
    else
      ! Radau right
      do i = 0, q-1
        w(i) = (1 + x(i)) / ((q+1) * JacobiPolynomial(ZERO, ZERO, q, x(i)))**2
      end do
      w(q) = TWO / (q+1)**2
    end if

  end function RadauWeights

  !-----------------------------------------------------------------------------
  !> Lagrange polynomial l_k to Gauss-Radau-Legendre points xi(0:n)

  pure function RadauPolynomial(k, xi, x) result(y)
    integer,   intent(in) :: k       !< polynomial ID, 0 <= k <= q = ubound(x,1)
    real(RNP), intent(in) :: xi(0:)  !< Gauss-Radau-Legendre points
    real(RNP), intent(in) :: x       !< position in [-1, 1]
    real(RNP)             :: y       !< y = l_k(x)

    real(RNP) :: x_, xi_k
    integer   :: q

    if (abs(x - xi(k)) < TOL) then

      y = 1

    else

      q = ubound(xi,1)

      if (abs(ONE + xi(0)) < TOL) then
        ! Radau left
        x_   =  x
        xi_k =  xi(k)
      else
        ! Radau right
        x_   = -x
        xi_k = -xi(k)
      end if

      y = (1 + x_) * JacobiPolynomial(ZERO, ONE, q, x_)                    &
        / ( ( (1 + xi_k) * JacobiPolynomialDerivative(ZERO, ONE, q, xi_k)  &
            + JacobiPolynomial(ZERO, ONE, q, xi_k)                         &
            ) * (x_ - xi_k)                                                &
          )

    end if

  end function RadauPolynomial

  !-----------------------------------------------------------------------------
  !> Derivative of the Radau Lagrange polynomial l_k at the p-th Radau point

  pure function RadauPolynomialDerivative(k, xi, p) result(dy)
    integer,   intent(in) :: k       !< polynomial ID, 0 <= k <= q = ubound(x,1)
    real(RNP), intent(in) :: xi(0:)  !< Gauss-Radau-Legendre points
    integer,   intent(in) :: p       !< point ID, 0 <= p <= q
    real(RNP)             :: dy      !< dy = pi'_k(xi(p))

    real(RNP) :: xi_k, xi_p
    integer   :: q, s
    logical   :: left

    left = abs(ONE + xi(0)) < TOL
    if (left) then
      s =  1
    else
      s = -1
    end if

    q = ubound(xi,1)

    if (s == 1 .and. k == 0 .and. p == 0) then

      dy = -q * (q+2) / FOUR

    else if (s == -1 .and. k == q .and. p == q) then

      dy =  q * (q+2) / FOUR

    else if (k == p) then

      dy =  s / (2 * (1 - s*xi(p)))

    else

      xi_k = s * xi(k)
      xi_p = s * xi(p)
      dy = s * JacobiPolynomial(ZERO, ZERO, q, xi_p) * (1-xi_k) &
         / (JacobiPolynomial(ZERO, ZERO, q, xi_k) * (1-xi_p) * (xi_p-xi_k))

    end if

  end function RadauPolynomialDerivative

  !-----------------------------------------------------------------------------
  !> Returns the Gauss-Radau-Legendre differentiation matrix
  !>
  !> Let xi(0:q) denote the Lobatto points of degree q and l_j(x) the Lagrange
  !> polynomial associated with xi(j). The colllocation differentiation matrix
  !> is then defined as
  !>
  !>     D(i,j) = l_j'(xi(i))    for 0 <= i,j <= q

  pure function RadauDiffMatrix(xi) result(D)
    real(RNP), intent(in) :: xi(0:)                !< Radau points
    real(RNP) :: D(0:ubound(xi,1),0:ubound(xi,1))  !< Radau diff matrix

    integer :: i, j

    forall(i = 0:ubound(xi,1), j = 0:ubound(xi,1))
      D(i,j) = RadauPolynomialDerivative(j, xi, i)
    end forall

  end function RadauDiffMatrix

  !=============================================================================

end module Gauss_Jacobi
