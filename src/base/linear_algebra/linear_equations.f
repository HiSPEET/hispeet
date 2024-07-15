!> summary:   Solvers for small linear equation systems.
!> author:    Joerg Stiller, Immo Huismann, Lars Haupt
!> date:      2013/06/26
!> license:   Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Solvers for small linear equation systems
!==============================================================================

module Linear_Equations
  use Kind_Parameters, only: RNP, RDP, RHP
  implicit none
  private

  public :: TridiagonalSolver
  public :: CyclicTridiagonalSolver
  public :: LU_Decomposition
  public :: LU_Solver
  public :: CholeskyDecomposition
  public :: CholeskySolver

  interface LU_Decomposition
    module procedure LU_Decomposition_Pivot
    module procedure LU_Decomposition_NoPivot
  end interface LU_Decomposition

  interface LU_Solver
    module procedure LU_Solver1
  end interface LU_Solver

  interface CholeskySolver
    module procedure CholeskySolver1
    module procedure CholeskySolver2
  end interface CholeskySolver

contains

  !-----------------------------------------------------------------------------
  !> Gauss method for tridiagonal systems
  !>
  !> Computes the solution of the tridiagonal system
  !>
  !>                     b(1)*x(1)  +  c(1)*x(2)    =  f(1),
  !>     a(i)*x(i-1)  +  b(i)*x(i)  +  c(i)*x(i+1)  =  f(i),  1 < i < n
  !>     a(n)*x(n-1)  +  b(n)*x(n)                  =  f(n)
  !>
  !> Input
  !>     a, b, c   : coefficients
  !>     x         : RHS, i.e., x=f
  !>     regular   : .true. to indicate that the system matrix is regular, \n
  !>                 .false. for singular system of rank n-1
  !>
  !> Output
  !>     x        : solution
  !>     a, b, c  : destroyed

  pure subroutine TridiagonalSolver(x, a, b, c, regular)
    real(RNP), intent(inout) :: x(:)    !< f on input, x on output
    real(RNP), intent(inout) :: a(:)    !< coefficients a(i), will be destroyed!
    real(RNP), intent(inout) :: b(:)    !< coefficients b(i), will be destroyed!
    real(RNP), intent(inout) :: c(:)    !< coefficients c(i), will be destroyed!
    logical,   intent(in)    :: regular !< true/false for regular/irregular case

    real(RNP) :: h, d, p, q, z
    integer   :: i, n

    n = size(x)

    ! factorization ............................................................

    if (abs(b(1)) < epsilon(b)) return
    h = 1 / b(1)
    d = h*x(1)
    p =-h*c(1)
    a(1) = d
    b(1) = p
    do i = 2, n-1
       q = b(i) + a(i)*b(i-1)
       if (abs(q) < epsilon(q)) then
         error stop 'TridiagonalSolver: failed'
       end if
       h    = 1 / q
       b(i) =-h*c(i)
       a(i) = h*(x(i) - a(i)*a(i-1))
       d    = p*a(i) + d
       p    = p*b(i)
    end do
    h = b(n-1)

    ! solution .................................................................

    if (regular) then
       q = b(n) + a(n)*h
       if (abs(q) < epsilon(q)) return
       z = (x(n) - a(n)*a(n-1)) / q
    else
       z = 0
    end if
    x(n) = z
    x(n-1) = a(n-1) + h*z
    do i = n-2, 1, -1
       x(i) = a(i) + b(i)*x(i+1)
    end do

  end subroutine TridiagonalSolver

  !-----------------------------------------------------------------------------
  !> Gauss method for cyclic tridiagonal systems
  !>
  !> Computes the solution of the tridiagonal system
  !>
  !>     a(1)*x(n)    +  b(1)*x(1)  +  c(1)*x(2)    =  f(1),
  !>     a(i)*x(i-1)  +  b(i)*x(i)  +  c(i)*x(i+1)  =  f(i),  1 < i < n
  !>     a(n)*x(n-1)  +  b(n)*x(n)  +  c(n)*x(1)    =  f(n)
  !>
  !> Input
  !>     a, b, c  : coefficients
  !>     x        : RHS, i.e., x=f
  !>     regular  : .true.  to indicate that the system matrix is regular,
  !>                .false. for singular system of rank n-1
  !>
  !> Output
  !>     x        : solution
  !>     a, b, c  : destroyed

  pure subroutine CyclicTridiagonalSolver(x, a, b, c, regular)
    real(RNP), intent(inout) :: x(:)    !< f on input, x on output
    real(RNP), intent(inout) :: a(:)    !< coefficients a(i), will be destroyed!
    real(RNP), intent(inout) :: b(:)    !< coefficients b(i), will be destroyed!
    real(RNP), intent(inout) :: c(:)    !< coefficients c(i), will be destroyed!
    logical,   intent(in)    :: regular !< true/false for regular/irregular case

    real(RNP) :: h, d, e, p, q, z
    integer   :: i, n

    n = size(x)

    ! factorization ............................................................

    if (abs(b(1)) < epsilon(b)) return
    h = 1 / b(1)
    d = h*x(1)
    e =-h*a(1)
    p =-h*c(1)
    a(1) = d
    b(1) = p
    c(1) = e
    do i = 2, n-1
       q = b(i) + a(i)*b(i-1)
       if (abs(q) < epsilon(q)) then
         error stop 'CyclicTridiagonalSolver: failed'
       end if
       h    = 1 / q
       b(i) =-h*c(i)
       c(i) =-h*a(i)*c(i-1)
       a(i) = h*(x(i) - a(i)*a(i-1))
       d    = p*a(i) + d
       e    = p*c(i) + e
       p    = p*b(i)
    end do
    e = e + p
    h = b(n-1) + c(n-1)

    ! solution .................................................................

    if (regular) then
       q = b(n) + a(n)*h + c(n)*e
       if (abs(q) < epsilon(q)) return
       z = (x(n) - a(n)*a(n-1) - c(n)*d) / q
    else
       z = 0
    end if
    x(n) = z
    x(n-1) = a(n-1) + h*z
    do i = n-2, 1, -1
       x(i) = a(i) + b(i)*x(i+1) + c(i)*z
    end do

  end subroutine CyclicTridiagonalSolver

  !-----------------------------------------------------------------------------
  !> LU decomposition of a regular matrix with pivoting.
  !>
  !> Returns the LU decomposition C of the given matrix A with respect to a
  !> suitable permutation which is recorded in p. This routine is used in
  !> combination with LU_Solver to solve linear equations or invert a matrix.
  !>
  !> Note that for minimization of round-off errors a higher precision is used
  !> for C and related arithmetic operations.

  pure subroutine LU_Decomposition_Pivot(A, C, p)
    real(RNP), intent(in)  :: A(:,:)  !< regular matrix, A(n,n)
    real(RHP), intent(out) :: C(:,:)  !< LU decompostion of A, C(n,n)
    integer,   intent(out) :: p(:)    !< permutation vector, p(n)

    real(RHP), dimension(size(p)) :: vv, swap
    integer :: i, n, imax

    C  = A
    n  = size(p)
    vv = maxval(abs(C),dim=2)

    if (any(vv == 0)) then
      error stop 'LU_Decomposition_Pivot: failed'
    end if

    vv = 1/vv
    do i = 1, n
      imax = transfer(maxloc(vv(i:n)*abs(C(i:n,i))),imax) + i - 1
      if (i /= imax) then
         swap = C(imax,:);  C(imax,:) = C(i,:);  C(i,:) = swap
         vv(imax) = vv(i)
      end if
      p(i) = imax
      if (C(i,i) == 0) C(i,i) = 1000 * tiny(C)
      C(i+1:n,i) = C(i+1:n,i) / C(i,i)
      C(i+1:n,i+1:n) = C(i+1:n,i+1:n) &
                     - spread(C(i+1:n,i),dim=2,ncopies=n-i) &
                     * spread(C(i,i+1:n),dim=1,ncopies=n-i)
    end do

  end subroutine LU_Decomposition_Pivot

  !-----------------------------------------------------------------------------
  !> LU decomposition of a regular matrix without pivoting.
  !>
  !> Returns the LU decomposition C of the given matrix A .
  !>
  !> Note that for minimization of round-off errors a higher precision is used
  !> for C and related arithmetic operations.

  pure subroutine LU_Decomposition_NoPivot(A, C)
    real(RNP), intent(in)  :: A(:,:)  !< regular matrix, A(n,n)
    real(RHP), intent(out) :: C(:,:)  !< LU decompostion of A, C(n,n)

    real(RHP), dimension(size(A,1)) :: vv
    integer :: i, n

    C  = A
    n  = size(vv)
    vv = maxval(abs(C),dim=2)

    if (any(vv == 0)) then
      error stop 'LU_Decomposition_NoPivot: failed'
    end if

    vv = 1/vv
    do i = 1, n
      if (C(i,i) == 0) C(i,i) = 1000 * tiny(C)
      C(i+1:n,i) = C(i+1:n,i) / C(i,i)
      C(i+1:n,i+1:n) = C(i+1:n,i+1:n) &
                     - spread(C(i+1:n,i),dim=2,ncopies=n-i) &
                     * spread(C(i,i+1:n),dim=1,ncopies=n-i)
    end do

  end subroutine LU_Decomposition_NoPivot

  !-----------------------------------------------------------------------------
  !> Solves a the linear system, AX = B using the LU decomposition of A.
  !>
  !> Solves a set of n linear equations, AX = B. The input matrix C is the LU
  !> decomposition of the nxn coefficient matrix A provided by LU_Decomposition.
  !> P is the permutation vector returned by that routine. b is input as the
  !> right-hand-side vector, and returns with the solution vector X.
  !> This routine takes into acount the possibility that b will begin with many
  !> zero elements, so it is efficient for use in matrix inversion.
  !>
  !> Note that for minimization of round-off errors a higher precision is used
  !> for C and internal arithmetic operations.

  pure subroutine LU_Solver1(C, p, B)
    real(RHP), intent(in)    :: C(:,:) !< LU decompostion of A, C(n,n)
    integer,   intent(in)    :: p(:)   !< permutation vector, p(n)
    real(RNP), intent(inout) :: B(:)   !< RHS (in), solution (out)

    real(RHP) :: X(size(B)), s
    integer   :: i, n, ii, ll

    n = size(p)
    X = B
    ii = 0
    do i = 1, n
       ll = p(i)
       s  = X(ll)
       X(ll) = X(i)
       if (ii /= 0) then
          s = s - dot_product(C(i,ii:i-1), X(ii:i-1))
       else if (s /= 0) then
          ii = i
       end if
       X(i) = s
    end do
    do i = n, 1, -1
       X(i) = (X(i) - dot_product(C(i,i+1:n), X(i+1:n)))/C(i,i)
    end do
    B = real(X, kind=RNP)

  end subroutine LU_Solver1

  !-----------------------------------------------------------------------------
  !> Cholesky factorization of a real symmetric definite matrix
  !>
  !> In-place Cholesky factorization of a definite matrix using DPOTRF from
  !> LAPACK.
  !>
  !> @note
  !> Thread-safe, can be used with OpenMP.
  !> @endnote

  subroutine CholeskyDecomposition(A)
    real(RDP), intent(inout) :: A(:,:)  !< symmetric definite matrix A

    integer :: n, info

    n = size(A,1)
    call DPOTRF('L', n, A, n, info)

    if (info /= 0) then
      error stop 'CholeskyDecomposition: DPOTRF failed'
    end if

  end subroutine CholeskyDecomposition

  !---------------------------------------------------------------------------
  !> Solves AX = B for a single RHS using the Cholesky decomcosition of A
  !>
  !> Solves the linear system AX = B for the real symmetric definite matrix A
  !> using the Cholesky decomposition computed with CholeskyDecomposition.
  !>
  !> @note
  !> Thread-safe, can be used with OpenMP.
  !> @endnote

  subroutine CholeskySolver1(A, B)
    real(RDP), intent(in)    :: A(:,:)  !< Cholesky factor L
    real(RDP), intent(inout) :: B(:)    !< RHS (in), solution (out)

    integer :: n, info

    n = size(A,1)
    call DPOTRS('L', n, 1, A, n, b, n, info)

    if (info /= 0) then
      error stop 'CholeskySolver1: DPOTRS failed'
    end if

  end subroutine CholeskySolver1

  !---------------------------------------------------------------------------
  !> Solves AX = B for multiple RHS using the Cholesky decomcosition of A
  !>
  !> Solves the linear system AX = B for the real symmetric definite matrix A
  !> using the Cholesky decomposition computed with CholeskyDecomposition.
  !>
  !> @note
  !> Thread-safe, can be used with OpenMP.
  !> @endnote

  subroutine CholeskySolver2(A, B)
    real(RDP), intent(in)    :: A(:,:)  !< Cholesky factor L
    real(RDP), intent(inout) :: B(:,:)  !< RHS (in), solution (out)

    integer :: n, info

    n = size(A,1)
    call DPOTRS('L', n, size(B,2), A, n, b, n, info)

    if (info /= 0) then
      error stop 'CholeskySolver2: DPOTRS failed'
    end if

  end subroutine CholeskySolver2

  !============================================================================

end module Linear_Equations
