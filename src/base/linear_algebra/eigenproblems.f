!> summary:   Solvers for eigenproblems
!> author:    Joerg Stiller, Immo Huismann
!> date:      2013/06/26
!> license:   Institute of Fluid Mechanics &
!>
!>### Solvers for eigenproblems
!>
!> @note
!>   * Choosing the parameter `abstol` too small in SolveSymmetricEigenproblem
!>     can trigger a divide by zero in DSYEVR, which may remain undetected in
!>     Release buildi and cause a lot of trouble when debugging.
!>
!==============================================================================

module Eigenproblems
  use Kind_Parameters, only: RDP
  implicit none
  private

  public :: SolveSymmetricEigenproblem
  public :: SolveNonsymmetricEigenproblem
  public :: SolveGeneralizedEigenproblem

contains

!-------------------------------------------------------------------------------
!> Solves the eigenproblem for a real symmetric matrix
!>
!> The real symmetric matrix A(n,n) possesses n real eigenvalues {λᵢ} that are
!> sorted in ascending order. Given the range i = [i₁,i₂] the routine computes
!> the eigenvalues λᵢ₁ to λᵢ₂. If omitted, i₁ = 1 and i₂ = n are assumed.
!> The computed eigenvalues are returned as
!>
!>     lambda(1:m) = [(λᵢ, i=i₁,i₂)],   m = i₂ - i₁ + 1
!>
!> If V is present, V(:,m) contains the normalized eigenvector to lambda(m).
!>
!> @note
!> Thread-safe, can be used with OpenMP.
!> @endnote

subroutine SolveSymmetricEigenproblem(A, lambda, V, i, success)
  real(RDP), intent(in)  :: A(:,:)    !< real symmetric matrix, A(n,n)
  real(RDP), intent(out) :: lambda(:) !< computed eigenvalues in ascending order
  real(RDP), intent(out) :: V(:,:)    !< orthonormal eigenvectors to lambda
  integer,   intent(in)  :: i(2)      !< range of eigenvalues to be computed
  logical,   intent(out) :: success   !< true upon successful execution
  optional :: V, i, success

  ! interface to LAPACK and auxiliary parameters
  character :: jobv
  real(RDP) :: A_(size(A,1),size(A,2))
  real(RDP) :: lambda_(size(A,1))
  real(RDP) :: V_(size(A,1),size(A,1))
  real(RDP) :: abstol, vl, vu
  real(RDP) :: work(100*size(A,1))
  integer   :: iwork(10*size(A,1)), isuppz(2*size(A,1))
  integer   :: il, iu, m, n, info

!#ifdef CHECK
!  print '(A)', '*** SolveSymmetricEigenproblem *** #0'
!  print '(A,2(1X,G0))', '*** shape(A)  =', shape(A)
!  print '(A,2(1X,G0))', '*** max|A-Aᵀ| =', maxval(abs(A - transpose(A)))
!  print '(A,2(1X,G0))', '*** shape(lambda) =', shape(lambda)
!  if (present(V)) then
!    print '(A,2(1X,G0))', '*** shape(V) =', shape(V)
!  end if
!#endif

  ! intialization ..............................................................

  n  = size(A,1)
  A_ = A

  if (present(i)) then
     il = max(1, min(i(1), n))
     iu = min(n, max(i(2), i(1)))
  else
     il = 1
     iu = n
  end if

  if (present(V)) then
     jobv = 'V'
  else
     jobv = 'N'
  end if

  vl = 0
  vu = 1
  abstol = 1000 * epsilon(1D0)

  ! solve eigenproblem .........................................................

  call DSYEVR(jobv,'I','U', n, A_, n, vl, vu, il, iu, abstol, m, lambda_, V_, &
              n, isuppz, work, size(work), iwork, size(iwork), info)

  ! assign results .............................................................

  lambda = lambda_

  if (present(V)) then
     V(:,1:m) = V_(:,1:m)
  end if

  if (present(success)) success = info == 0

end subroutine SolveSymmetricEigenproblem

!-------------------------------------------------------------------------------
!> Solves the eigenproblem for a nonsymmetric matrix
!>
!> Computes the eigenvalues λᵢ of a given real matrix A(n,n) and, optionally,
!> the left and right eigenvectors vlᵢ and vrᵢ.                             <br>
!>
!> Eigenvalues
!>
!>   *  Are stored in the complex array lambda.
!>   *  Complex conjugate pairs of eigenvalues appear consecutively with
!>      pairs of eigenvalues appear consecutively with the eigenvalue
!>      having positive imaginary part first.
!>
!> Eigenvectors
!>
!>   * If the eigenvalue lambda(m) is real, then Vl(:,m) and Vr(:,m) contain
!>     the corresponding left and right eigenvectors.
!>   * If lambda(m:m+1) form a complex conjugate pair, the corresponding
!>     eigenvectors are Vl(:,m) ± i*Vl(:,m+1) and Vr(:,m) ± i*Vr(:,m+1),
!>     respectively.
!>   * The computed eigenvectors are normalized to have Euclidean norm
!>     equal to 1 and largest component real.
!>
!> @note
!> Thread-safe, can be used with OpenMP.
!> @endnote

subroutine SolveNonsymmetricEigenproblem(A, lambda, Vl, Vr, success)
  real(RDP),    intent(in)  :: A(:,:)
  complex(RDP), intent(out) :: lambda(:)
  real(RDP),    intent(out) :: Vl(:,:)
  real(RDP),    intent(out) :: Vr(:,:)
  logical,      intent(out) :: success
  optional :: Vl, Vr, success

  ! interface to LAPACK and auxiliary parameters
  real(RDP), dimension(size(A,1),size(A,2)) :: A_, Vl_, Vr_
  real(RDP), dimension(size(lambda)) :: wr, wi
  real(RDP) :: work(100*size(A,1))
  character :: jobvl, jobvr
  integer   :: n, info

  ! intialization ..............................................................

  n  = size(A,1)
  A_ = A

  if (present(Vl)) then
     jobvl = 'V'
  else
     jobvl = 'N'
  end if
  if (present(Vr)) then
     jobvr = 'V'
  else
     jobvr = 'N'
  end if

  ! solve eigenproblem .........................................................

  call DGEEV(jobvl, jobvr, n, A_, n, wr, wi, Vl_, n, Vr_, n, work, size(work), &
             info)

  ! assign results .............................................................

  lambda = cmplx(wr, wi, kind=RDP)

  if (present(Vl)) Vl = Vl_
  if (present(Vr)) Vr = Vr_

  if (present(success)) success = info == 0

end subroutine SolveNonsymmetricEigenproblem

!-------------------------------------------------------------------------------
!> Solves a generalized eigenproblem for matrix A
!>
!> Solves the generalized eigenproblem
!>
!>    A V = lambda D V
!>
!> where
!>
!>   * A is a symmetric definite matrix,
!>   * D a positive diagonal matrix,
!>   * lambda the diagonal matrix of eigenvalues, and
!>   * V the matrix comprising the eigen vectors.
!>
!> The eigensytem satisfies the conditions:
!>
!>   *  V^T D V  is the identity matrix, and
!>   *  V^T A V  is diagonal.
!>
!> @note
!> Thread-safe, can be used with OpenMP.
!> @endnote

subroutine SolveGeneralizedEigenproblem(A, d, lambda, V, success)
  real(RDP), intent(in)  :: A(:,:)    !< symmetric definite matrix
  real(RDP), intent(in)  :: d(:)      !< diagonal of d
  real(RDP), intent(out) :: lambda(:) !< generalized eigenvalues
  real(RDP), intent(out) :: V(:,:)    !< generalized eigenvectors in columns
  logical,   intent(out) :: success   !< true upon successful execution
  optional :: success

  real(RDP) :: C(size(A,1),size(A,2))
  integer   :: i, j, n
  logical   :: passed

  n = size(d)

  ! scale the matrix A to get an ordinary eigenproblem
  forall (i = 1:n, j = 1:n)
    C(i,j) = A(i,j) / (sqrt(d(i)) * sqrt(d(j)))
  end forall

  ! solve the eigenproblem for the scaled matrix
  call SolveSymmetricEigenproblem(C, lambda, V, success=passed)

  ! scale the eigenvectors to ensure V^T D V = I
  forall (i = 1:n, j = 1:n)
    V(i,j) = V(i,j) / sqrt(d(i))
  end forall

  if (present(success)) success = passed

end subroutine SolveGeneralizedEigenproblem

!==============================================================================

end module Eigenproblems
