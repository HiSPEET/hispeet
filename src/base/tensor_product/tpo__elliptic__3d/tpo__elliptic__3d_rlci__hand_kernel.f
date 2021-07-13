!-------------------------------------------------------------------------------
!> Parametrized 3d elliptic kernel using hand-crafted suboperators (RLCI)

subroutine PROC(TPO_Elliptic_RLCI_Hand__,_NP_)(ne, Ms, Ls, dx, lambda, nu, u, v)
  integer,   intent(in)  :: ne                   !< num elements
  real(RWP), intent(in)  :: Ms(_NP_)             !< standard mass matrix
  real(RWP), intent(in)  :: Ls(_NP_,_NP_)        !< standard stiffness matrix
  real(RWP), intent(in)  :: dx(3)                !< element extensions
  real(RWP), intent(in)  :: lambda               !< Helmholtz parameter λ
  real(RWP), intent(in)  :: nu                   !< diffusivity nu
  real(RWP), intent(in)  :: u(_NP_,_NP_,_NP_,ne) !< operand
  real(RWP), intent(out) :: v(_NP_,_NP_,_NP_,ne) !< result

  real(RWP), parameter :: beta = 1
  real(RWP) :: M(_NP_,_NP_,_NP_), M_u(_NP_,_NP_,_NP_), Lm(_NP_,_NP_)
  real(RWP) :: c(3), tmp

  integer, parameter :: NP_E = _NP_**3
  integer :: e, i, j, k

  !-----------------------------------------------------------------------------
  ! initialization

  ! element mass matrix
  tmp = product(dx) / 8
  do k = 1, _NP_
  do j = 1, _NP_
  do i = 1, _NP_
    M(i,j,k) = tmp * Ms(k) * Ms(j) * Ms(i)
  end do
  end do
  end do

  ! mass-weighted stiffness matrix: Lm = Ms^-1 Ls = (Ls Ms^-1)^T
  do j = 1, _NP_
  do i = 1, _NP_
    Lm(i,j) = Ls(i,j) / Ms(i)
  end do
  end do

  ! coefficients
  c = 4 * nu / dx**2

  !---------------------------------------------------------------------------
  ! evaluation

  !$acc data present(u,v) copyin(c,M,Lm)
  !$acc parallel
  !$acc loop gang worker private(M_u)

  !$omp do
  do e = 1, ne

    call SetOperands(lambda, M, u(:,:,:,e), M_u, v(:,:,:,e))
    call PROC(IxIxQt__,_NP_)(Lm, c(1), beta, M_u, v(:,:,:,e))
    call PROC(IxQtxI__,_NP_)(Lm, c(2), beta, M_u, v(:,:,:,e))
    call PROC(QtxIxI__,_NP_)(Lm, c(3), beta, M_u, v(:,:,:,e))

  end do

  !$acc end parallel
  !$acc end data

contains

  !-----------------------------------------------------------------------------
  !> Projection to element mass matrix and computation of the Helmholtz term

  subroutine SetOperands(lambda, M, u, M_u, v)
    !$acc routine vector
    real(RWP), intent(in)  :: lambda
    real(RWP), intent(in)  :: M   (NP_E)
    real(RWP), intent(in)  :: u   (NP_E)
    real(RWP), intent(out) :: M_u (NP_E)
    real(RWP), intent(out) :: v   (NP_E)

    integer :: l

    !$acc loop vector
    !DIR$ SIMD
    do l = 1, NP_E
      M_u(l) = M(l) * u(l)
      v(l) = lambda * M_u(l)
    end do

  end subroutine SetOperands

  !-----------------------------------------------------------------------------

end subroutine PROC(TPO_Elliptic_RLCI_Hand__,_NP_)
