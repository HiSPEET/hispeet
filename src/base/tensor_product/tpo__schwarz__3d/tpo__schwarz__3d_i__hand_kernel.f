!-------------------------------------------------------------------------------
!> Parametrized 3d isotropic Schwarz kernel using hand-crafted suboperators

subroutine PROC(TPO_Schwarz_I_Hand__,_NP_)(nc, nd, S, V, W, g, cfg, lambda, nu, f, u)
  integer,   intent(in)  :: nc                       !< num configurations
  integer,   intent(in)  :: nd                       !< num subdomains
  real(RWP), intent(in)  :: S(_NP_,_NP_,nc)          !< 1D eigenvectors
  real(RWP), intent(in)  :: V(_NP_,nc)               !< 1D eigenvalues
  real(RWP), intent(in)  :: W(_NP_,nc)               !< 1D weights
  real(RWP), intent(in)  :: g(4,nd)                  !< metrics
  integer,   intent(in)  :: cfg(3,nd)                !< subdomain configurations
  real(RWP), intent(in)  :: lambda                   !< λ
  real(RWP), intent(in)  :: nu(nd)                   !< ν + νˢ
  real(RWP), intent(in)  :: f(_NP_,_NP_,_NP_,nd)     !< RHS
  real(RWP), intent(out) :: u(_NP_,_NP_,_NP_,nd)     !< solution

  real(RWP), parameter :: alpha = 1
  real(RWP), parameter :: beta  = 0
  real(RWP), parameter :: eps   = epsilon(lambda)

  real(RWP) :: WS_t ( _NP_, _NP_,  nc  )
  real(RWP) :: y    ( _NP_, _NP_, _NP_ )
  real(RWP) :: z    ( _NP_, _NP_, _NP_ )
  real(RWP) :: a1, a2, a3, a4, d

  integer :: i, j, k, l
  integer :: c, c1, c2, c3

  !---------------------------------------------------------------------------
  ! initialization

  ! WS_t(:,:,c) = (WS)ᵀ(:,:,c)
  do c = 1, nc
  do j = 1, _NP_
  do i = 1, _NP_
    WS_t(j,i,c) = W(i,c) * S(i,j,c)
  end do
  end do
  end do

  ! assure that aux arrays do not contain NaNs
  y = 0
  z = 0

  !---------------------------------------------------------------------------
  ! evaluation

  !$acc data async present(cfg, D_inv, f, u) copyin(S, WS_t)
  !$acc parallel
  !$acc loop gang worker private(y,z)

  !$omp do private(l)
  subdomains: do l = 1, nd

    ! element-boundary configurations
    c1 = cfg(1,l)
    c2 = cfg(2,l)
    c3 = cfg(3,l)

    ! coefficients
    a1 = g(1,l) * nu(l)
    a2 = g(2,l) * nu(l)
    a3 = g(3,l) * nu(l)
    a4 = g(4,l) * lambda

    ! y = Sᵀ x I x I uᵉ
    call PROC(QtxIxI__,_NP_)(S(:,:,c3), alpha, beta, f(:,:,:,l), y)

    ! z = I x Sᵀ x I y
    call PROC(IxQtxI__,_NP_)(S(:,:,c2), alpha, beta, y, z)

    ! y = I x I x Sᵀ z
    call PROC(IxIxQt__,_NP_)(S(:,:,c1), alpha, beta, z, y)

    ! y = D⁻¹ y
    do k = 1, _NP_
    do j = 1, _NP_
    do i = 1, _NP_
      d = a1 * V(i,c1) + a2 * V(j,c2) + a3 * V(k,c3) + a4
      y(i,j,k) = y(i,j,k) / sign(max(abs(d),eps), d)
    end do
    end do
    end do

    ! z = I x I x WS y
    call PROC(IxIxQt__,_NP_)(WS_t(:,:,c1), alpha, beta, y, z)

    ! y  = I x WS x I z
    call PROC(IxQtxI__,_NP_)(WS_t(:,:,c2), alpha, beta, z, y)

    ! uᵉ = WS x I x I y
    u(:,:,:,l) = 0  ! avoid trouble with NaNs
    call PROC(QtxIxI__,_NP_)(WS_t(:,:,c3), alpha, beta, y, u(:,:,:,l))

  end do subdomains

  !$acc end parallel
  !$acc end data

end subroutine PROC(TPO_Schwarz_I_Hand__,_NP_)
