!-------------------------------------------------------------------------------
!> Parametrized 3d isotropic Schwarz kernel using hand-crafted suboperators

subroutine PROC(TPO_Schwarz_CI_Hand__,_NP_)(nc, nd, S, W, cfg, D_inv, f, u)
  integer,   intent(in)  :: nc                       !< num configurations
  integer,   intent(in)  :: nd                       !< num subdomains
  real(RNP), intent(in)  :: S(_NP_,_NP_,nc)          !< 1D eigenvectors
  real(RNP), intent(in)  :: W(_NP_,nc)               !< 1D weights
  integer,   intent(in)  :: cfg(:,:)                 !< subdomain configurations
  real(RNP), intent(in)  :: D_inv(_NP_,_NP_,_NP_,nd) !< inverse 3D eigenvalues
  real(RNP), intent(in)  :: f(_NP_,_NP_,_NP_,nd)     !< RHS
  real(RNP), intent(out) :: u(_NP_,_NP_,_NP_,nd)     !< solution

  real(RNP), parameter :: alpha = 1
  real(RNP), parameter :: beta  = 0

  real(RNP) :: WS_t ( _NP_, _NP_,  nc  )
  real(RNP) :: y    ( _NP_, _NP_, _NP_ )
  real(RNP) :: z    ( _NP_, _NP_, _NP_ )

  integer :: i, j, l
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

    ! y = Sᵀ x I x I uᵉ
    call PROC(QtxIxI__,_NP_)(S(:,:,c3), alpha, beta, f(:,:,:,l), y)

    ! z = I x Sᵀ x I y
    call PROC(IxQtxI__,_NP_)(S(:,:,c2), alpha, beta, y, z)

    ! y = I x I x Sᵀ z
    call PROC(IxIxQt__,_NP_)(S(:,:,c1), alpha, beta, z, y)

    ! y = D⁻¹ y
    y = D_inv(:,:,:,l) * y

    ! z = I x I x WS y
    call PROC(IxIxQt__,_NP_)(WS_t(:,:,c1), alpha, beta, y, z)

    ! y  = I x WS x I z
    call PROC(IxQtxI__,_NP_)(WS_t(:,:,c2), alpha, beta, z, y)

    ! vᵉ = WS x I x I y
    call PROC(QtxIxI__,_NP_)(WS_t(:,:,c3), alpha, beta, y, u(:,:,:,l))

  end do subdomains

  !$acc end parallel
  !$acc end data

end subroutine PROC(TPO_Schwarz_CI_Hand__,_NP_)
