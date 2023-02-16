!> summary:  3d generic isotropic Schwarz operator
!> author:   Joerg Stiller
!> date:     2020/06/02
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

subroutine TPO_Schwarz_I_Gen_RWP(S, V, W, g, cfg, lambda, nu, f, u)
  real(RWP), intent(in)  :: S(:,:,:)       !< 1D eigenvectors   (np,np,nc)
  real(RWP), intent(in)  :: V(:,:)         !< 1D eigenvalues    (np,nc)
  real(RWP), intent(in)  :: W(:,:)         !< 1D weights        (np,nc)
  real(RWP), intent(in)  :: g(:,:)         !< metrics           (4,nd)
  integer,   intent(in)  :: cfg(:,:)       !< subdomain configs (3,nd)
  real(RWP), intent(in)  :: lambda         !< λ
  real(RWP), intent(in)  :: nu(:)          !< ν + νˢ            (nd)
  real(RWP), intent(in)  :: f(:,:,:,:)     !< RHS               (np,np,np,nd)
  real(RWP), intent(out) :: u(:,:,:,:)     !< solution          (np,np,np,nd)

  real(RWP), parameter :: eps = epsilon(lambda)

  !-----------------------------------------------------------------------------
  ! local variables

  real(RWP) :: WS_t ( size(S,1), size(S,2), size(S,3) )
  real(RWP) :: y    ( size(f,1), size(f,2), size(f,3) )
  real(RWP) :: z    ( size(f,1), size(f,2), size(f,3) )
  real(RWP) :: a1, a2, a3, a4, d, tmp

  integer :: nc, nd, np
  integer :: i, j, k, l, p
  integer :: c, c1, c2, c3

  !-----------------------------------------------------------------------------
  ! initialization

  np = size(S,1)
  nc = size(S,3)
  nd = size(f,4)

  ! WS_t(:,:,c) = (WS)ᵀ(:,:,c)
  do c = 1, nc
  do j = 1, np
  do i = 1, np
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

    ! coefficients
    a1 = g(1,l) * nu(l)
    a2 = g(2,l) * nu(l)
    a3 = g(3,l) * nu(l)
    a4 = g(4,l) * lambda

    ! y = Sᵀ x I x I f .......................................................

    !$acc loop collapse(3) vector
    do k = 1, np
    do j = 1, np
      !DIR$ SIMD
      do i = 1, np
        tmp = 0
        do p = 1, np
          tmp = tmp + S(p,k,c3) * f(i,j,p,l)
        end do
        y(i,j,k) = tmp
      end do
    end do
    end do

    ! z = I x Sᵀ x I y .......................................................

    !$acc loop collapse(3) vector
    do k = 1, np
    do j = 1, np
      !DIR$ SIMD
      do i = 1, np
        tmp = 0
        do p = 1, np
          tmp = tmp + S(p,j,c2) * y(i,p,k)
        end do
        z(i,j,k) = tmp
      end do
    end do
    end do

    ! y = D⁻¹ (I x I x Sᵀ) z .................................................

    !$acc loop collapse(3) vector
    do k = 1, np
    do j = 1, np
      !DIR$ SIMD
      do i = 1, np
        tmp = 0
        do p = 1, np
          tmp = tmp + S(p,i,c1) * z(p,j,k)
        end do
        d = a1 * V(i,c1) + a2 * V(j,c2) + a3 * V(k,c3) + a4
        y(i,j,k) = tmp / sign(max(abs(d),eps), d)
      end do
    end do
    end do

    ! z = I x I x WS y .......................................................

    !$acc loop collapse(3) vector
    do k = 1, np
    do j = 1, np
      !DIR$ SIMD
      do i = 1, np
        tmp = 0
        do p = 1, np
          tmp = tmp + WS_t(p,i,c1) * y(p,j,k)
        end do
        z(i,j,k) = tmp
      end do
    end do
    end do

    ! y = I x WS x I z .......................................................

    !$acc loop collapse(3) vector
    do k = 1, np
    do j = 1, np
      !DIR$ SIMD
      do i = 1, np
        tmp = 0
        do p = 1, np
          tmp = tmp + WS_t(p,j,c2) * z(i,p,k)
        end do
        y(i,j,k) = tmp
      end do
    end do
    end do

    ! u = WS x I x I y .......................................................

    !$acc loop collapse(3) vector
    do k = 1, np
    do j = 1, np
      !DIR$ SIMD
      do i = 1, np
        tmp = 0
        do p = 1, np
          tmp = tmp + WS_t(p,k,c3) * y(i,j,p)
        end do
        u(i,j,k,l) = tmp
      end do
    end do
    end do

  end do subdomains

  !$acc end parallel
  !$acc end data

end subroutine TPO_Schwarz_I_Gen_RWP
