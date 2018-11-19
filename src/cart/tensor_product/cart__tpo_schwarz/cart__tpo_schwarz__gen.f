!> summary:  Schwarz operator: bf-marching, separate, generic
!> author:   Joerg Stiller
!> date:     2018/11/17
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Schwarz operator: bf-marching, separate, generic
!===============================================================================

subroutine CART__TPO_Schwarz__gen( n1, n2, n3, nc, nd, S1, S2, S3, &
                                   W1, W2, W3, cfg, D_inv, f, u    )

  !-----------------------------------------------------------------------------
  ! modules

  use Kind_Parameters, only: RNP
  implicit none

  !-----------------------------------------------------------------------------
  ! arguments

  integer,   intent(in)  :: n1                 !< dimension of S1
  integer,   intent(in)  :: n2                 !< dimension of S2
  integer,   intent(in)  :: n3                 !< dimension of S3
  integer,   intent(in)  :: nc                 !< number of configurations
  integer,   intent(in)  :: nd                 !< number of subdomains
  real(RNP), intent(in)  :: S1(n1,n1,nc)       !< 1D eigenvectors for dir 1
  real(RNP), intent(in)  :: S2(n2,n2,nc)       !< 1D eigenvectors for dir 2
  real(RNP), intent(in)  :: S3(n3,n3,nc)       !< 1D eigenvectors for dir 3
  real(RNP), intent(in)  :: W1(n1,nc)          !< 1D weights for dir 1
  real(RNP), intent(in)  :: W2(n2,nc)          !< 1D weights for dir 2
  real(RNP), intent(in)  :: W3(n3,nc)          !< 1D weights for dir 3
  integer,   intent(in)  :: cfg(3,nd)          !< subdomain configurations
  real(RNP), intent(in)  :: D_inv(n1,n2,n3,nd) !< inverse 3D eigenvalues
  real(RNP), intent(in)  :: f(n1,n2,n3,nd)     !< RHS
  real(RNP), intent(out) :: u(n1,n2,n3,nd)     !< solution

  !-----------------------------------------------------------------------------
  ! local variables

  real(RNP), dimension(:,:,:), allocatable :: WS1_t, WS2_t, WS3_t, y, z
  real(RNP) :: tmp

  integer :: i, j, k, l, m
  integer :: c, c1, c2, c3
  integer :: vec_len

  !-----------------------------------------------------------------------------
  ! initialization

  ! OpenACC vector length
  if (min(n1,n2,n3) < 8) then
    vec_len = 128
  else
    vec_len = 256
  end if

  allocate(WS1_t(n1,n1,nc))
  do c = 1, nc
  do j = 1, n1
  do i = 1, n1
    WS1_t(j,i,c) = W1(i,c) * S1(i,j,c)
  end do
  end do
  end do

  allocate(WS2_t(n2,n2,nc))
  do c = 1, nc
  do j = 1, n2
  do i = 1, n2
    WS2_t(j,i,c) = W2(i,c) * S2(i,j,c)
  end do
  end do
  end do

  allocate(WS3_t(n3,n3,nc))
  do c = 1, nc
  do j = 1, n3
  do i = 1, n3
    WS3_t(j,i,c) = W3(i,c) * S3(i,j,c)
  end do
  end do
  end do

  allocate(y(n1, n2, n3), z(n1, n2, n3))

  !-----------------------------------------------------------------------------
  ! evaluation

  !$acc data async &
  !$acc & present(cfg, D_inv, f, u) &
  !$acc & copyin(S1, S2, S3, WS1_t, WS2_t, WS3_t)

  !$acc parallel async &
  !$acc & device_type(nvidia) num_workers(1024/vec_len) vector_length(vec_len)

  !$acc loop gang worker private(y,z)

  !$omp do private(l)
  subdomains: do l = 1, nd

    ! element-boundary configurations
    c1 = cfg(1,l)
    c2 = cfg(2,l)
    c3 = cfg(3,l)

    ! y = S₃ᵀ x I x I f ........................................................

    !$acc loop collapse(3) vector
    do k = 1, n3
    do j = 1, n2
    do i = 1, n1
      tmp = 0
      do m = 1, n3
        tmp = tmp + S3(m,k,c3) * f(i,j,m,l)
      end do
      y(i,j,k) = tmp
    end do
    end do
    end do

    ! z = I x S₂ᵀ x I y ........................................................

    !$acc loop vector collapse(3)
    do k = 1, n3
    do j = 1, n2
    do i = 1, n1
      tmp = 0
      do m = 1, n2
        tmp = tmp + S2(m,j,c2) * y(i,m,k)
      end do
      z(i,j,k) = tmp
    end do
    end do
    end do

    ! y = D⁻¹ (I x I x S₁ᵀ) z ..................................................

    !$acc loop vector collapse(3)
    do k = 1, n3
    do j = 1, n2
    do i = 1, n1
      tmp = 0
      do m = 1, n1
        tmp = tmp + S1(m,i,c1) * z(m,j,k)
      end do
      y(i,j,k) = tmp * D_inv(i,j,k,l)
    end do
    end do
    end do

    ! z = I x I x WS₁ y ........................................................

    !$acc loop vector collapse(3)
    do k = 1, n3
    do j = 1, n2
    do i = 1, n1
      tmp = 0
      do m = 1, n1
        tmp = tmp + WS1_t(m,i,c1) * y(m,j,k)
      end do
      z(i,j,k) = tmp
    end do
    end do
    end do

    ! y = I x WS₂ x I z ........................................................

    !$acc loop vector collapse(3)
    do k = 1, n3
    do j = 1, n2
    do i = 1, n1
      tmp = 0
      do m = 1, n2
        tmp = tmp + WS2_t(m,j,c2) * z(i,m,k)
      end do
      y(i,j,k) = tmp
    end do
    end do
    end do

    ! u = WS₃ x I x I y ........................................................

    !$acc loop collapse(3) vector
    do k = 1, n3
    do j = 1, n2
    do i = 1, n1
      tmp = 0
      do m = 1, n3
        tmp = tmp + WS3_t(m,k,c3) * y(i,j,m)
      end do
      u(i,j,k,l) = tmp
    end do
    end do
    end do

  end do subdomains
  !$omp end do

  !$acc end parallel
  !$acc end data

!===============================================================================

end subroutine  CART__TPO_Schwarz__gen
