











program Vec_Test
  implicit none
  integer, parameter :: RNP = kind(1D0)
  integer   :: nb, nc, ne
  integer   :: j, k, e
  real(RNP) :: A(8,16)
  real(RNP) :: alpha = 1
  real(RNP) :: beta  = 1
  real(RNP), allocatable :: u(:,:,:,:), v(:,:,:,:), w(:,:,:,:)
  real(RNP) :: error, nope, mflops, t0, t

  nb = 8
  nc = 8
  ne = 2000

! FLOPs per element
  nope = (2*8 + 3) * 16 * nb * nc

  allocate(u(8,nb,nc,ne))
  allocate(v(16,nb,nc,ne))
  allocate(w, mold=v)

  call random_number(A)
  call random_number(u)
  call random_number(v)

! reference
  do e = 1, ne
    do k = 1, nc
    do j = 1, nb
      w(:,j,k,e) = alpha * matmul(u(:,j,k,e), A) + beta * v(:,j,k,e)
    end do
    end do
  end do

! test
  call cpu_time(t0)
  do e = 1, ne
    call IxIxAt__8x16(nb, nc, A, alpha, beta, u(:,:,:,e), v(:,:,:,e))
  end do
  call cpu_time(t)

! report
  t       =  t - t0
  error   =  maxval(abs(v - w))
  mflops  =  1e-6 * nope * ne / t
  print '(A,ES12.5)', 'cpu time = ', t
  print '(A,ES12.5)', 'MFLOPs   = ', mflops
  print '(A,ES12.5)', 'error    = ', error

contains

!Computes  v = α I⊗I⊗Aᵀ u + β v
!>
!>   * k:  block 1,  unroll 1
!>   * j:  block 1,  unroll 1
!>   * i:  block 1,  unroll 2
!>   * p:  block 1,  unroll 1
!>   * using Intel SIMD directive

subroutine IxIxAt__8x16(nb, nc, A, alpha, beta, u, v)
!$acc routine vector
  integer,   intent(in)    :: nb              !< 2nd dimension of u,v
  integer,   intent(in)    :: nc              !< 3rd dimension of u,v
  real(RNP), intent(in)    :: A(8,16)  !< rectangular matrix
  real(RNP), intent(in)    :: alpha           !< factor α
  real(RNP), intent(in)    :: beta            !< factor β
  real(RNP), intent(in)    :: u(8,nb,nc)  !< operand
  real(RNP), intent(inout) :: v(16,nb,nc)  !< result

  real(RNP) :: tmp0, tmp1
  integer   :: i, j, k, p



  do k = 1, nc
  do j = 1, nb
!DIR$ SIMD
    do i = 1, (16 / 2) * 2, 2
     tmp0 = 0
     tmp1 = 0
      do p = 1, 8
        tmp0 = tmp0 + A(p,i  ) * u(p,j,k)
        tmp1 = tmp1 + A(p,i+1) * u(p,j,k)
      end do
      v(i  ,j,k) = alpha * tmp0 + beta * v(i  ,j,k)
      v(i+1,j,k) = alpha * tmp1 + beta * v(i+1,j,k)
    end do
  end do
  end do





end subroutine IxIxAt__8x16

end program Vec_Test
