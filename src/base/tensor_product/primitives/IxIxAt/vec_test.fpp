
# ifdef __GFORTRAN__
#   define PASTE(base_name) base_name
#   define PROC(base_name,na1,na2) PASTE(base_name)PASTE(na1)PASTE(x)na2
# else
#  define PASTE(base_name,na1,na2) base_name ## na1 ## x ##na2
#  define PROC(base_name,na1,na2) PASTE(base_name,na1,na2)
# endif

# define _NA1_ 8
# define _NA2_ 16

program Vec_Test
  implicit none
  integer, parameter :: RNP = kind(1D0)
  integer   :: nb, nc, ne
  integer   :: j, k, e
  real(RNP) :: A(_NA1_,_NA2_)
  real(RNP) :: alpha = 1
  real(RNP) :: beta  = 1
  real(RNP), allocatable :: u(:,:,:,:), v(:,:,:,:), w(:,:,:,:)
  real(RNP) :: error, nope, mflops, t0, t

  nb = 8
  nc = 8
  ne = 2000

  ! FLOPs per element
  nope = (2*_NA1_ + 3) * _NA2_ * nb * nc

  allocate(u(_NA1_,nb,nc,ne))
  allocate(v(_NA2_,nb,nc,ne))
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
    call PROC(IxIxAt__,_NA1_,_NA2_)(nb, nc, A, alpha, beta, u(:,:,:,e), v(:,:,:,e))
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

# include "IxIxAt_01_01_01_01.f"

end program Vec_Test
