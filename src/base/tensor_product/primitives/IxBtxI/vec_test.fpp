
# ifdef __GFORTRAN__
#   define PASTE(base_name) base_name
#   define PROC(base_name,nb1,nb2) PASTE(base_name)PASTE(nb1)PASTE(x)nb2
# else
#  define PASTE(base_name,nb1,nb2) base_name ## nb1 ## x ##nb2
#  define PROC(base_name,nb1,nb2) PASTE(base_name,nb1,nb2)
# endif

# define _RWP_ RNP
# define _NB1_ 8
# define _NB1_T4_ (_NB1_ / 4) * 4
# define _NB2_ 16
# define _NB2_T2_ (_NB2_ / 2) * 2
# define _NB2_T4_ (_NB2_ / 4) * 4
# define _NB2_T8_ (_NB2_ / 8) * 8


program Vec_Test
  implicit none
  integer, parameter :: RNP = kind(1D0)
  integer   :: na, nc, ne
  integer   :: i, k, e
  real(RNP) :: B(_NB1_,_NB2_)
  real(RNP) :: alpha = 1
  real(RNP) :: beta  = 1
  real(RNP), allocatable :: u(:,:,:,:), v(:,:,:,:), w(:,:,:,:)
  real(RNP) :: error, nope, mflops, t0, t

  na = 8
  nc = 8
  ne = 2000

  ! FLOPs per element
  nope = (2*_NB1_ + 3) * _NB2_ * na * nc

  allocate(u(na,_NB1_,nc,ne))
  allocate(v(na,_NB2_,nc,ne))
  allocate(w, mold=v)

  call random_number(B)
  call random_number(u)
  call random_number(v)

  ! reference
  do e = 1, ne
    do k = 1, nc
    do i = 1, na
      w(i,:,k,e) = alpha * matmul(u(i,:,k,e), B) + beta * v(i,:,k,e)
    end do
    end do
  end do

  ! test
  call cpu_time(t0)
  do e = 1, ne
    call PROC(IxBtxI__,_NB1_,_NB2_)(na, nc, B, alpha, beta, u(:,:,:,e), v(:,:,:,e))
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

# include "IxBtxI_01_04_01_41.f"

end program Vec_Test
