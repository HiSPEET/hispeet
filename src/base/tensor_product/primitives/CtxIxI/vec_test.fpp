
# ifdef __GFORTRAN__
#   define PASTE(base_name) base_name
#   define PROC(base_name,nc1,nc2) PASTE(base_name)PASTE(nc1)PASTE(x)nc2
# else
#  define PASTE(base_name,nc1,nc2) base_name ## nc1 ## x ##nc2
#  define PROC(base_name,nc1,nc2) PASTE(base_name,nc1,nc2)
# endif

# define _RWP_ RNP
# define _NC1_ 8
# define _NC1_T4_ (_NC1_ / 4) * 4
# define _NC2_ 16
# define _NC2_T2_ (_NC2_ / 2) * 2
# define _NC2_T4_ (_NC2_ / 4) * 4
# define _NC2_T8_ (_NC2_ / 8) * 8

program Vec_Test
  implicit none
  integer, parameter :: RNP = kind(1D0)
  integer   :: na, nb, ne
  integer   :: i, j, e
  real(RNP) :: C(_NC1_,_NC2_)
  real(RNP) :: alpha = 1
  real(RNP) :: beta  = 1
  real(RNP), allocatable :: u(:,:,:,:), v(:,:,:,:), w(:,:,:,:)
  real(RNP) :: error, nope, mflops, t0, t

  na = 8
  nb = 8
  ne = 2000

  ! FLOPs per element
  nope = (2*_NC1_ + 3) * _NC2_ * na * nb

  allocate(u(na,nb,_NC1_,ne))
  allocate(v(na,nb,_NC2_,ne))
  allocate(w, mold=v)

  call random_number(C)
  call random_number(u)
  call random_number(v)

  ! reference
  do e = 1, ne
    do j = 1, nb
    do i = 1, na
      w(i,j,:,e) = alpha * matmul(u(i,j,:,e), C) + beta * v(i,j,:,e)
    end do
    end do
  end do

  ! test
  call cpu_time(t0)
  do e = 1, ne
    call PROC(CtxIxI__,_NC1_,_NC2_)(na*nb, C, alpha, beta, u(:,:,:,e), v(:,:,:,e))
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

# include "CtxIxI_04_jj_41.f"

end program Vec_Test
