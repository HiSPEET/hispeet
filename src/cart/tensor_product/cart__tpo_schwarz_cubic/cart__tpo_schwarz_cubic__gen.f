!> summary:  Schwarz operator cubic subdomains: bf-marching, separate, generic
!> author:   Joerg Stiller
!> date:     2017/10/16
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Schwarz operator cubic subdomains: bf-marching, separate, generic
!===============================================================================

subroutine CART__TPO_Schwarz_Cubic__gen(ns, nc, nd, S, W, cfg, D_inv, f, u)

  !-----------------------------------------------------------------------------
  ! modules

  use Kind_Parameters, only: RNP
  implicit none

  !-----------------------------------------------------------------------------
  ! arguments

  integer,   intent(in)  :: ns                 !< dimension of S
  integer,   intent(in)  :: nc                 !< number of configurations
  integer,   intent(in)  :: nd                 !< number of subdomains
  real(RNP), intent(in)  :: S(ns,ns,nc)        !< 1D eigenvectors
  real(RNP), intent(in)  :: W(ns,nc)           !< 1D weights
  integer,   intent(in)  :: cfg(3,nd)          !< subdomain configurations
  real(RNP), intent(in)  :: D_inv(ns,ns,ns,nd) !< inverse 3D eigenvalues
  real(RNP), intent(in)  :: f(ns,ns,ns,nd)     !< RHS
  real(RNP), intent(out) :: u(ns,ns,ns,nd)     !< solution

  !-----------------------------------------------------------------------------
  ! local variables

  real(RNP), allocatable :: WS_t(:,:,:), z(:,:,:)
  real(RNP) :: tmp

  integer :: i, j, k, l, m
  integer :: c, c1, c2, c3
  integer :: vec_len

  !-----------------------------------------------------------------------------
  ! initialization

  ! OpenACC vector length
  if (ns < 8) then
    vec_len = 128
  else
    vec_len = 256
  end if

  allocate(WS_t(ns,ns,nc), z(ns,ns,ns))

  ! transposed 1D operator
  do c = 1, nc
  do j = 1, ns
  do i = 1, ns
    WS_t(j,i,c) = W(i,c) * S(i,j,c)
  end do
  end do
  end do

  !-----------------------------------------------------------------------------
  ! evaluation

  !$acc data present(cfg,D_inv,f,u) copyin(S,WS_t) async
  !$acc parallel async &
  !$acc & device_type(nvidia) num_workers(1024/vec_len) vector_length(vec_len)
  !$acc loop gang worker private(z)

  !$omp do private(l)
  elements: do l = 1, nd

    ! element-boundary configurations
    c1 = cfg(1,l)
    c2 = cfg(2,l)
    c3 = cfg(3,l)

    ! z = S^t x I x I f ........................................................

    !$acc loop collapse(3) vector
    do k = 1, ns
    do j = 1, ns
    do i = 1, ns
      tmp = 0
      do m = 1, ns
        tmp = tmp + S(m,k,c3) * f(i,j,m,l)
      end do
      z(i,j,k) = tmp
    end do
    end do
    end do

    ! u = I x S^t x I z ........................................................

    !$acc loop collapse(3) vector
    do k = 1, ns
    do j = 1, ns
    do i = 1, ns
      tmp = 0
      do m = 1, ns
        tmp = tmp + S(m,j,c2) * z(i,m,k)
      end do
      u(i,j,k,l) = tmp
    end do
    end do
    end do

    ! z = D⁻¹ (I x I x S^t) u ..................................................

    !$acc loop collapse(3) vector
    do k = 1, ns
    do j = 1, ns
    do i = 1, ns
      tmp = 0
      do m = 1, ns
        tmp = tmp + S(m,i,c1) * u(m,j,k,l)
      end do
      z(i,j,k) = D_inv(i,j,k,l) * tmp
    end do
    end do
    end do

    ! u = I x I x WS z .........................................................

    !$acc loop collapse(3) vector
    do k = 1, ns
    do j = 1, ns
    do i = 1, ns
      tmp = 0
      do m = 1, ns
        tmp = tmp + WS_t(m,i,c1) * z(m,j,k)
      end do
      u(i,j,k,l) = tmp
    end do
    end do
    end do

    ! z = I x WS x I u .........................................................

    !$acc loop collapse(3) vector
    do k = 1, ns
    do j = 1, ns
    do i = 1, ns
      tmp = 0
      do m = 1, ns
        tmp = tmp + WS_t(m,j,c2) * u(i,m,k,l)
      end do
      z(i,j,k) = tmp
    end do
    end do
    end do

    ! u = WS x I x I z .........................................................

    !$acc loop collapse(3) vector
    do k = 1, ns
    do j = 1, ns
    do i = 1, ns
      tmp = 0
      do m = 1, ns
        tmp = tmp + WS_t(m,k,c3) * z(i,j,m)
      end do
      u(i,j,k,l) = tmp
    end do
    end do
    end do

  end do elements
  !$omp end do

  !$acc end parallel
  !$acc end data

!===============================================================================

end subroutine CART__TPO_Schwarz_Cubic__gen
