!> summary:  Spectral operator backward-forward marching, separate, generic
!> author:   Joerg Stiller
!> date:     2017/07/19
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Spectral operator backward-forward marching, separate, generic
!===============================================================================

subroutine CART__TPO_Spectral_Iso__gen(na, ne, A, Lambda, u, v)

  !-----------------------------------------------------------------------------
  ! modules

  use Kind_Parameters, only: RNP
  implicit none

  !-----------------------------------------------------------------------------
  ! arguments

  integer,   intent(in)  :: na               !< dimension of A
  integer,   intent(in)  :: ne               !< number of elements or regions
  real(RNP), intent(in)  :: A(na,na)         !< 1D eigenvectors
  real(RNP), intent(in)  :: Lambda(na,na,na) !< 3D eigenvalues
  real(RNP), intent(in)  :: u(na,na,na,ne)   !< operand
  real(RNP), intent(out) :: v(na,na,na,ne)   !< result

  !-----------------------------------------------------------------------------
  ! local variables

  real(RNP), allocatable :: At(:,:), z(:,:,:)
  real(RNP) :: tmp

  integer :: e, i, j, k, m
  integer :: vec_len

  !-----------------------------------------------------------------------------
  ! initialization

  ! OpenACC vector length
  if (na < 8) then
    vec_len = 128
  else
    vec_len = 256
  end if

  allocate(At(na,na), z(na,na,na))

  ! transposed 1D operator
  do j = 1, na
  do i = 1, na
    At(j,i) = A(i,j)
  end do
  end do

  !-----------------------------------------------------------------------------
  ! evaluation

  !$acc data present(u,v) copyin(A,At,Lambda) async
  !$acc parallel async &
  !$acc & device_type(nvidia) num_workers(1024/vec_len) vector_length(vec_len)
  !$acc loop gang worker private(z)

  !$omp do private(e)
  elements: do e = 1, ne

    ! z = A^t x I x I u ........................................................

    !$acc loop collapse(3) vector
    do k = 1, na
    do j = 1, na
    do i = 1, na
      tmp = 0
      do m = 1, na
        tmp = tmp + A(m,k) * u(i,j,m,e)
      end do
      z(i,j,k) = tmp
    end do
    end do
    end do

    ! v = I x A^t x I z ........................................................

    !$acc loop collapse(3) vector
    do k = 1, na
    do j = 1, na
    do i = 1, na
      tmp = 0
      do m = 1, na
        tmp = tmp + A(m,j) * z(i,m,k)
      end do
      v(i,j,k,e) = tmp
    end do
    end do
    end do

    ! z = Lambda I x I x A^t v .................................................

    !$acc loop collapse(3) vector
    do k = 1, na
    do j = 1, na
    do i = 1, na
      tmp = 0
      do m = 1, na
        tmp = tmp + A(m,i) * v(m,j,k,e)
      end do
      z(i,j,k) = tmp * Lambda(i,j,k)
    end do
    end do
    end do

    ! v = I x I x A z ..........................................................

    !$acc loop collapse(3) vector
    do k = 1, na
    do j = 1, na
    do i = 1, na
      tmp = 0
      do m = 1, na
        tmp = tmp + At(m,i) * z(m,j,k)
      end do
      v(i,j,k,e) = tmp
    end do
    end do
    end do

    ! z = I x A x I v ..........................................................

    !$acc loop collapse(3) vector
    do k = 1, na
    do j = 1, na
    do i = 1, na
      tmp = 0
      do m = 1, na
        tmp = tmp + At(m,j) * v(i,m,k,e)
      end do
      z(i,j,k) = tmp
    end do
    end do
    end do

    ! v = A x I x I z ..........................................................

    !$acc loop collapse(3) vector
    do k = 1, na
    do j = 1, na
    do i = 1, na
      tmp = 0
      do m = 1, na
        tmp = tmp + At(m,k) * z(i,j,m)
      end do
      v(i,j,k,e) = tmp
    end do
    end do
    end do

  end do elements
  !$omp end do

  !$acc end parallel
  !$acc end data

!===============================================================================

end subroutine CART__TPO_Spectral_Iso__gen
