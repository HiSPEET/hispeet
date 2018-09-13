!> summary:  Generic AxAxA operator with separate loop nests
!> author:   Joerg Stiller
!> date:     2017/05/16
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Generic AxAxA operator with separate loop nests
!===============================================================================

subroutine TPO_AAA__gen(na1, na2, ne, A, u, v)

  !-----------------------------------------------------------------------------
  ! used modules

  use Kind_Parameters, only: RNP
  implicit none

  !-----------------------------------------------------------------------------
  ! arguments

  integer,   intent(in)  :: na1               !< first dimension of A
  integer,   intent(in)  :: na2               !< second dimension of A
  integer,   intent(in)  :: ne                !< number of elements
  real(RNP), intent(in)  :: A(na1,na2)        !< 1D operator
  real(RNP), intent(in)  :: u(na2,na2,na2,ne) !< 3D operand
  real(RNP), intent(out) :: v(na1,na1,na1,ne) !< 3D result

  !-----------------------------------------------------------------------------
  ! local variables

  real(RNP), parameter :: ZERO = 0
  real(RNP), allocatable :: At(:,:), z2(:,:,:), z3(:,:,:)
  real(RNP) :: tmp

  integer :: e, i, j, k, m
  integer :: vec_len

  !-----------------------------------------------------------------------------
  ! initialization

  ! OpenACC vector length
  if (min(na1,na2) < 8) then
    vec_len = 128
  else
    vec_len = 256
  end if

  allocate(At(na2,na1), z2(na2,na1,na1), z3(na2,na2,na1))
  At = transpose(A)

  !-----------------------------------------------------------------------------
  ! evaluation

  !$acc data present(u,v) copyin(At) async
  !$acc parallel async &
  !$acc & device_type(nvidia) num_workers(1024/vec_len) vector_length(vec_len)
  !$acc loop gang worker private(z2,z3)

  !$omp do private(e)
  do e = 1, ne

    ! z3 = AxIxI u^e ...........................................................

    !$acc loop collapse(3) vector
    do k = 1, na1
    do j = 1, na2
    !DIR$ SIMD
    do i = 1, na2
      tmp = ZERO
      do m = 1, na2
        tmp = tmp + At(m,k) * u(i,j,m,e)
      end do
      z3(i,j,k) = tmp
    end do
    end do
    end do

    ! z2 = IxAxI z3 ............................................................

    !$acc loop collapse(3) vector
    do k = 1, na1
    do j = 1, na1
    !DIR$ SIMD
    do i = 1, na2
      tmp = ZERO
      do m = 1, na2
        tmp = tmp + At(m,j) * z3(i,m,k)
      end do
      z2(i,j,k) = tmp
    end do
    end do
    end do

    ! v^e = IxIxA z2 ...........................................................

    !$acc loop collapse(3) vector
    do k = 1, na1
    do j = 1, na1
    !DIR$ SIMD
    do i = 1, na1
      tmp = ZERO
      do m = 1, na2
        tmp = tmp + At(m,i) * z2(m,j,k)
      end do
      v(i,j,k,e) = tmp
    end do
    end do
    end do

  end do
  !$omp end do

  !$acc end parallel
  !$acc end data

!===============================================================================

end subroutine TPO_AAA__gen
