!-------------------------------------------------------------------------------
!> Performs `v = D (I x I x A^T) u` for na = 12
!>
!>   * k: simple loop (l)
!>   * j: simple loop (l)
!>   * i: blocked with length of 4 (b4)
!>   * p: unrolled and jammed with length of 4 (u4)
!>   * no remainder handling

subroutine SubOp_1d(A, D, u, v)
  !$acc routine vector
  real(RNP), intent(in)  :: D(12,12,12) !< 3D diagonal operator
  real(RNP), intent(in)  :: A(12,12)     !< A
  real(RNP), intent(in)  :: u(12,12,12) !< operand
  real(RNP), intent(out) :: v(12,12,12) !< result

  real(RNP) :: tmp(0:3)
  integer   :: i, j, k

  !$acc loop collapse(2) independent vector
  do k = 1, 12
  do j = 1, 12

    !$acc loop independent vector private(tmp)
    do i = 1, 12, 4
      tmp =       A( 1,i:i+3) * u( 1,j,k)
      tmp = tmp + A( 2,i:i+3) * u( 2,j,k)
      tmp = tmp + A( 3,i:i+3) * u( 3,j,k)
      tmp = tmp + A( 4,i:i+3) * u( 4,j,k)

      v(i:i+3,j,k) = tmp
    end do

    !$acc loop independent vector private(tmp)
    do i = 1, 12, 4
      tmp =       A( 5,i:i+3) * u( 5,j,k)
      tmp = tmp + A( 6,i:i+3) * u( 6,j,k)
      tmp = tmp + A( 7,i:i+3) * u( 7,j,k)
      tmp = tmp + A( 8,i:i+3) * u( 8,j,k)

      v(i:i+3,j,k) = v(i:i+3,j,k) + tmp
    end do

    !$acc loop independent vector private(tmp)
    do i = 1, 12, 4
      tmp =       A( 9,i:i+3) * u( 9,j,k)
      tmp = tmp + A(10,i:i+3) * u(10,j,k)
      tmp = tmp + A(11,i:i+3) * u(11,j,k)
      tmp = tmp + A(12,i:i+3) * u(12,j,k)

      v(i:i+3,j,k) = (v(i:i+3,j,k) + tmp) * D(i:i+3,j,k)
    end do

  end do
  end do

end subroutine SubOp_1d
