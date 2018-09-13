!-------------------------------------------------------------------------------
!> Performs `v = D (I x I x A^T) u`
!>
!>   * k: simple loop (l)
!>   * j: simple loop (l)
!>   * i: blocked with length of 4 (b4)
!>   * p: unrolled and jammed with length of 4 (u4)
!>   * no remainder handling

subroutine SubOp_1d(A, D, u, v)
  !$acc routine vector
  real(RNP), intent(in)  :: D(8,8,8) !< 3D diagonal operator
  real(RNP), intent(in)  :: A(8,8)   !< A
  real(RNP), intent(in)  :: u(8,8,8) !< operand
  real(RNP), intent(out) :: v(8,8,8) !< result

  real(RNP) :: tmp(0:3)
  integer   :: i, j, k

  !$acc loop collapse(3) independent vector private(tmp)
  do k = 1, 8
    do j = 1, 8
    do i = 1, 8, 4

      tmp =       A(1,i:i+3) * u(1,j,k)
      tmp = tmp + A(2,i:i+3) * u(2,j,k)
      tmp = tmp + A(3,i:i+3) * u(3,j,k)
      tmp = tmp + A(4,i:i+3) * u(4,j,k)
      tmp = tmp + A(5,i:i+3) * u(5,j,k)
      tmp = tmp + A(6,i:i+3) * u(6,j,k)
      tmp = tmp + A(7,i:i+3) * u(7,j,k)
      tmp = tmp + A(8,i:i+3) * u(8,j,k)

      v(i:i+3,j,k) = tmp * D(i:i+3,j,k)

    end do
    end do
  end do

end subroutine SubOp_1d
