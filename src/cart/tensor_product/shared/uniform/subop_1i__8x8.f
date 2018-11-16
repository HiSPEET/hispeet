!-------------------------------------------------------------------------------
!> Performs `v = g IxIxA u` for na = 8
!>
!>   * k: simple loop (l)
!>   * j: simple loop (l)
!>   * i: blocked with length of 4 (b4)
!>   * p: unrolled and jammed with length of 4 (u4)
!>   * no remainder handling

subroutine SubOp_1i(g, At, u, v)
  !$acc routine vector
  real(RNP), intent(in)  :: g                       !< metric factor
  real(RNP), intent(in)  :: At(__NA__,__NA__)       !< transpose of A
  real(RNP), intent(in)  :: u(__NA__,__NA__,__NA__) !< operand
  real(RNP), intent(out) :: v(__NA__,__NA__,__NA__) !< result

  real(RNP) :: tmp(0:3)
  integer   :: i, j, k

  !$acc loop collapse(3) independent vector private(tmp)
  do k = 1, __NA__
    do j = 1, __NA__
    do i = 1, __NA_T4__, 4

      tmp =       At(1,i:i+3) * u(1,j,k)
      tmp = tmp + At(2,i:i+3) * u(2,j,k)
      tmp = tmp + At(3,i:i+3) * u(3,j,k)
      tmp = tmp + At(4,i:i+3) * u(4,j,k)
      tmp = tmp + At(5,i:i+3) * u(5,j,k)
      tmp = tmp + At(6,i:i+3) * u(6,j,k)
      tmp = tmp + At(7,i:i+3) * u(7,j,k)
      tmp = tmp + At(8,i:i+3) * u(8,j,k)

      v(i:i+3,j,k) = g * tmp

    end do
    end do
  end do

end subroutine SubOp_1i
