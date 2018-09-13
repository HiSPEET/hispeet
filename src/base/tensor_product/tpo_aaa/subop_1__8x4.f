!-------------------------------------------------------------------------------
!> Performs `z1 = IxA u` for fixed k and na1 = 8, na2 = 4
!>
!>   * j: simple loop (l)
!>   * i: blocked with length of 4 (b4)
!>   * p: unrolled and jammed with length of 4 (u4)
!>   * no remainder handling

subroutine SubOp_1(At, u, z1)
  !$acc routine vector
  real(RNP), intent(in)  :: At(4,8)  !< transpose of A
  real(RNP), intent(in)  :: u (4,4)  !< k-slice of u
  real(RNP), intent(out) :: z1(8,4)  !< intermediate result

  real(RNP) :: tmp(0:3)
  integer   :: i, j

  !$acc loop collapse(2) independent vector private(tmp)
  do j = 1, 4
  do i = 1, 8, 4

    tmp =       At(1,i:i+3) * u(1,j)
    tmp = tmp + At(2,i:i+3) * u(2,j)
    tmp = tmp + At(3,i:i+3) * u(3,j)
    tmp = tmp + At(4,i:i+3) * u(4,j)

    z1(i:i+3,j) = tmp

  end do
  end do

end subroutine SubOp_1
