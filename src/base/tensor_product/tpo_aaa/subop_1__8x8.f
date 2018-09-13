!-------------------------------------------------------------------------------
!> Performs `z1 = IxIxA u` for na1 = na2 = 8
!>
!>   * k: simple loop (l)
!>   * j: simple loop (l)
!>   * i: blocked with length of 4 (b4)
!>   * p: unrolled and jammed with length of 4 (u4)
!>   * explicit remainder handling (r)

subroutine SubOp_1(At, u, z1)
  !$acc routine vector
  real(RNP), intent(in)  :: At(8,8)   !< transpose of A
  real(RNP), intent(in)  :: u (8,8,8) !< u
  real(RNP), intent(out) :: z1(8,8,8) !< intermediate result

  real(RNP) :: tmp(0:3)
  integer   :: i, j, k

  !$acc collapse(3) independent vector private(tmp)
  do k = 1, 8
  do j = 1, 8
  do i = 1, 8, 4

    tmp =       At(1,i:i+3) * u(1,j,k)
    tmp = tmp + At(2,i:i+3) * u(2,j,k)
    tmp = tmp + At(3,i:i+3) * u(3,j,k)
    tmp = tmp + At(4,i:i+3) * u(4,j,k)
    tmp = tmp + At(5,i:i+3) * u(5,j,k)
    tmp = tmp + At(6,i:i+3) * u(6,j,k)
    tmp = tmp + At(7,i:i+3) * u(7,j,k)
    tmp = tmp + At(8,i:i+3) * u(8,j,k)

    z1(i:i+3,j,k) = tmp

  end do
  end do
  end do

end subroutine SubOp_1
