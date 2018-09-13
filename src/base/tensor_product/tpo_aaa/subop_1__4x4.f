!-------------------------------------------------------------------------------
!> Performs `z1 = IxIxA u` for na1 = na2 = 4
!>
!>   * k: simple loop (l)
!>   * j: simple loop (l)
!>   * i: blocked with length of 4 (b4)
!>   * p: unrolled and jammed with length of 4 (u4)
!>   * explicit remainder handling (r)

subroutine SubOp_1(At, u, z1)
  !$acc routine vector
  real(RNP), intent(in)  :: At(4,4)   !< transpose of A
  real(RNP), intent(in)  :: u (4,4,4) !< u
  real(RNP), intent(out) :: z1(4,4,4) !< intermediate result

  real(RNP) :: tmp(0:3)
  integer   :: j, k

  !$acc collapse(2) independent vector private(tmp)
  do k = 1, 4
  do j = 1, 4

    tmp =       At(1,1:4) * u(1,j,k)
    tmp = tmp + At(2,1:4) * u(2,j,k)
    tmp = tmp + At(3,1:4) * u(3,j,k)
    tmp = tmp + At(4,1:4) * u(4,j,k)

    z1(1:4,j,k) = tmp

  end do
  end do

end subroutine SubOp_1
