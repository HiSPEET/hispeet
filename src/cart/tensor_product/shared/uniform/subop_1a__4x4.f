!-------------------------------------------------------------------------------
!> Performs `v += g IxIxAᵀ u` for ns = 4
!>
!>   * k: simple loop (l)
!>   * j: simple loop (l)
!>   * i: blocked with length of 4 (b4)
!>   * p: unrolled and jammed with length of 4 (u4)
!>   * no remainder handling

subroutine SubOp_1a(g, A, u, v)
  !$acc routine vector
  real(RNP), intent(in)    :: g                       !< metric factor
  real(RNP), intent(in)    :: A(__NA__,__NA__)        !< A
  real(RNP), intent(in)    :: u(__NA__,__NA__,__NA__) !< operand
  real(RNP), intent(inout) :: v(__NA__,__NA__,__NA__) !< result

  real(RNP) :: tmp(0:3)
  integer   :: i, j, k

  !$acc loop collapse(3) independent vector private(tmp)
  do k = 1, __NA__
    do j = 1, __NA__
    do i = 1, __NA_T4__, 4

      tmp =       A(1,i:i+3) * u(1,j,k)
      tmp = tmp + A(2,i:i+3) * u(2,j,k)
      tmp = tmp + A(3,i:i+3) * u(3,j,k)
      tmp = tmp + A(4,i:i+3) * u(4,j,k)

      v(i:i+3,j,k) = v(i:i+3,j,k) + g * tmp

    end do
    end do
  end do

end subroutine SubOp_1a
