!-------------------------------------------------------------------------------
!> Performs `v = g IxIxAᵀ u` for na = 20
!>
!>   * k: simple loop (l)
!>   * j: simple loop (l)
!>   * i: blocked with length of 4 (b4)
!>   * p: unrolled and jammed with length of 4 (u4)
!>   * no remainder handling

subroutine SubOp_1i(g, A, u, v)
  !$acc routine vector
  real(RNP), intent(in)  :: g                       !< metric factor
  real(RNP), intent(in)  :: A(__NA__,__NA__)        !< A
  real(RNP), intent(in)  :: u(__NA__,__NA__,__NA__) !< operand
  real(RNP), intent(out) :: v(__NA__,__NA__,__NA__) !< result

  real(RNP) :: tmp(0:3)
  integer   :: i, j, k

  !$acc loop collapse(2) independent vector
  do k = 1, __NA__
  do j = 1, __NA__

    !$acc loop independent vector private(tmp)
    do i = 1, __NA_T4__, 4
      tmp =       A( 1,i:i+3) * u( 1,j,k)
      tmp = tmp + A( 2,i:i+3) * u( 2,j,k)
      tmp = tmp + A( 3,i:i+3) * u( 3,j,k)
      tmp = tmp + A( 4,i:i+3) * u( 4,j,k)

      v(i:i+3,j,k) = g * tmp
    end do

    !$acc loop independent vector private(tmp)
    do i = 1, __NA_T4__, 4
      tmp =       A( 5,i:i+3) * u( 5,j,k)
      tmp = tmp + A( 6,i:i+3) * u( 6,j,k)
      tmp = tmp + A( 7,i:i+3) * u( 7,j,k)
      tmp = tmp + A( 8,i:i+3) * u( 8,j,k)

      v(i:i+3,j,k) = v(i:i+3,j,k) + g * tmp
    end do

    !$acc loop independent vector private(tmp)
    do i = 1, __NA_T4__, 4
      tmp =       A( 9,i:i+3) * u( 9,j,k)
      tmp = tmp + A(10,i:i+3) * u(10,j,k)
      tmp = tmp + A(11,i:i+3) * u(11,j,k)
      tmp = tmp + A(12,i:i+3) * u(12,j,k)

      v(i:i+3,j,k) = v(i:i+3,j,k) + g * tmp
    end do

    !$acc loop independent vector private(tmp)
    do i = 1, __NA_T4__, 4
      tmp =       A(13,i:i+3) * u(13,j,k)
      tmp = tmp + A(14,i:i+3) * u(14,j,k)
      tmp = tmp + A(15,i:i+3) * u(15,j,k)
      tmp = tmp + A(16,i:i+3) * u(16,j,k)

      v(i:i+3,j,k) = v(i:i+3,j,k) + g * tmp
    end do

    !$acc loop independent vector private(tmp)
    do i = 1, __NA_T4__, 4
      tmp =       A(17,i:i+3) * u(17,j,k)
      tmp = tmp + A(18,i:i+3) * u(18,j,k)
      tmp = tmp + A(19,i:i+3) * u(19,j,k)
      tmp = tmp + A(20,i:i+3) * u(20,j,k)

      v(i:i+3,j,k) = v(i:i+3,j,k) + g * tmp
    end do

  end do
  end do

end subroutine SubOp_1i
