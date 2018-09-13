!-------------------------------------------------------------------------------
!> Performs `v = D (I x I x A^T) u`
!>
!>   * k: simple loop (l)
!>   * j: simple loop (l)
!>   * i: simple loop (l)
!>   * p: simple loop (l)
!>   * using Intel SIMD directive

subroutine SubOp_1d(A, D, u, v)
  !$acc routine vector
  real(RNP), intent(in)  :: D(__NA__,__NA__,__NA__) !< 3D diagonal operator
  real(RNP), intent(in)  :: A(__NA__,__NA__)        !< A
  real(RNP), intent(in)  :: u(__NA__,__NA__,__NA__) !< operand
  real(RNP), intent(out) :: v(__NA__,__NA__,__NA__) !< result

  real(RNP) :: tmp
  integer   :: i, j, k, p

  !$acc loop collapse(3) independent vector
  do k = 1, __NA__
  do j = 1, __NA__
    !DIR$ SIMD VECREMAINDER
    do i = 1, __NA__
      tmp = 0
      do p = 1, __NA__
        tmp = tmp + A(p,i) * u(p,j,k)
      end do
      v(i,j,k) = tmp * D(i,j,k)
    end do
  end do
  end do

end subroutine SubOp_1d
