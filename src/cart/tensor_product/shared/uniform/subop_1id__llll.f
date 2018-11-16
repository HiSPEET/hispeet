!-------------------------------------------------------------------------------
!> Performs `v = g D IxIxA u`
!>
!>   * k: simple loop (l)
!>   * j: simple loop (l)
!>   * i: simple loop (l)
!>   * p: simple loop (l)
!>   * using Intel SIMD directive

subroutine SubOp_1id(g, At, D, u, v)
  !$acc routine vector
  real(RNP), intent(in)  :: g                       !< metric factor
  real(RNP), intent(in)  :: At(__NA__,__NA__)       !< transpose of A
  real(RNP), intent(in)  :: D(__NA__,__NA__,__NA__) !< diagonal operator
  real(RNP), intent(in)  :: u(__NA__,__NA__,__NA__) !< operand
  real(RNP), intent(out) :: v(__NA__,__NA__,__NA__) !< result

  real(RNP) :: tmp
  integer   :: i, j, k, p

  !$acc loop collapse(3) vector
  do k = 1, __NA__
  do j = 1, __NA__
    !DIR$ SIMD VECREMAINDER
    do i = 1, __NA__
      tmp = 0
      do p = 1, __NA__
        tmp = tmp + At(p,i) * u(p,j,k)
      end do
      v(i,j,k) = g * D(i,j,k) * tmp
    end do
  end do
  end do

end subroutine SubOp_1id
