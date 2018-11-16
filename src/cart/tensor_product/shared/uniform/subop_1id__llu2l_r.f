!-------------------------------------------------------------------------------
!> Performs `v = g D IxIxA u`
!>
!>   * k: simple loop (l)
!>   * j: simple loop (l)
!>   * i: unrolled and jammed with length of 2 (u2)
!>   * p: simple loop (l)
!>   * explicit remainder handling (r)

subroutine SubOp_1id(g, At, D, u, v)
  !$acc routine vector
  real(RNP), intent(in)  :: g                       !< metric factor
  real(RNP), intent(in)  :: At(__NA__,__NA__)       !< transpose of A
  real(RNP), intent(in)  :: D(__NA__,__NA__,__NA__) !< diagonal operator
  real(RNP), intent(in)  :: u(__NA__,__NA__,__NA__) !< operand
  real(RNP), intent(out) :: v(__NA__,__NA__,__NA__) !< result

  real(RNP) :: tmp0, tmp1
  integer   :: i, j, k, p

  !$acc loop collapse(3) vector
  do k = 1, __NA__
  do j = 1, __NA__
  do i = 1, __NA_T2__, 2
    tmp0 = 0
    tmp1 = 0
    do p = 1, __NA__
      tmp0 = tmp0 + At(p,i  ) * u(p,j,k)
      tmp1 = tmp1 + At(p,i+1) * u(p,j,k)
    end do
    v(i  ,j,k) = g * D(i  ,j,k) * tmp0
    v(i+1,j,k) = g * D(i+1,j,k) * tmp1
  end do
  end do
  end do

#if __NA__ == __NA_T2__ + 1
  !$acc loop collapse(2) vector
  do k = 1, __NA__
  do j = 1, __NA__
    tmp0 = 0
    do p = 1, __NA__
      tmp0 = tmp0 + At(p,__NA__) * u(p,j,k)
    end do
    v(__NA__,j,k) = g * D(__NA__,j,k) * tmp0
  end do
  end do
#endif

end subroutine SubOp_1id
