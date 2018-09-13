!-------------------------------------------------------------------------------
!> Performs `v = g2 IxAxI u`
!>
!>   * k: simple loop (l)
!>   * j: unrolled and jammed with length of 2 (u2)
!>   * i: simple loop (l)
!>   * p: simple loop (l)
!>   * explicit remainder handling (r)

subroutine SubOp_2i(g2, At, u, v)
  !$acc routine vector
  real(RNP), intent(in)  :: g2                       !< metric factor
  real(RNP), intent(in)  :: At(__NA__,__NA__)        !< transpose of A
  real(RNP), intent(in)  :: u(__NA__,__NA__,__NA__)  !< operand
  real(RNP), intent(out) :: v(__NA__,__NA__,__NA__)  !< result

  real(RNP) :: tmp0, tmp1
  integer   :: i, j, k, p

  !$acc loop collapse(3) vector
  do k = 1, __NA__
  do j = 1, __NA_T2__, 2
  do i = 1, __NA__
    tmp0 = 0
    tmp1 = 0
    do p = 1, __NA__
      tmp0 = tmp0 + At(p,j  ) * u(i,p,k)
      tmp1 = tmp1 + At(p,j+1) * u(i,p,k)
    end do
    v(i,j  ,k) = g2 * tmp0
    v(i,j+1,k) = g2 * tmp1
  end do
  end do
  end do

#if __NA__ == __NA_T2__ + 1

  !$acc loop collapse(2) vector
  do k = 1, __NA__
  do i = 1, __NA__
    tmp0 = 0
    do p = 1, __NA__
      tmp0 = tmp0 + At(p,__NA__) * u(i,p,k)
    end do
    v(i,__NA__,k) = g2 * tmp0
  end do
  end do

#endif

end subroutine SubOp_2i
