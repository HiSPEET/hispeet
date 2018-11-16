!-------------------------------------------------------------------------------
!> Performs `v = g IxIxA u`
!>
!>   * k: simple loop (l)
!>   * j: unrolled and jammed with length of 2 (u2)
!>   * i: unrolled and jammed with length of 2 (u2)
!>   * p: simple loop (l)
!>   * explicit remainder handling (r)

subroutine SubOp_1i(g, At, u, v)
  !$acc routine vector
  real(RNP), intent(in)  :: g                       !< metric factor
  real(RNP), intent(in)  :: At(__NA__,__NA__)       !< transpose of A
  real(RNP), intent(in)  :: u(__NA__,__NA__,__NA__) !< operand
  real(RNP), intent(out) :: v(__NA__,__NA__,__NA__) !< result

  real(RNP) :: tmp0, tmp1, tmp2, tmp3
  integer   :: i, j, k, p

  !$acc loop collapse(3) vector
  do k = 1, __NA__
  do j = 1, __NA_T2__, 2
  do i = 1, __NA_T2__, 2
    tmp0 = 0
    tmp1 = 0
    tmp2 = 0
    tmp3 = 0
    do p = 1, __NA__
      tmp0 = tmp0 + At(p,i  ) * u(p,j  ,k)
      tmp1 = tmp1 + At(p,i+1) * u(p,j  ,k)
      tmp2 = tmp2 + At(p,i  ) * u(p,j+1,k)
      tmp3 = tmp3 + At(p,i+1) * u(p,j+1,k)
    end do
    v(i  ,j  ,k) = g * tmp0
    v(i+1,j  ,k) = g * tmp1
    v(i  ,j+1,k) = g * tmp2
    v(i+1,j+1,k) = g * tmp3
  end do
  end do
  end do

#if __NA__ == __NA_T2__ + 1

  ! remainder: i = 1:na-1, j = na ............................................

  j = __NA__

  !$acc loop collapse(2) vector
  do k = 1, __NA__
  do i = 1, __NA_T2__, 2
    tmp0 = 0
    tmp1 = 0
    do p = 1, __NA__
      tmp0 = tmp0 + At(p,i  ) * u(p,j,k)
      tmp1 = tmp1 + At(p,i+1) * u(p,j,k)
    end do
    v(i  ,j,k) = g * tmp0
    v(i+1,j,k) = g * tmp1
  end do
  end do

  ! remainder: i = na, j = 1:na ..............................................

  i = __NA__

  !$acc loop collapse(2) vector
  do k = 1, __NA__
  do j = 1, __NA__
    tmp0 = 0
    do p = 1, __NA__
      tmp0 = tmp0 + At(p,i) * u(p,j,k)
    end do
    v(i,j,k) = g * tmp0
  end do
  end do

#endif

end subroutine SubOp_1i
