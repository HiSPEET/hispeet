!-------------------------------------------------------------------------------
!> Performs `v = g IxAxI u`
!>
!>   * k: simple loop (l)
!>   * j: unrolled and jammed with length of 8 (u8)
!>   * i: simple loop (l)
!>   * p: simple loop (l)
!>   * explicit remainder handling (r)

subroutine SubOp_2i(g, At, u, v)
  !$acc routine vector
  real(RNP), intent(in)  :: g                        !< metric factor
  real(RNP), intent(in)  :: At(__NA__,__NA__)        !< transpose of A
  real(RNP), intent(in)  :: u(__NA__,__NA__,__NA__)  !< operand
  real(RNP), intent(out) :: v(__NA__,__NA__,__NA__)  !< result

  real(RNP) :: tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7
  integer   :: i, j, k, p

#if __NA_T8__ > 0

  !$acc loop collapse(3) vector
  do k = 1, __NA__
  do j = 1, __NA_T8__, 8
  do i = 1, __NA__
    tmp0 = 0
    tmp1 = 0
    tmp2 = 0
    tmp3 = 0
    tmp4 = 0
    tmp5 = 0
    tmp6 = 0
    tmp7 = 0
    do p = 1, __NA__
      tmp0 = tmp0 + At(p,j  ) * u(i,p,k)
      tmp1 = tmp1 + At(p,j+1) * u(i,p,k)
      tmp2 = tmp2 + At(p,j+2) * u(i,p,k)
      tmp3 = tmp3 + At(p,j+3) * u(i,p,k)
      tmp4 = tmp4 + At(p,j+4) * u(i,p,k)
      tmp5 = tmp5 + At(p,j+5) * u(i,p,k)
      tmp6 = tmp6 + At(p,j+6) * u(i,p,k)
      tmp7 = tmp7 + At(p,j+7) * u(i,p,k)
    end do
    v(i,j  ,k) = g * tmp0
    v(i,j+1,k) = g * tmp1
    v(i,j+2,k) = g * tmp2
    v(i,j+3,k) = g * tmp3
    v(i,j+4,k) = g * tmp4
    v(i,j+5,k) = g * tmp5
    v(i,j+6,k) = g * tmp6
    v(i,j+7,k) = g * tmp7
  end do
  end do
  end do

#endif

#if __NA__ == __NA_T8__ + 1

  ! remainder: j = na ........................................................

  j = __NA_T8__ + 1

  !$acc loop collapse(2) vector
  do k = 1, __NA__
  do i = 1, __NA__
    tmp0 = 0
    do p = 1, __NA__
      tmp0 = tmp0 + At(p,j  ) * u(i,p,k)
    end do
    v(i,j,k) = g * tmp0
  end do
  end do

#elif __NA__ == __NA_T8__ + 2

  ! remainder: j = na-1:na ...................................................

  j = __NA_T8__ + 1

  !$acc loop collapse(2) vector
  do k = 1, __NA__
  do i = 1, __NA__
    tmp0 = 0
    tmp1 = 0
    do p = 1, __NA__
      tmp0 = tmp0 + At(p,j  ) * u(i,p,k)
      tmp1 = tmp1 + At(p,j+1) * u(i,p,k)
    end do
    v(i,j  ,k) = g * tmp0
    v(i,j+1,k) = g * tmp1
  end do
  end do

#elif __NA__ == __NA_T8__ + 3

  ! remainder: j = na-2:na ...................................................

  j = __NA_T8__ + 1

  !$acc loop collapse(2) vector
  do k = 1, __NA__
  do i = 1, __NA__
    tmp0 = 0
    tmp1 = 0
    tmp2 = 0
    do p = 1, __NA__
      tmp0 = tmp0 + At(p,j  ) * u(i,p,k)
      tmp1 = tmp1 + At(p,j+1) * u(i,p,k)
      tmp2 = tmp2 + At(p,j+2) * u(i,p,k)
    end do
    v(i,j  ,k) = g * tmp0
    v(i,j+1,k) = g * tmp1
    v(i,j+2,k) = g * tmp2
  end do
  end do

#elif __NA__ == __NA_T8__ + 4

  ! remainder: j = na-3:na ...................................................

  j = __NA_T8__ + 1

  !$acc loop collapse(2) vector
  do k = 1, __NA__
  do i = 1, __NA__
    tmp0 = 0
    tmp1 = 0
    tmp2 = 0
    tmp3 = 0
    do p = 1, __NA__
      tmp0 = tmp0 + At(p,j  ) * u(i,p,k)
      tmp1 = tmp1 + At(p,j+1) * u(i,p,k)
      tmp2 = tmp2 + At(p,j+2) * u(i,p,k)
      tmp3 = tmp3 + At(p,j+3) * u(i,p,k)
    end do
    v(i,j  ,k) = g * tmp0
    v(i,j+1,k) = g * tmp1
    v(i,j+2,k) = g * tmp2
    v(i,j+3,k) = g * tmp3
  end do
  end do

#elif __NA__ == __NA_T8__ + 5

  ! remainder: j = na-4:na ...................................................

  j = __NA_T8__ + 1

  !$acc loop collapse(2) vector
  do k = 1, __NA__
  do i = 1, __NA__
    tmp0 = 0
    tmp1 = 0
    tmp2 = 0
    tmp3 = 0
    tmp4 = 0
    do p = 1, __NA__
      tmp0 = tmp0 + At(p,j  ) * u(i,p,k)
      tmp1 = tmp1 + At(p,j+1) * u(i,p,k)
      tmp2 = tmp2 + At(p,j+2) * u(i,p,k)
      tmp3 = tmp3 + At(p,j+3) * u(i,p,k)
      tmp4 = tmp4 + At(p,j+4) * u(i,p,k)
    end do
    v(i,j  ,k) = g * tmp0
    v(i,j+1,k) = g * tmp1
    v(i,j+2,k) = g * tmp2
    v(i,j+3,k) = g * tmp3
    v(i,j+4,k) = g * tmp4
  end do
  end do

#elif __NA__ == __NA_T8__ + 6

  ! remainder: j = na-5:na ...................................................

  j = __NA_T8__ + 1

  !$acc loop collapse(2) vector
  do k = 1, __NA__
  do i = 1, __NA__
    tmp0 = 0
    tmp1 = 0
    tmp2 = 0
    tmp3 = 0
    tmp4 = 0
    tmp5 = 0
    do p = 1, __NA__
      tmp0 = tmp0 + At(p,j  ) * u(i,p,k)
      tmp1 = tmp1 + At(p,j+1) * u(i,p,k)
      tmp2 = tmp2 + At(p,j+2) * u(i,p,k)
      tmp3 = tmp3 + At(p,j+3) * u(i,p,k)
      tmp4 = tmp4 + At(p,j+4) * u(i,p,k)
      tmp5 = tmp5 + At(p,j+5) * u(i,p,k)
    end do
    v(i,j  ,k) = g * tmp0
    v(i,j+1,k) = g * tmp1
    v(i,j+2,k) = g * tmp2
    v(i,j+3,k) = g * tmp3
    v(i,j+4,k) = g * tmp4
    v(i,j+5,k) = g * tmp5
  end do
  end do

#elif __NA__ == __NA_T8__ + 7

  ! remainder: j = na-6:na ...................................................

  j = __NA_T8__ + 1

  !$acc loop collapse(2) vector
  do k = 1, __NA__
  do i = 1, __NA__
    tmp0 = 0
    tmp1 = 0
    tmp2 = 0
    tmp3 = 0
    tmp4 = 0
    tmp5 = 0
    tmp6 = 0
    do p = 1, __NA__
      tmp0 = tmp0 + At(p,j  ) * u(i,p,k)
      tmp1 = tmp1 + At(p,j+1) * u(i,p,k)
      tmp2 = tmp2 + At(p,j+2) * u(i,p,k)
      tmp3 = tmp3 + At(p,j+3) * u(i,p,k)
      tmp4 = tmp4 + At(p,j+4) * u(i,p,k)
      tmp5 = tmp5 + At(p,j+5) * u(i,p,k)
      tmp6 = tmp6 + At(p,j+6) * u(i,p,k)
    end do
    v(i,j  ,k) = g * tmp0
    v(i,j+1,k) = g * tmp1
    v(i,j+2,k) = g * tmp2
    v(i,j+3,k) = g * tmp3
    v(i,j+4,k) = g * tmp4
    v(i,j+5,k) = g * tmp5
    v(i,j+6,k) = g * tmp6
  end do
  end do

#endif

end subroutine SubOp_2i
