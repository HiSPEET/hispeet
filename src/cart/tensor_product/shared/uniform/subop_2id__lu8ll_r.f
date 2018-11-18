!-------------------------------------------------------------------------------
!> Performs `v = g D IxAᵀxI u`
!>
!>   * k: simple loop (l)
!>   * j: unrolled and jammed with length of 8 (u8)
!>   * i: simple loop (l)
!>   * p: simple loop (l)
!>   * explicit remainder handling (r)

subroutine SubOp_2id(g, A, D, u, v)
  !$acc routine vector
  real(RNP), intent(in)  :: g                        !< metric factor
  real(RNP), intent(in)  :: A(__NA__,__NA__)         !< A
  real(RNP), intent(in)  :: D(__NA__,__NA__,__NA__)  !< diagonal operator
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
      tmp0 = tmp0 + A(p,j  ) * u(i,p,k)
      tmp1 = tmp1 + A(p,j+1) * u(i,p,k)
      tmp2 = tmp2 + A(p,j+2) * u(i,p,k)
      tmp3 = tmp3 + A(p,j+3) * u(i,p,k)
      tmp4 = tmp4 + A(p,j+4) * u(i,p,k)
      tmp5 = tmp5 + A(p,j+5) * u(i,p,k)
      tmp6 = tmp6 + A(p,j+6) * u(i,p,k)
      tmp7 = tmp7 + A(p,j+7) * u(i,p,k)
    end do
    v(i,j  ,k) = g * D(i,j  ,k) * tmp0
    v(i,j+1,k) = g * D(i,j+1,k) * tmp1
    v(i,j+2,k) = g * D(i,j+2,k) * tmp2
    v(i,j+3,k) = g * D(i,j+3,k) * tmp3
    v(i,j+4,k) = g * D(i,j+4,k) * tmp4
    v(i,j+5,k) = g * D(i,j+5,k) * tmp5
    v(i,j+6,k) = g * D(i,j+6,k) * tmp6
    v(i,j+7,k) = g * D(i,j+7,k) * tmp7
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
      tmp0 = tmp0 + A(p,j  ) * u(i,p,k)
    end do
    v(i,j,k) = g * D(i,j,k) * tmp0
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
      tmp0 = tmp0 + A(p,j  ) * u(i,p,k)
      tmp1 = tmp1 + A(p,j+1) * u(i,p,k)
    end do
    v(i,j  ,k) = g * D(i,j  ,k) * tmp0
    v(i,j+1,k) = g * D(i,j+1,k) * tmp1
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
      tmp0 = tmp0 + A(p,j  ) * u(i,p,k)
      tmp1 = tmp1 + A(p,j+1) * u(i,p,k)
      tmp2 = tmp2 + A(p,j+2) * u(i,p,k)
    end do
    v(i,j  ,k) = g * D(i,j  ,k) * tmp0
    v(i,j+1,k) = g * D(i,j+1,k) * tmp1
    v(i,j+2,k) = g * D(i,j+2,k) * tmp2
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
      tmp0 = tmp0 + A(p,j  ) * u(i,p,k)
      tmp1 = tmp1 + A(p,j+1) * u(i,p,k)
      tmp2 = tmp2 + A(p,j+2) * u(i,p,k)
      tmp3 = tmp3 + A(p,j+3) * u(i,p,k)
    end do
    v(i,j  ,k) = g * D(i,j  ,k) * tmp0
    v(i,j+1,k) = g * D(i,j+1,k) * tmp1
    v(i,j+2,k) = g * D(i,j+2,k) * tmp2
    v(i,j+3,k) = g * D(i,j+3,k) * tmp3
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
      tmp0 = tmp0 + A(p,j  ) * u(i,p,k)
      tmp1 = tmp1 + A(p,j+1) * u(i,p,k)
      tmp2 = tmp2 + A(p,j+2) * u(i,p,k)
      tmp3 = tmp3 + A(p,j+3) * u(i,p,k)
      tmp4 = tmp4 + A(p,j+4) * u(i,p,k)
    end do
    v(i,j  ,k) = g * D(i,j  ,k) * tmp0
    v(i,j+1,k) = g * D(i,j+1,k) * tmp1
    v(i,j+2,k) = g * D(i,j+2,k) * tmp2
    v(i,j+3,k) = g * D(i,j+3,k) * tmp3
    v(i,j+4,k) = g * D(i,j+4,k) * tmp4
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
      tmp0 = tmp0 + A(p,j  ) * u(i,p,k)
      tmp1 = tmp1 + A(p,j+1) * u(i,p,k)
      tmp2 = tmp2 + A(p,j+2) * u(i,p,k)
      tmp3 = tmp3 + A(p,j+3) * u(i,p,k)
      tmp4 = tmp4 + A(p,j+4) * u(i,p,k)
      tmp5 = tmp5 + A(p,j+5) * u(i,p,k)
    end do
    v(i,j  ,k) = g * D(i,j  ,k) * tmp0
    v(i,j+1,k) = g * D(i,j+1,k) * tmp1
    v(i,j+2,k) = g * D(i,j+2,k) * tmp2
    v(i,j+3,k) = g * D(i,j+3,k) * tmp3
    v(i,j+4,k) = g * D(i,j+4,k) * tmp4
    v(i,j+5,k) = g * D(i,j+5,k) * tmp5
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
      tmp0 = tmp0 + A(p,j  ) * u(i,p,k)
      tmp1 = tmp1 + A(p,j+1) * u(i,p,k)
      tmp2 = tmp2 + A(p,j+2) * u(i,p,k)
      tmp3 = tmp3 + A(p,j+3) * u(i,p,k)
      tmp4 = tmp4 + A(p,j+4) * u(i,p,k)
      tmp5 = tmp5 + A(p,j+5) * u(i,p,k)
      tmp6 = tmp6 + A(p,j+6) * u(i,p,k)
    end do
    v(i,j  ,k) = g * D(i,j  ,k) * tmp0
    v(i,j+1,k) = g * D(i,j+1,k) * tmp1
    v(i,j+2,k) = g * D(i,j+2,k) * tmp2
    v(i,j+3,k) = g * D(i,j+3,k) * tmp3
    v(i,j+4,k) = g * D(i,j+4,k) * tmp4
    v(i,j+5,k) = g * D(i,j+5,k) * tmp5
    v(i,j+6,k) = g * D(i,j+6,k) * tmp6
  end do
  end do

#endif

end subroutine SubOp_2id
