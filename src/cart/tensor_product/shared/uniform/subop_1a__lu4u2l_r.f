!-------------------------------------------------------------------------------
!> Performs `v += g1 IxIxA u`
!>
!>   * k: simple loop (l)
!>   * j: unrolled and jammed with length of 4 (u4)
!>   * i: unrolled and jammed with length of 2 (u2)
!>   * p: simple loop (l)
!>   * explicit remainder handling (r)

subroutine SubOp_1a(g1, At, u, v)
  !$acc routine vector
  real(RNP), intent(in)    :: g1                      !< metric factor
  real(RNP), intent(in)    :: At(__NA__,__NA__)       !< transpose of A
  real(RNP), intent(in)    :: u(__NA__,__NA__,__NA__) !< operand
  real(RNP), intent(inout) :: v(__NA__,__NA__,__NA__) !< result

  real(RNP) :: tmp00, tmp10
  real(RNP) :: tmp01, tmp11
  real(RNP) :: tmp02, tmp12
  real(RNP) :: tmp03, tmp13
  integer   :: i, j, k, p

#if __NA_T4__ > 0

  !$acc loop collapse(3) vector
  do k = 1, __NA__
  do j = 1, __NA_T4__, 4
  do i = 1, __NA_T2__, 2
    tmp00 = 0
    tmp10 = 0
    tmp01 = 0
    tmp11 = 0
    tmp02 = 0
    tmp12 = 0
    tmp03 = 0
    tmp13 = 0
    do p = 1, __NA__
      tmp00 = tmp00 + At(p,i  ) * u(p,j  ,k)
      tmp10 = tmp10 + At(p,i+1) * u(p,j  ,k)
      tmp01 = tmp01 + At(p,i  ) * u(p,j+1,k)
      tmp11 = tmp11 + At(p,i+1) * u(p,j+1,k)
      tmp02 = tmp02 + At(p,i  ) * u(p,j+2,k)
      tmp12 = tmp12 + At(p,i+1) * u(p,j+2,k)
      tmp03 = tmp03 + At(p,i  ) * u(p,j+3,k)
      tmp13 = tmp13 + At(p,i+1) * u(p,j+3,k)
    end do
    v(i  ,j  ,k) = v(i  ,j  ,k) + g1 * tmp00
    v(i+1,j  ,k) = v(i+1,j  ,k) + g1 * tmp10
    v(i  ,j+1,k) = v(i  ,j+1,k) + g1 * tmp01
    v(i+1,j+1,k) = v(i+1,j+1,k) + g1 * tmp11
    v(i  ,j+2,k) = v(i  ,j+2,k) + g1 * tmp02
    v(i+1,j+2,k) = v(i+1,j+2,k) + g1 * tmp12
    v(i  ,j+3,k) = v(i  ,j+3,k) + g1 * tmp03
    v(i+1,j+3,k) = v(i+1,j+3,k) + g1 * tmp13
  end do
  end do
  end do

#endif

#if __NA__ == __NA_T4__ + 1

  ! remainder: i = 1:na-1, j = na ............................................

  j = __NA__

  !$acc loop collapse(2) vector
  do k = 1, __NA__
  do i = 1, __NA_T2__, 2
    tmp00 = 0
    tmp10 = 0
    do p = 1, __NA__
      tmp00 = tmp00 + At(p,i  ) * u(p,j,k)
      tmp10 = tmp10 + At(p,i+1) * u(p,j,k)
    end do
    v(i  ,j,k) = v(i  ,j,k) + g1 * tmp00
    v(i+1,j,k) = v(i+1,j,k) + g1 * tmp10
  end do
  end do

  ! remainder: i = na, j = 1:na ..............................................

  i = __NA__

  !$acc loop collapse(2) vector
  do k = 1, __NA__
  do j = 1, __NA__
    tmp00 = 0
    do p = 1, __NA__
      tmp00 = tmp00 + At(p,i) * u(p,j,k)
    end do
    v(i,j,k) = v(i,j,k) + g1 * tmp00
  end do
  end do

#elif __NA__ == __NA_T4__ + 2

  ! remainder: i = 1:na, j = na-1:na .........................................

  j = __NA__ - 1

  !$acc loop collapse(2) vector
  do k = 1, __NA__
  do i = 1, __NA_T2__, 2
    tmp00 = 0
    tmp10 = 0
    tmp01 = 0
    tmp11 = 0
    do p = 1, __NA__
      tmp00 = tmp00 + At(p,i  ) * u(p,j  ,k)
      tmp10 = tmp10 + At(p,i+1) * u(p,j  ,k)
      tmp01 = tmp01 + At(p,i  ) * u(p,j+1,k)
      tmp11 = tmp11 + At(p,i+1) * u(p,j+1,k)
    end do
    v(i  ,j  ,k) = v(i  ,j  ,k) + g1 * tmp00
    v(i+1,j  ,k) = v(i+1,j  ,k) + g1 * tmp10
    v(i  ,j+1,k) = v(i  ,j+1,k) + g1 * tmp01
    v(i+1,j+1,k) = v(i+1,j+1,k) + g1 * tmp11
  end do
  end do

#elif __NA__ == __NA_T4__ + 3

  ! remainder: i = 1:na-1, j = na-2:na .......................................

  j = __NA__ - 2

  !$acc loop collapse(2) vector
  do k = 1, __NA__
  do i = 1, __NA_T2__, 2
     tmp00 = 0
     tmp10 = 0
     tmp01 = 0
     tmp11 = 0
     tmp02 = 0
     tmp12 = 0
     do p = 1, __NA__
       tmp00 = tmp00 + At(p,i  ) * u(p,j  ,k)
       tmp10 = tmp10 + At(p,i+1) * u(p,j  ,k)
       tmp01 = tmp01 + At(p,i  ) * u(p,j+1,k)
       tmp11 = tmp11 + At(p,i+1) * u(p,j+1,k)
       tmp02 = tmp02 + At(p,i  ) * u(p,j+2,k)
       tmp12 = tmp12 + At(p,i+1) * u(p,j+2,k)
     end do
     v(i  ,j  ,k) = v(i  ,j  ,k) + g1 * tmp00
     v(i+1,j  ,k) = v(i+1,j  ,k) + g1 * tmp10
     v(i  ,j+1,k) = v(i  ,j+1,k) + g1 * tmp01
     v(i+1,j+1,k) = v(i+1,j+1,k) + g1 * tmp11
     v(i  ,j+2,k) = v(i  ,j+2,k) + g1 * tmp02
     v(i+1,j+2,k) = v(i+1,j+2,k) + g1 * tmp12
  end do
  end do

  ! remainder: i = na, j = 1:na ..............................................

  i = __NA__

  !$acc loop collapse(2) vector
  do k = 1, __NA__
  do j = 1, __NA__
    tmp00 = 0
    do p = 1, __NA__
      tmp00 = tmp00 + At(p,i) * u(p,j,k)
    end do
    v(i,j,k) = v(i,j,k) + g1 * tmp00
  end do
  end do

#endif

end subroutine SubOp_1a
