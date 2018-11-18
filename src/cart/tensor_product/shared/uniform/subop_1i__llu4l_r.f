!-------------------------------------------------------------------------------
!> Performs `v = g IxIxAᵀ u`
!>
!>   * k: simple loop (l)
!>   * j: simple loop (l)
!>   * i: unrolled and jammed with length of 4 (u4)
!>   * p: simple loop (l)
!>   * explicit remainder handling (r)

subroutine SubOp_1i(g, A, u, v)
  !$acc routine vector
  real(RNP), intent(in)  :: g                       !< metric factor
  real(RNP), intent(in)  :: A(__NA__,__NA__)        !< A
  real(RNP), intent(in)  :: u(__NA__,__NA__,__NA__) !< operand
  real(RNP), intent(out) :: v(__NA__,__NA__,__NA__) !< result

  real(RNP) :: tmp0, tmp1, tmp2, tmp3
  integer   :: i, j, k, p

#if __NA_T4__ > 0

  !$acc loop collapse(3) vector
  do k = 1, __NA__
  do j = 1, __NA__
  do i = 1, __NA_T4__, 4
    tmp0 = 0
    tmp1 = 0
    tmp2 = 0
    tmp3 = 0
    do p = 1, __NA__
      tmp0 = tmp0 + A(p,i  ) * u(p,j,k)
      tmp1 = tmp1 + A(p,i+1) * u(p,j,k)
      tmp2 = tmp2 + A(p,i+2) * u(p,j,k)
      tmp3 = tmp3 + A(p,i+3) * u(p,j,k)
    end do
    v(i  ,j,k) = g * tmp0
    v(i+1,j,k) = g * tmp1
    v(i+2,j,k) = g * tmp2
    v(i+3,j,k) = g * tmp3
  end do
  end do
  end do

#endif

#if __NA__ == __NA_T4__ + 1

  ! remainder: i = np ........................................................

  !$acc loop collapse(2) vector
  do k = 1, __NA__
  do j = 1, __NA__
    tmp0 = 0
    do p = 1, __NA__
      tmp0 = tmp0 + A(p,__NA__) * u(p,j,k)
    end do
    v(__NA__,j,k) = g * tmp0
  end do
  end do

#elif __NA__ == __NA_T4__ + 2

  ! remainder: i = np-1:np ...................................................

  !$acc loop collapse(2) vector
  do k = 1, __NA__
  do j = 1, __NA__
    tmp0 = 0
    tmp1 = 0
    do p = 1, __NA__
      tmp0 = tmp0 + A(p,__NA__-1) * u(p,j,k)
      tmp1 = tmp1 + A(p,__NA__  ) * u(p,j,k)
    end do
    v(__NA__-1,j,k) = g * tmp0
    v(__NA__  ,j,k) = g * tmp1
  end do
  end do

#elif __NA__ == __NA_T4__ + 3

  ! remainder: i = np-2:np ...................................................

  !$acc loop collapse(2) vector
  do k = 1, __NA__
  do j = 1, __NA__
    tmp0 = 0
    tmp1 = 0
    tmp2 = 0
    do p = 1, __NA__
      tmp0 = tmp0 + A(p,__NA__-2) * u(p,j,k)
      tmp1 = tmp1 + A(p,__NA__-1) * u(p,j,k)
      tmp2 = tmp2 + A(p,__NA__  ) * u(p,j,k)
    end do
    v(__NA__-2,j,k) = g * tmp0
    v(__NA__-1,j,k) = g * tmp1
    v(__NA__  ,j,k) = g * tmp2
  end do
  end do

#endif

end subroutine SubOp_1i
