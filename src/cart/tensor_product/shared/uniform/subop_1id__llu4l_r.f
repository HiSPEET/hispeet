!-------------------------------------------------------------------------------
!> Performs `v = g D IxIxA u`
!>
!>   * k: simple loop (l)
!>   * j: simple loop (l)
!>   * i: unrolled and jammed with length of 4 (u4)
!>   * p: simple loop (l)
!>   * explicit remainder handling (r)

subroutine SubOp_1id(g, At, D, u, v)
  !$acc routine vector
  real(RNP), intent(in)  :: g                       !< metric factor
  real(RNP), intent(in)  :: At(__NA__,__NA__)       !< transpose of A
  real(RNP), intent(in)  :: D(__NA__,__NA__,__NA__) !< diagonal operator
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
      tmp0 = tmp0 + At(p,i  ) * u(p,j,k)
      tmp1 = tmp1 + At(p,i+1) * u(p,j,k)
      tmp2 = tmp2 + At(p,i+2) * u(p,j,k)
      tmp3 = tmp3 + At(p,i+3) * u(p,j,k)
    end do
    v(i  ,j,k) = g * D(i  ,j,k) * tmp0
    v(i+1,j,k) = g * D(i+1,j,k) * tmp1
    v(i+2,j,k) = g * D(i+2,j,k) * tmp2
    v(i+3,j,k) = g * D(i+3,j,k) * tmp3
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
      tmp0 = tmp0 + At(p,__NA__) * u(p,j,k)
    end do
    v(__NA__,j,k) = g * D(__NA__,j,k) * tmp0
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
      tmp0 = tmp0 + At(p,__NA__-1) * u(p,j,k)
      tmp1 = tmp1 + At(p,__NA__  ) * u(p,j,k)
    end do
    v(__NA__-1,j,k) = g * D(__NA__-1,j,k) * tmp0
    v(__NA__  ,j,k) = g * D(__NA__  ,j,k) * tmp1
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
      tmp0 = tmp0 + At(p,__NA__-2) * u(p,j,k)
      tmp1 = tmp1 + At(p,__NA__-1) * u(p,j,k)
      tmp2 = tmp2 + At(p,__NA__  ) * u(p,j,k)
    end do
    v(__NA__-2,j,k) = g * D(__NA__-2,j,k) * tmp0
    v(__NA__-1,j,k) = g * D(__NA__-1,j,k) * tmp1
    v(__NA__  ,j,k) = g * D(__NA__  ,j,k) * tmp2
  end do
  end do

#endif

end subroutine SubOp_1id
