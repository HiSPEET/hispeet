!-------------------------------------------------------------------------------
!> Performs `v = D (I x I x A^T) u`
!>
!>   * k: simple loop (l)
!>   * j: unrolled and jammed with length of 2 (u2)
!>   * i: unrolled and jammed with length of 2 (u2)
!>   * p: simple loop (l)
!>   * explicit remainder handling (r)

subroutine SubOp_1d(A, D, u, v)
  !$acc routine vector
  real(RNP), intent(in)  :: D(__NA__,__NA__,__NA__) !< 3D diagonal operator
  real(RNP), intent(in)  :: A(__NA__,__NA__)        !< A
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
      tmp0 = tmp0 + A(p,i  ) * u(p,j  ,k)
      tmp1 = tmp1 + A(p,i+1) * u(p,j  ,k)
      tmp2 = tmp2 + A(p,i  ) * u(p,j+1,k)
      tmp3 = tmp3 + A(p,i+1) * u(p,j+1,k)
    end do
    v(i  ,j  ,k) = tmp0 * D(i  ,j  ,k)
    v(i+1,j  ,k) = tmp1 * D(i+1,j  ,k)
    v(i  ,j+1,k) = tmp2 * D(i  ,j+1,k)
    v(i+1,j+1,k) = tmp3 * D(i+1,j+1,k)
  end do
  end do
  end do

#if __NA__ == __NA_T2__ + 1

  !$acc loop independent
  do k = 1, __NA__

    ! remainder: i = 1:np-1, j = np ..........................................

    !$acc loop vector
    do i = 1, __NA_T2__, 2
      tmp0 = 0
      tmp1 = 0
      do p = 1, __NA__
        tmp0 = tmp0 + A(p,i  ) * u(p,__NA__,k)
        tmp1 = tmp1 + A(p,i+1) * u(p,__NA__,k)
      end do
      v(i  ,__NA__,k) = tmp0 * D(i  ,__NA__,k)
      v(i+1,__NA__,k) = tmp1 * D(i+1,__NA__,k)
    end do

    ! remainder: i = np, j = 1:np ............................................

    !$acc loop vector
    do j = 1, __NA__
      tmp0 = 0
      do p = 1, __NA__
        tmp0 = tmp0 + A(p,__NA__) * u(p,j,k)
      end do
      v(__NA__,j,k) = tmp0 * D(__NA__,j,k)
    end do

  end do

#endif

end subroutine SubOp_1d
