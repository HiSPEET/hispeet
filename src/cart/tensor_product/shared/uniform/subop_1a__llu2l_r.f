!-------------------------------------------------------------------------------
!> Performs `v += g IxIxAᵀ u`
!>
!>   * k: simple loop (l)
!>   * j: simple loop (l)
!>   * i: unrolled and jammed with length of 2 (u2)
!>   * p: simple loop (l)

subroutine SubOp_1a(g, A, u, v)
  !$acc routine vector
  real(RNP), intent(in)    :: g                       !< metric factor
  real(RNP), intent(in)    :: A(__NA__,__NA__)        !< A
  real(RNP), intent(in)    :: u(__NA__,__NA__,__NA__) !< operand
  real(RNP), intent(inout) :: v(__NA__,__NA__,__NA__) !< result

  real(RNP) :: tmp0, tmp1
  integer   :: i, j, k, p

  !$acc loop collapse(3) vector
  do k = 1, __NA__
  do j = 1, __NA__
  do i = 1, __NA_T2__, 2
    tmp0 = 0
    tmp1 = 0
    do p = 1, __NA__
      tmp0 = tmp0 + A(p,i  ) * u(p,j,k)
      tmp1 = tmp1 + A(p,i+1) * u(p,j,k)
    end do
    v(i  ,j,k) = v(i  ,j,k) + g * tmp0
    v(i+1,j,k) = v(i+1,j,k) + g * tmp1
  end do
  end do
  end do

#if __NA__ == __NA_T2__ + 1
  !$acc loop collapse(2)  vector
  do k = 1, __NA__
    !DIR$ SIMD
    do j = 1, __NA__
      tmp0 = 0
      do p = 1, __NA__
        tmp0 = tmp0 + A(p,__NA__) * u(p,j,k)
      end do
      v(__NA__,j,k) = v(__NA__,j,k) + g * tmp0
    end do
  end do
#endif

end subroutine SubOp_1a
