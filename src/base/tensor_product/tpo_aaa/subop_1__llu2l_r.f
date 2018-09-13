!-------------------------------------------------------------------------------
!> Performs `z1 = IxIxA u`
!>
!>   * k: simple loop (l)
!>   * j: simple loop (l)
!>   * i: unrolled and jammed with length of 2 (u2)
!>   * p: simple loop (l)
!>   * explicit remainder handling (r)

subroutine SubOp_1(At, u, z1)
  !$acc routine vector
  real(RNP), intent(in)  :: At(__NA2__,__NA1__)         !< transpose of A
  real(RNP), intent(in)  :: u (__NA2__,__NA2__,__NA2__) !< u
  real(RNP), intent(out) :: z1(__NA1__,__NA2__,__NA2__) !< intermediate result

  real(RNP) :: tmp0, tmp1
  integer   :: i, j, k, p

  !$acc loop collapse(3) vector
  do k = 1, __NA2__
  do j = 1, __NA2__
    do i = 1, __NA1_T2__, 2
      tmp0 = 0
      tmp1 = 0
      do p = 1, __NA2__
        tmp0 = tmp0 + At(p,i  ) * u(p,j,k)
        tmp1 = tmp1 + At(p,i+1) * u(p,j,k)
      end do
      z1(i  ,j,k) = tmp0
      z1(i+1,j,k) = tmp1
    end do
  end do
  end do

#if __NA1__ == __NA1_T2__ + 1

  i = __NA1_T2__ + 1

  !$acc loop collapse(2) vector
  do k = 1, __NA2__
    !DIR$ SIMD
    do j = 1, __NA2__
      tmp0 = 0
      do p = 1, __NA2__
        tmp0 = tmp0 + At(p,i) * u(p,j,k)
      end do
      z1(i,j,k) = tmp0
    end do
  end do

#endif

end subroutine SubOp_1
