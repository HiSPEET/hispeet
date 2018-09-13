!-------------------------------------------------------------------------------
!> Performs `z2 = IxAxI z1`
!>
!>   * k: simple loop (l)
!>   * j: unrolled and jammed with length of 2 (u2)
!>   * i: simple loop (l)
!>   * p: simple loop (l)
!>   * explicit remainder handling (r)

subroutine SubOp_2(At, z1, z2)
  !$acc routine vector
  real(RNP), intent(in)  :: At(__NA2__,__NA1__)         !< transpose of A
  real(RNP), intent(in)  :: z1(__NA1__,__NA2__,__NA2__) !< z1 = IxIxA u
  real(RNP), intent(out) :: z2(__NA1__,__NA1__,__NA2__) !< z2

  real(RNP) :: tmp0, tmp1
  integer   :: i, j, k, p

  !$acc loop collapse(3) vector
  do k = 1, __NA2__
  do j = 1, __NA1_T2__, 2
  !DIR$ SIMD
  do i = 1, __NA1__
    tmp0 = 0
    tmp1 = 0
    do p = 1, __NA2__
      tmp0 = tmp0 + At(p,j  ) * z1(i,p,k)
      tmp1 = tmp1 + At(p,j+1) * z1(i,p,k)
    end do
    z2(i,j  ,k) = tmp0
    z2(i,j+1,k) = tmp1
  end do
  end do
  end do

#if __NA1__ == __NA1_T2__ + 1

  !$acc loop collapse(2) vector
  do k = 1, __NA2__
  !DIR$ SIMD
  do i = 1, __NA1__
    tmp0 = 0
    do p = 1, __NA2__
      tmp0 = tmp0 + At(p,__NA1__) * z1(i,p,k)
    end do
    z2(i,__NA1__,k) = tmp0
  end do
  end do

#endif

end subroutine SubOp_2
