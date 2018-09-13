!-------------------------------------------------------------------------------
!> Performs `z2 = IxAxI z1`
!>
!>   * k: simple loop (l)
!>   * j: unrolled and jammed with length of 4 (u4)
!>   * i: simple loop (l)
!>   * p: simple loop (l)
!>   * explicit remainder handling (r)

subroutine SubOp_2(At, z1, z2)
  !$acc routine vector
  real(RNP), intent(in)  :: At(__NA2__,__NA1__)         !< transpose of A
  real(RNP), intent(in)  :: z1(__NA1__,__NA2__,__NA2__) !< z1 = IxIxA u
  real(RNP), intent(out) :: z2(__NA1__,__NA1__,__NA2__) !< z2

  real(RNP) :: tmp0, tmp1, tmp2, tmp3
  integer   :: i, j, k, p

#if __NA1_T4__ > 0

  !$acc loop collapse(3) vector
  do k = 1, __NA2__
  do j = 1, __NA1_T4__, 4
  do i = 1, __NA1__
    tmp0 = 0
    tmp1 = 0
    tmp2 = 0
    tmp3 = 0
    do p = 1, __NA2__
      tmp0 = tmp0 + At(p,j  ) * z1(i,p,k)
      tmp1 = tmp1 + At(p,j+1) * z1(i,p,k)
      tmp2 = tmp2 + At(p,j+2) * z1(i,p,k)
      tmp3 = tmp3 + At(p,j+3) * z1(i,p,k)
    end do
    z2(i,j  ,k) = tmp0
    z2(i,j+1,k) = tmp1
    z2(i,j+2,k) = tmp2
    z2(i,j+3,k) = tmp3
  end do
  end do
  end do

#endif

#if __NA1__ == __NA1_T4__ + 1

  !$acc loop collapse(2) vector
  do k = 1, __NA2__
  do i = 1, __NA1__
    tmp0 = 0
    do p = 1, __NA2__
      tmp0 = tmp0 + At(p,__NA1__) * z1(i,p,k)
    end do
    z2(i,__NA1__,k) = tmp0
  end do
  end do

#elif __NA1__ == __NA1_T4__ + 2

  !$acc loop collapse(2) vector
  do k = 1, __NA2__
  do i = 1, __NA1__
    tmp0 = 0
    tmp1 = 0
    do p = 1, __NA2__
      tmp0 = tmp0 + At(p,__NA1__-1) * z1(i,p,k)
      tmp1 = tmp1 + At(p,__NA1__  ) * z1(i,p,k)
    end do
    z2(i,__NA1__-1,k) = tmp0
    z2(i,__NA1__  ,k) = tmp1
  end do
  end do

#elif __NA1__ == __NA1_T4__ + 3

  !$acc loop collapse(2) vector
  do k = 1, __NA2__
  do i = 1, __NA1__
    tmp0 = 0
    tmp1 = 0
    tmp2 = 0
    do p = 1, __NA2__
      tmp0 = tmp0 + At(p,__NA1__-2) * z1(i,p,k)
      tmp1 = tmp1 + At(p,__NA1__-1) * z1(i,p,k)
      tmp2 = tmp2 + At(p,__NA1__  ) * z1(i,p,k)
    end do
    z2(i,__NA1__-2,k) = tmp0
    z2(i,__NA1__-1,k) = tmp1
    z2(i,__NA1__  ,k) = tmp2
  end do
  end do

#endif

end subroutine SubOp_2
