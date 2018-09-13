!-------------------------------------------------------------------------------
!> Performs `z1 = IxIxA u`
!>
!>   * k: simple loop (l)
!>   * j: unrolled and jammed with length of 4 (u4)
!>   * i: unrolled and jammed with length of 2 (u2)
!>   * p: simple loop (l)
!>   * explicit remainder handling (r)

subroutine SubOp_1(At, u, z1)
  !$acc routine vector
  real(RNP), intent(in)  :: At(__NA2__,__NA1__)         !< transpose of A
  real(RNP), intent(in)  :: u (__NA2__,__NA2__,__NA2__) !< u
  real(RNP), intent(out) :: z1(__NA1__,__NA2__,__NA2__) !< intermediate result

  real(RNP) :: tmp00, tmp10
  real(RNP) :: tmp01, tmp11
  real(RNP) :: tmp02, tmp12
  real(RNP) :: tmp03, tmp13
  integer   :: i, j, k, p

#if __NA2_T4__ > 0

  !$acc loop collapse(3) vector
  do k = 1, __NA2__
    do j = 1, __NA2_T4__, 4
    do i = 1, __NA1_T2__, 2
      tmp00 = 0
      tmp10 = 0
      tmp01 = 0
      tmp11 = 0
      tmp02 = 0
      tmp12 = 0
      tmp03 = 0
      tmp13 = 0
      do p = 1, __NA2__
        tmp00 = tmp00 + At(p,i  ) * u(p,j  ,k)
        tmp10 = tmp10 + At(p,i+1) * u(p,j  ,k)
        tmp01 = tmp01 + At(p,i  ) * u(p,j+1,k)
        tmp11 = tmp11 + At(p,i+1) * u(p,j+1,k)
        tmp02 = tmp02 + At(p,i  ) * u(p,j+2,k)
        tmp12 = tmp12 + At(p,i+1) * u(p,j+2,k)
        tmp03 = tmp03 + At(p,i  ) * u(p,j+3,k)
        tmp13 = tmp13 + At(p,i+1) * u(p,j+3,k)
      end do
      z1(i  ,j  ,k) = tmp00
      z1(i+1,j  ,k) = tmp10
      z1(i  ,j+1,k) = tmp01
      z1(i+1,j+1,k) = tmp11
      z1(i  ,j+2,k) = tmp02
      z1(i+1,j+2,k) = tmp12
      z1(i  ,j+3,k) = tmp03
      z1(i+1,j+3,k) = tmp13
    end do
  end do
  end do

#endif

#if __NA2__ == __NA2_T4__ + 1

  ! remainder: i = 1:na1_t2, j = na2 -------------------------------------------

  j = __NA2_T4__ + 1

  !$acc loop collapse(2) vector
  do k = 1, __NA2__
    do i = 1, __NA1_T2__, 2
      tmp00 = 0
      tmp10 = 0
      do p = 1, __NA2__
        tmp00 = tmp00 + At(p,i  ) * u(p,j  ,k)
        tmp10 = tmp10 + At(p,i+1) * u(p,j  ,k)
      end do
      z1(i  ,j  ,k) = tmp00
      z1(i+1,j  ,k) = tmp10
    end do
  end do

#elif __NA2__ == __NA2_T4__ + 2

  ! remainder: i = 1:na1_t2, j = na2-1:na2 -------------------------------------

  j = __NA2_T4__ + 1

  !$acc loop collapse(2) vector
  do k = 1, __NA2__
    do i = 1, __NA1_T2__, 2
      tmp00 = 0
      tmp10 = 0
      tmp01 = 0
      tmp11 = 0
      do p = 1, __NA2__
        tmp00 = tmp00 + At(p,i  ) * u(p,j  ,k)
        tmp10 = tmp10 + At(p,i+1) * u(p,j  ,k)
        tmp01 = tmp01 + At(p,i  ) * u(p,j+1,k)
        tmp11 = tmp11 + At(p,i+1) * u(p,j+1,k)
      end do
      z1(i  ,j  ,k) = tmp00
      z1(i+1,j  ,k) = tmp10
      z1(i  ,j+1,k) = tmp01
      z1(i+1,j+1,k) = tmp11
    end do
  end do

#elif __NA2__ == __NA2_T4__ + 3

  ! remainder: i = 1:na1_t2, j = na2-2:na2 -------------------------------------

  j = __NA2_T4__ + 1

  !$acc loop collapse(2) vector
  do k = 1, __NA2__
    do i = 1, __NA1_T2__, 2
      tmp00 = 0
      tmp10 = 0
      tmp01 = 0
      tmp11 = 0
      tmp02 = 0
      tmp12 = 0
      do p = 1, __NA2__
        tmp00 = tmp00 + At(p,i  ) * u(p,j  ,k)
        tmp10 = tmp10 + At(p,i+1) * u(p,j  ,k)
        tmp01 = tmp01 + At(p,i  ) * u(p,j+1,k)
        tmp11 = tmp11 + At(p,i+1) * u(p,j+1,k)
        tmp02 = tmp02 + At(p,i  ) * u(p,j+2,k)
        tmp12 = tmp12 + At(p,i+1) * u(p,j+2,k)
      end do
      z1(i  ,j  ,k) = tmp00
      z1(i+1,j  ,k) = tmp10
      z1(i  ,j+1,k) = tmp01
      z1(i+1,j+1,k) = tmp11
      z1(i  ,j+2,k) = tmp02
      z1(i+1,j+2,k) = tmp12
    end do
  end do

#endif

#if __NA1__ == __NA1_T2__ + 1

  ! remainder: i = na, j = 1:na2 -----------------------------------------------

  i = __NA1__

  !$acc loop collapse(2) vector
  do k = 1, __NA2__
    !DIR$ SIMD
    do j = 1, __NA2__
      tmp00 = 0
      do p = 1, __NA2__
        tmp00 = tmp00 + At(p,i) * u(p,j,k)
      end do
      z1(i,j,k) = tmp00
    end do
  end do

#endif

end subroutine SubOp_1
