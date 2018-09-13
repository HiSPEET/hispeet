!-------------------------------------------------------------------------------
!> Performs `v = AxIxI z2` for single element
!>
!>   * k: unrolled and jammed with length of 8 (u8)
!>   * j: joined with i and flattened (f)
!>   * i: joined with j and flattened (f)
!>   * p: simple loop (l)
!>   * explicit remainder handling (r)

subroutine SubOp_3(At, z2, v)
  !$acc routine vector
  real(RNP), intent(in)  :: At(__NA2__,__NA1__)    !< transpose of A
  real(RNP), intent(in)  :: z2(__NA1__**2,__NA2__) !< z2 = IxIxA u
  real(RNP), intent(out) :: v (__NA1__**2,__NA1__) !< result

  real(RNP) :: tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7
  integer   :: ij, k, p

#if __NA1_T8__ > 0

  !$acc loop collapse(2) vector
  do k = 1, __NA1_T8__, 8
  !DIR$ SIMD
  do ij = 1, __NA1__**2
    tmp0 = 0
    tmp1 = 0
    tmp2 = 0
    tmp3 = 0
    tmp4 = 0
    tmp5 = 0
    tmp6 = 0
    tmp7 = 0
    do p = 1, __NA2__
      tmp0 = tmp0 + At(p,k  ) * z2(ij,p)
      tmp1 = tmp1 + At(p,k+1) * z2(ij,p)
      tmp2 = tmp2 + At(p,k+2) * z2(ij,p)
      tmp3 = tmp3 + At(p,k+3) * z2(ij,p)
      tmp4 = tmp4 + At(p,k+4) * z2(ij,p)
      tmp5 = tmp5 + At(p,k+5) * z2(ij,p)
      tmp6 = tmp6 + At(p,k+6) * z2(ij,p)
      tmp7 = tmp7 + At(p,k+7) * z2(ij,p)
    end do
    v(ij,k  ) = tmp0
    v(ij,k+1) = tmp1
    v(ij,k+2) = tmp2
    v(ij,k+3) = tmp3
    v(ij,k+4) = tmp4
    v(ij,k+5) = tmp5
    v(ij,k+6) = tmp6
    v(ij,k+7) = tmp7
  end do
  end do

#endif

#if __NA1__ == __NA1_T8__ + 1

  k =  __NA1_T8__ + 1

  !$acc loop vector
  !DIR$ SIMD
  do ij = 1, __NA1__**2
    tmp0 = 0
    do p = 1, __NA2__
      tmp0 = tmp0 + At(p,k) * z2(ij,p)
    end do
    v(ij,k) = tmp0
  end do

#elif __NA1__ == __NA1_T8__ + 2

  k =  __NA1_T8__ + 1

  !$acc loop vector
  !DIR$ SIMD
  do ij = 1, __NA1__**2
    tmp0 = 0
    tmp1 = 0
    do p = 1, __NA2__
      tmp0 = tmp0 + At(p,k  ) * z2(ij,p)
      tmp1 = tmp1 + At(p,k+1) * z2(ij,p)
    end do
    v(ij,k  ) = tmp0
    v(ij,k+1) = tmp1
  end do

#elif __NA1__ == __NA1_T8__ + 3

  k =  __NA1_T8__ + 1

  !$acc loop vector
  !DIR$ SIMD
  do ij = 1, __NA1__**2
    tmp0 = 0
    tmp1 = 0
    tmp2 = 0
    do p = 1, __NA2__
      tmp0 = tmp0 + At(p,k  ) * z2(ij,p)
      tmp1 = tmp1 + At(p,k+1) * z2(ij,p)
      tmp2 = tmp2 + At(p,k+2) * z2(ij,p)
    end do
    v(ij,k  ) = tmp0
    v(ij,k+1) = tmp1
    v(ij,k+2) = tmp2
  end do

#elif __NA1__ == __NA1_T8__ + 4

  k =  __NA1_T8__ + 1

  !$acc loop vector
  !DIR$ SIMD
  do ij = 1, __NA1__**2
    tmp0 = 0
    tmp1 = 0
    tmp2 = 0
    tmp3 = 0
    do p = 1, __NA2__
      tmp0 = tmp0 + At(p,k  ) * z2(ij,p)
      tmp1 = tmp1 + At(p,k+1) * z2(ij,p)
      tmp2 = tmp2 + At(p,k+2) * z2(ij,p)
      tmp3 = tmp3 + At(p,k+3) * z2(ij,p)
    end do
    v(ij,k  ) = tmp0
    v(ij,k+1) = tmp1
    v(ij,k+2) = tmp2
    v(ij,k+3) = tmp3
  end do

#elif __NA1__ == __NA1_T8__ + 5

  k =  __NA1_T8__ + 1

  !$acc loop vector
  !DIR$ SIMD
  do ij = 1, __NA1__**2
    tmp0 = 0
    tmp1 = 0
    tmp2 = 0
    tmp3 = 0
    tmp4 = 0
    do p = 1, __NA2__
      tmp0 = tmp0 + At(p,k  ) * z2(ij,p)
      tmp1 = tmp1 + At(p,k+1) * z2(ij,p)
      tmp2 = tmp2 + At(p,k+2) * z2(ij,p)
      tmp3 = tmp3 + At(p,k+3) * z2(ij,p)
      tmp4 = tmp4 + At(p,k+4) * z2(ij,p)
    end do
    v(ij,k  ) = tmp0
    v(ij,k+1) = tmp1
    v(ij,k+2) = tmp2
    v(ij,k+3) = tmp3
    v(ij,k+4) = tmp4
  end do

#elif __NA1__ == __NA1_T8__ + 6

  k =  __NA1_T8__ + 1

  !$acc loop vector
  !DIR$ SIMD
  do ij = 1, __NA1__**2
    tmp0 = 0
    tmp1 = 0
    tmp2 = 0
    tmp3 = 0
    tmp4 = 0
    tmp5 = 0
    do p = 1, __NA2__
      tmp0 = tmp0 + At(p,k  ) * z2(ij,p)
      tmp1 = tmp1 + At(p,k+1) * z2(ij,p)
      tmp2 = tmp2 + At(p,k+2) * z2(ij,p)
      tmp3 = tmp3 + At(p,k+3) * z2(ij,p)
      tmp4 = tmp4 + At(p,k+4) * z2(ij,p)
      tmp5 = tmp5 + At(p,k+5) * z2(ij,p)
    end do
    v(ij,k  ) = tmp0
    v(ij,k+1) = tmp1
    v(ij,k+2) = tmp2
    v(ij,k+3) = tmp3
    v(ij,k+4) = tmp4
    v(ij,k+5) = tmp5
  end do

#elif __NA1__ == __NA1_T8__ + 7

  k =  __NA1_T8__ + 1

  !$acc loop vector
  !DIR$ SIMD
  do ij = 1, __NA1__**2
    tmp0 = 0
    tmp1 = 0
    tmp2 = 0
    tmp3 = 0
    tmp4 = 0
    tmp5 = 0
    tmp6 = 0
    do p = 1, __NA2__
      tmp0 = tmp0 + At(p,k  ) * z2(ij,p)
      tmp1 = tmp1 + At(p,k+1) * z2(ij,p)
      tmp2 = tmp2 + At(p,k+2) * z2(ij,p)
      tmp3 = tmp3 + At(p,k+3) * z2(ij,p)
      tmp4 = tmp4 + At(p,k+4) * z2(ij,p)
      tmp5 = tmp5 + At(p,k+5) * z2(ij,p)
      tmp6 = tmp6 + At(p,k+6) * z2(ij,p)
    end do
    v(ij,k  ) = tmp0
    v(ij,k+1) = tmp1
    v(ij,k+2) = tmp2
    v(ij,k+3) = tmp3
    v(ij,k+4) = tmp4
    v(ij,k+5) = tmp5
    v(ij,k+6) = tmp6
  end do

#endif

end subroutine SubOp_3
