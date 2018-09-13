!-------------------------------------------------------------------------------
!> Performs `v = AxIxI z2` for single element
!>
!>   * k: unrolled and jammed with length of 4 (u4)
!>   * j: joined with i and flattened (f)
!>   * i: joined with j and flattened (f)
!>   * p: simple loop (l)
!>   * explicit remainder handling (r)

subroutine SubOp_3(At, z2, v)
  !$acc routine vector
  real(RNP), intent(in)  :: At(__NA2__,__NA1__)    !< transpose of A
  real(RNP), intent(in)  :: z2(__NA1__**2,__NA2__) !< z2 = IxIxA u
  real(RNP), intent(out) :: v (__NA1__**2,__NA1__) !< result

  real(RNP) :: tmp0, tmp1, tmp2, tmp3
  integer   :: ij, k, p

#if __NA1_T4__ > 0

  !$acc loop collapse(2) vector
  do k = 1, __NA1_T4__, 4
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
  end do

#endif

#if __NA1__ == __NA1_T4__ + 1

  !$acc loop vector
  !DIR$ SIMD
  do ij = 1, __NA1__**2
    tmp0 = 0
    do p = 1, __NA2__
      tmp0 = tmp0 + At(p,__NA1__) * z2(ij,p)
    end do
    v(ij,__NA1__) = tmp0
  end do

#elif __NA1__ == __NA1_T4__ + 2

  !$acc loop vector
  !DIR$ SIMD
  do ij = 1, __NA1__**2
    tmp0 = 0
    tmp1 = 0
    do p = 1, __NA2__
      tmp0 = tmp0 + At(p,__NA1__-1) * z2(ij,p)
      tmp1 = tmp1 + At(p,__NA1__  ) * z2(ij,p)
    end do
    v(ij,__NA1__-1) = tmp0
    v(ij,__NA1__  ) = tmp1
  end do

#elif __NA1__ == __NA1_T4__ + 3

  !$acc loop vector
  !DIR$ SIMD
  do ij = 1, __NA1__**2
    tmp0 = 0
    tmp1 = 0
    tmp2 = 0
    do p = 1, __NA2__
      tmp0 = tmp0 + At(p,__NA1__-2) * z2(ij,p)
      tmp1 = tmp1 + At(p,__NA1__-1) * z2(ij,p)
      tmp2 = tmp2 + At(p,__NA1__  ) * z2(ij,p)
    end do
    v(ij,__NA1__-2) = tmp0
    v(ij,__NA1__-1) = tmp1
    v(ij,__NA1__  ) = tmp2
  end do

#endif

end subroutine SubOp_3
