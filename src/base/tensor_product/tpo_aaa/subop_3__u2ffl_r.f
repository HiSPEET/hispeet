!-------------------------------------------------------------------------------
!> Performs `v = AxIxI z2` for single element
!>
!>   * k: unrolled and jammed with length of 2 (u2)
!>   * j: joined with i and flattened (f)
!>   * i: joined with j and flattened (f)
!>   * p: simple loop (l)
!>   * explicit remainder handling (r)

subroutine SubOp_3(At, z2, v)
  !$acc routine vector
  real(RNP), intent(in)  :: At(__NA2__,__NA1__)    !< transpose of A
  real(RNP), intent(in)  :: z2(__NA1__**2,__NA2__) !< z2 = IxIxA u
  real(RNP), intent(out) :: v (__NA1__**2,__NA1__) !< result

  real(RNP) :: tmp0, tmp1
  integer   :: ij, k, p

  !$acc loop collapse(2) independent vector
  do k = 1, __NA1_T2__, 2
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
  end do

#if __NA1__ == __NA1_T2__ + 1

  !$acc loop independent vector
  !DIR$ SIMD
  do ij = 1, __NA1__**2
    tmp0 = 0
    do p = 1, __NA2__
      tmp0 = tmp0 + At(p,__NA1__) * z2(ij,p)
    end do
    v(ij,__NA1__) = tmp0
  end do

#endif

end subroutine SubOp_3
