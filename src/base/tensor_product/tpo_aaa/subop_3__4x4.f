!-------------------------------------------------------------------------------
!> Performs `v = AxIxI z2` for single element and na1 = na2 = 4
!>
!>   * k: unrolled and jammed with length of 4 (u4)
!>   * j: joined with i and flattened (f)
!>   * i: joined with j and flattened (f)
!>   * p: blocked with length of 4 (b4)
!>   * explicit remainder handling (r)

subroutine SubOp_3(At, z2, v)
  !$acc routine vector
  real(RNP), intent(in)  :: At(4,4)    !< transpose of A
  real(RNP), intent(in)  :: z2(4**2,4) !< z2 = IxIxA u
  real(RNP), intent(out) :: v (4**2,4) !< result

  real(RNP) :: tmp0, tmp1, tmp2, tmp3
  integer   :: ij, p

  !$acc loop vector
  !DIR$ SIMD
  do ij = 1, 4**2
    tmp0 = 0
    tmp1 = 0
    tmp2 = 0
    tmp3 = 0
    do p = 1, 4
      tmp0 = tmp0 + At(p,1) * z2(ij,p)
      tmp1 = tmp1 + At(p,2) * z2(ij,p)
      tmp2 = tmp2 + At(p,3) * z2(ij,p)
      tmp3 = tmp3 + At(p,4) * z2(ij,p)
    end do
    v(ij,1) = tmp0
    v(ij,2) = tmp1
    v(ij,3) = tmp2
    v(ij,4) = tmp3
  end do

end subroutine SubOp_3
