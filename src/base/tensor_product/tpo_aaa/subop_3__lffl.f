!-------------------------------------------------------------------------------
!> Performs `v = AxIxI z2` for single element
!>
!>   * k: simple loop (l)
!>   * j: joined with i and flattened (f)
!>   * i: joined with j and flattened (f)
!>   * p: simple loop (l)

subroutine SubOp_3(At, z2, v)
  !$acc routine vector
  real(RNP), intent(in)  :: At(__NA2__,__NA1__)     !< transpose of A
  real(RNP), intent(in)  :: z2(__NA1__**2,__NA2__)  !< z2 = IxIxA u
  real(RNP), intent(out) :: v (__NA1__**2,__NA1__)  !< result

  real(RNP) :: tmp
  integer   :: ij, k, p

  !$acc loop collapse(2) vector
  do k = 1, __NA1__
    !DIR$ SIMD
    do ij = 1, __NA1__**2
      tmp = 0
      do p = 1, __NA2__
        tmp = tmp + At(p,k) * z2(ij,p)
      end do
      v(ij,k) = tmp
    end do
  end do

end subroutine SubOp_3