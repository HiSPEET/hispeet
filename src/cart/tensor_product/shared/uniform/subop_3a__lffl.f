!-------------------------------------------------------------------------------
!> Performs `v += g AxIxI u` for single element
!>
!>   * k: simple loop (l)
!>   * j: joined with i and flattened (f)
!>   * i: joined with j and flattened (f)
!>   * p: simple loop (l)

subroutine SubOp_3a(g, At, u, v)
  !$acc routine vector
  real(RNP), intent(in)    :: g                   !< metric factor
  real(RNP), intent(in)    :: At(__NA__,__NA__)   !< transpose of A
  real(RNP), intent(in)    :: u(__NA__**2,__NA__) !< operand
  real(RNP), intent(inout) :: v(__NA__**2,__NA__) !< result

  real(RNP) :: tmp
  integer   :: ij, k, p

  !$acc loop collapse(2) vector
  do k = 1, __NA__
  !DIR$ SIMD VECREMAINDER
  do ij = 1, __NA__**2
    tmp = 0
    do p = 1, __NA__
      tmp = tmp + At(p,k) * u(ij,p)
    end do
    v(ij,k) = v(ij,k) + g * tmp
  end do
  end do

end subroutine SubOp_3a
