!-------------------------------------------------------------------------------
!> Performs `v = g D AᵀxIxI u` for single element
!>
!>   * k: simple loop (l)
!>   * j: joined with i and flattened (f)
!>   * i: joined with j and flattened (f)
!>   * p: simple loop (l)

subroutine SubOp_3id(g, A, D, u, v)
  !$acc routine vector
  real(RNP), intent(in)  :: g                   !< metric factor
  real(RNP), intent(in)  :: A(__NA__,__NA__)    !< A
  real(RNP), intent(in)  :: D(__NA__**2,__NA__) !< diagonal operator
  real(RNP), intent(in)  :: u(__NA__**2,__NA__) !< operand
  real(RNP), intent(out) :: v(__NA__**2,__NA__) !< result

  real(RNP) :: tmp
  integer   :: ij, k, p

  !$acc loop collapse(2) vector
  do k = 1, __NA__
  !DIR$ SIMD VECREMAINDER
  do ij = 1, __NA__**2
    tmp = 0
    do p = 1, __NA__
      tmp = tmp + A(p,k) * u(ij,p)
    end do
    v(ij,k) = g * D(ij,k) * tmp
  end do
  end do

end subroutine SubOp_3id
