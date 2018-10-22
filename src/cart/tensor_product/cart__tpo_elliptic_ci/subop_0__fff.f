!-------------------------------------------------------------------------------
!> Projection to element mass matrix and computation of the Helmholtz term
!>
!>   * k:  joined with i,j and flattened  (f)
!>   * j:  joined with i,k and flattened  (f)
!>   * i:  joined with j,k and flattened  (f)

subroutine SubOp_0(lambda, M, u, M_u, v)
  !$acc routine vector
  real(RNP), intent(in)  :: lambda
  real(RNP), intent(in)  :: M(__NA__**3)
  real(RNP), intent(in)  :: u(__NA__**3)
  real(RNP), intent(out) :: M_u(__NA__**3)
  real(RNP), intent(out) :: v(__NA__**3)

  integer :: l

  !$acc loop vector
  !DIR$ SIMD
  do l = 1, __NA__**3
    M_u(l) = M(l) * u(l)
    v(l) = lambda * M_u(l)
  end do

end subroutine SubOp_0
