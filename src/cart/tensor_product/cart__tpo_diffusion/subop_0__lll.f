!-------------------------------------------------------------------------------
!> Projection to element mass matrix and computation of the Helmholtz term
!>
!>   * k:  simple loop (l)
!>   * j:  simple loop (l)
!>   * i:  simple loop (l)

subroutine SubOp_0(lambda, M, u, M_u, v)
  !$acc routine vector
  real(RNP), intent(in)  :: lambda
  real(RNP), intent(in)  :: M(__NA__,__NA__,__NA__)
  real(RNP), intent(in)  :: u(__NA__,__NA__,__NA__)
  real(RNP), intent(out) :: M_u(__NA__,__NA__,__NA__)
  real(RNP), intent(out) :: v(__NA__,__NA__,__NA__)

  integer :: i, j, k

  !$acc loop collapse(3) vector
  do k = 1, __NA__
  do j = 1, __NA__
  !DIR$ SIMD
  do i = 1, __NA__
    M_u(i,j,k) = M(i,j,k) * u(i,j,k)
    v(i,j,k) = lambda * M_u(i,j,k)
  end do
  end do
  end do

end subroutine SubOp_0
