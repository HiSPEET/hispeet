!-------------------------------------------------------------------------------
!> Performs `v = (d₀IxIxI + d₁IxIxV₁ + d₂IxIxV₂ + d₃IxIxV₃)⁻¹ (IxIxA') u`
!>
!>   * k: simple loop (l)
!>   * j: simple loop (l)
!>   * i: simple loop (l)
!>   * p: simple loop (l)
!>   * using Intel SIMD directive

subroutine SubOp_1d(d0, d1, d2, d3, V1, V2, V3, A, u, v)
  !$acc routine vector
  real(RNP), intent(in)  :: d0                       !< λ Δx₁ Δx₂ Δx₃
  real(RNP), intent(in)  :: d1                       !< ν Δx₂ Δx₃ / Δx₁
  real(RNP), intent(in)  :: d2                       !< ν Δx₃ Δx₁ / Δx₂
  real(RNP), intent(in)  :: d3                       !< ν Δx₁ Δx₂ / Δx₃
  real(RNP), intent(in)  :: V1(__NA__)               !< 1D eigenvalues in dir 1
  real(RNP), intent(in)  :: V2(__NA__)               !< 1D eigenvalues in dir 2
  real(RNP), intent(in)  :: V3(__NA__)               !< 1D eigenvalues in dir 3
  real(RNP), intent(in)  :: A(__NA__,__NA__)         !< A
  real(RNP), intent(in)  :: u(__NA__,__NA__,__NA__)  !< operand
  real(RNP), intent(out) :: v(__NA__,__NA__,__NA__)  !< result

  real(RNP) :: tmp
  integer   :: i, j, k, p

  !$acc loop collapse(3) vector
  do k = 1, __NA__
  do j = 1, __NA__
    !DIR$ SIMD VECREMAINDER
    do i = 1, __NA__
      tmp = 0
      do p = 1, __NA__
        tmp = tmp + A(p,i) * u(p,j,k)
      end do
      v(i,j,k) = tmp / (d0 + d1*V1(i) + d2*V2(j) + d3*V3(k))
    end do
  end do
  end do

end subroutine SubOp_1d
