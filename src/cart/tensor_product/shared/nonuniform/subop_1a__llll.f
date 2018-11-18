!-------------------------------------------------------------------------------
!> Performs `v += g IxIxAᵀ u`
!>
!>   * k: simple loop (l)
!>   * j: simple loop (l)
!>   * i: simple loop (l)
!>   * p: simple loop (l)
!>   * using Intel SIMD directive

subroutine SubOp_1a(g, A, u, v)
  !$acc routine vector
  real(RNP), intent(in)    :: g                       !< metric factor
  real(RNP), intent(in)    :: A(__NA__,__NA__)        !< transpose of A
  real(RNP), intent(in)    :: u(__NA__,__NB__,__NC__) !< operand
  real(RNP), intent(inout) :: v(__NA__,__NB__,__NC__) !< result

  real(RNP) :: tmp
  integer   :: i, j, k, p

  !$acc loop collapse(3) vector
  do k = 1, __NC__
  do j = 1, __NB__
    !DIR$ SIMD VECREMAINDER
    do i = 1, __NA__
      tmp = 0
      do p = 1, __NA__
        tmp = tmp + A(p,i) * u(p,j,k)
      end do
      v(i,j,k) = v(i,j,k) + g * tmp
    end do
  end do
  end do

end subroutine SubOp_1a
