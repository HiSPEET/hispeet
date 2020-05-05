!-------------------------------------------------------------------------------
!> Computes  v = α I⊗I⊗Aᵀ u + β v
!>
!>   * k:  block 1,  unroll 1
!>   * j:  block 1,  unroll 1
!>   * i:  block 1,  unroll 1
!>   * p:  block 1,  unroll 1
!>   * using Intel SIMD directive

subroutine SubOp_1(A, u, v)
  !$acc routine vector
  real(RNP), intent(in)    :: A(__NA2__,__NA1__)         !< rectangular matrix
  real(RNP), intent(in)    :: alpha                      !< factor α
  real(RNP), intent(in)    :: beta                       !< factor β
  real(RNP), intent(in)    :: u(__NA2__,__NA2__,__NA2__) !< operand
  real(RNP), intent(inout) :: v(__NA1__,__NA2__,__NA2__) !< result

  real(RNP) :: tmp
  integer   :: i, j, k, p

  !$acc loop collapse(3) independent vector
  do k = 1, __NA2__
  do j = 1, __NA2__
    !DIR$ SIMD VECREMAINDER
    do i = 1, __NA1__
      tmp = 0
      do p = 1, __NA2__
        tmp = tmp + A(p,i) * u(p,j,k)
      end do
      v(i,j,k) = alpha * tmp + beta * v(i,j,k)
    end do
  end do
  end do

end subroutine SubOp_1
