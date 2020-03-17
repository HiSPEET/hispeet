!-----------------------------------------------------------------------------
!> Computes  v = α I⊗I⊗Qᵀ u + β v
!>
!>   * k:  block 1,  unroll 1
!>   * j:  block 1,  unroll 1
!>   * i:  block 1,  unroll 1
!>   * p:  block 1,  unroll 1
!>   * using Intel SIMD directive

subroutine PROC(IxIxQt__,_NQ_)(Q, alpha, beta, u, v)
  !$acc routine vector
  real(RNP), intent(in)    :: Q(_NQ_,_NQ_)      !< square matrix
  real(RNP), intent(in)    :: alpha             !< factor α
  real(RNP), intent(in)    :: beta              !< factor β
  real(RNP), intent(in)    :: u(_NQ_,_NQ_,_NQ_) !< operand
  real(RNP), intent(inout) :: v(_NQ_,_NQ_,_NQ_) !< result

  real(RNP) :: tmp
  integer   :: i, j, k, p

  !$acc loop collapse(3) vector
  do k = 1, _NQ_
  do j = 1, _NQ_
    !DIR$ SIMD VECREMAINDER
    do i = 1, _NQ_
      tmp = 0
      do p = 1, _NQ_
        tmp = tmp + Q(p,i) * u(p,j,k)
      end do
      v(i,j,k) = alpha * tmp + beta * v(i,j,k)
    end do
  end do
  end do

end subroutine PROC(IxIxQt__,_NQ_)
