!-----------------------------------------------------------------------------
!> Computes  v = α I⊗Qᵀ⊗I u + β v
!>
!>   * k:  block 1,  unroll 1
!>   * j:  block 1,  unroll 1
!>   * i:  block 1,  unroll 1
!>   * p:  block 1,  unroll 1
!>   * using Intel SIMD directive

subroutine PROC(IxQtxI__,_NQ_)(Q, alpha, beta, u, v)
  !$acc routine vector
  real(_RWP_), intent(in)    :: Q(_NQ_,_NQ_)      !< square matrix
  real(_RWP_), intent(in)    :: alpha             !< factor α
  real(_RWP_), intent(in)    :: beta              !< factor β
  real(_RWP_), intent(in)    :: u(_NQ_,_NQ_,_NQ_) !< operand
  real(_RWP_), intent(inout) :: v(_NQ_,_NQ_,_NQ_) !< result

  real(_RWP_) :: tmp
  integer   :: i, j, k, p

  !$acc loop collapse(3) vector
  do k = 1, _NQ_
  do j = 1, _NQ_
    !DIR$ SIMD VECREMAINDER
    do i = 1, _NQ_
      tmp = 0
      do p = 1, _NQ_
        tmp = tmp + Q(p,j) * u(i,p,k)
      end do
      v(i,j,k) = alpha * tmp + beta * v(i,j,k)
    end do
  end do
  end do

end subroutine PROC(IxQtxI__,_NQ_)
