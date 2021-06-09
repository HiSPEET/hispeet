!-------------------------------------------------------------------------------
!> Computes  v = α I⊗Bᵀ⊗I u + β v
!>
!>   * k:  block 1,  unroll 1
!>   * j:  block 1,  unroll 1
!>   * i:  block 1,  unroll 1
!>   * p:  block 1,  unroll 1
!>   * using Intel SIMD directive

subroutine PROC(IxBtxI__,_NB1_,_NB2_)(na, nc, B, alpha, beta, u, v)
  !$acc routine vector
  integer,   intent(in)    :: na              !< 1st dimension of u,v
  integer,   intent(in)    :: nc              !< 3rd dimension of u,v
  real(_RWP_), intent(in)    :: B(_NB1_,_NB2_)  !< rectangular matrix
  real(_RWP_), intent(in)    :: alpha           !< factor α
  real(_RWP_), intent(in)    :: beta            !< factor β
  real(_RWP_), intent(in)    :: u(na,_NB1_,nc)  !< operand
  real(_RWP_), intent(inout) :: v(na,_NB2_,nc)  !< result

  real(_RWP_) :: tmp
  integer   :: i, j, k, p

  !$acc loop collapse(3) independent vector
  do k = 1, nc
  do j = 1, _NB2_
    !DIR$ SIMD
    do i = 1, na
      tmp = 0
      do p = 1, _NB1_
        tmp = tmp + B(p,j) * u(i,p,k)
      end do
      v(i,j,k) = alpha * tmp + beta * v(i,j,k)
    end do
  end do
  end do

end subroutine PROC(IxBtxI__,_NB1_,_NB2_)
