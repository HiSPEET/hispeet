!-------------------------------------------------------------------------------
!> Computes  v = α I⊗I⊗Aᵀ u + β v
!>
!>   * k:  block 1,  unroll 1
!>   * j:  block 1,  unroll 1
!>   * i:  block 1,  unroll 1
!>   * p:  block 1,  unroll 1
!>   * using Intel SIMD directive

subroutine PROC(IxIxAt__,_NA1_,_NA2_)(nb, nc, A, alpha, beta, u, v)
  !$acc routine vector
  integer,   intent(in)    :: nb              !< 2nd dimension of u,v
  integer,   intent(in)    :: nc              !< 3rd dimension of u,v
  real(RNP), intent(in)    :: A(_NA1_,_NA2_)  !< rectangular matrix
  real(RNP), intent(in)    :: alpha           !< factor α
  real(RNP), intent(in)    :: beta            !< factor β
  real(RNP), intent(in)    :: u(_NA1_,nb,nc)  !< operand
  real(RNP), intent(inout) :: v(_NA2_,nb,nc)  !< result

  real(RNP) :: tmp
  integer   :: i, j, k, p

  !$acc loop collapse(3) independent vector
  do k = 1, nc
  do j = 1, nb
    !DIR$ SIMD
    do i = 1, _NA2_
      tmp = 0
      do p = 1, _NA1_
        tmp = tmp + A(p,i) * u(p,j,k)
      end do
      v(i,j,k) = alpha * tmp + beta * v(i,j,k)
    end do
  end do
  end do

end subroutine PROC(IxIxAt__,_NA1_,_NA2_)
