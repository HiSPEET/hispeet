!-------------------------------------------------------------------------------
!> Computes  v = α Cᵀ⊗I⊗I u + β v
!>
!>   * k:  block 1,  unroll 1
!>   * j:  joined with i
!>   * i:  joined with j
!>   * p:  block 1,  unroll 1
!>   * using Intel SIMD directive

subroutine PROC(CtxIxI__,_NC1_,_NC2_)(na_nb, C, alpha, beta, u, v)
  !$acc routine vector
  integer,   intent(in)    :: na_nb           !< size of dimensions 1+2 of u,v
  real(RNP), intent(in)    :: C(_NC1_,_NC2_)  !< rectangular matrix
  real(RNP), intent(in)    :: alpha           !< factor α
  real(RNP), intent(in)    :: beta            !< factor β
  real(RNP), intent(in)    :: u(na_nb,_NC1_)  !< operand
  real(RNP), intent(inout) :: v(na_nb,_NC2_)  !< result

  real(RNP) :: tmp
  integer   :: ij, k, p

  !$acc loop collapse(2) vector
  do k = 1, _NC2_
    !DIR$ SIMD
    do ij = 1, na_nb
      tmp = 0
      do p = 1, _NC1_
        tmp = tmp + C(p,k) * u(ij,p)
      end do
      v(ij,k) = alpha * tmp + beta * v(ij,k)
    end do
  end do

end subroutine PROC(CtxIxI__,_NC1_,_NC2_)
