!-------------------------------------------------------------------------------
!> Computes  v = α Cᵀ⊗I⊗I u + β v
!>
!>   * k:  block 1,  unroll 2
!>   * j:  joined with i
!>   * i:  joined with j
!>   * p:  block 1,  unroll 1
!>   * using Intel SIMD directive

#define _NC2_T2_ (_NC2_ / 2) * 2

subroutine PROC(CtxIxI__,_NC1_,_NC2_)(na_nb, C, alpha, beta, u, v)
  !$acc routine vector
  integer,   intent(in)    :: na_nb           !< size of dimensions 1+2 of u,v
  real(RNP), intent(in)    :: C(_NC1_,_NC2_)  !< rectangular matrix
  real(RNP), intent(in)    :: alpha           !< factor α
  real(RNP), intent(in)    :: beta            !< factor β
  real(RNP), intent(in)    :: u(na_nb,_NC1_)  !< operand
  real(RNP), intent(inout) :: v(na_nb,_NC2_)  !< result

  real(RNP) :: tmp0, tmp1
  integer   :: ij, k, p

  do k = 1, _NC2_T2_, 2
  !DIR$ SIMD
  do ij = 1, na_nb
    tmp0 = 0
    tmp1 = 0
    do p = 1, _NC1_
      tmp0 = tmp0 + C(p,k  ) * u(ij,p)
      tmp1 = tmp1 + C(p,k+1) * u(ij,p)
    end do
    v(ij,k  ) = alpha * tmp0 + beta * v(ij,k  )
    v(ij,k+1) = alpha * tmp1 + beta * v(ij,k+1)
  end do
  end do

#if _NC2_ == _NC2_T2_ + 1

  !$acc loop independent vector
  !DIR$ SIMD
  do ij = 1, na_nb
    tmp0 = 0
    do p = 1, _NC1_
      tmp0 = tmp0 + C(p,_NC2_) * u(ij,p)
    end do
    v(ij,_NC2_) = alpha * tmp0 + beta * v(ij,_NC2_)
  end do

#endif

end subroutine PROC(CtxIxI__,_NC1_,_NC2_)

#undef _NC2_T2_
