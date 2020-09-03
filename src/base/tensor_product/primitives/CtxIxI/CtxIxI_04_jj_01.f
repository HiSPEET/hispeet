!-------------------------------------------------------------------------------
!> Computes  v = α Cᵀ⊗I⊗I u + β v
!>
!>   * k:  block 1,  unroll 4
!>   * j:  joined with i
!>   * i:  joined with j
!>   * p:  block 1,  unroll 1
!>   * using Intel SIMD directive

#define _NC2_T4_ (_NC2_ / 4) * 4

subroutine PROC(CtxIxI__,_NC1_,_NC2_)(na_nb, C, alpha, beta, u, v)
  !$acc routine vector
  integer,   intent(in)    :: na_nb           !< size of dimensions 1+2 of u,v
  real(RNP), intent(in)    :: C(_NC1_,_NC2_)  !< rectangular matrix
  real(RNP), intent(in)    :: alpha           !< factor α
  real(RNP), intent(in)    :: beta            !< factor β
  real(RNP), intent(in)    :: u(na_nb,_NC1_)  !< operand
  real(RNP), intent(inout) :: v(na_nb,_NC2_)  !< result

  real(RNP) :: tmp0, tmp1, tmp2, tmp3
  integer   :: ij, k, p

#if _NC2_T4_ > 0

  !$acc loop collapse(2) vector
  do k = 1, _NC2_T4_, 4
  !DIR$ SIMD
  do ij = 1, na_nb
    tmp0 = 0
    tmp1 = 0
    tmp2 = 0
    tmp3 = 0
    do p = 1, _NC1_
      tmp0 = tmp0 + C(p,k  ) * u(ij,p)
      tmp1 = tmp1 + C(p,k+1) * u(ij,p)
      tmp2 = tmp2 + C(p,k+2) * u(ij,p)
      tmp3 = tmp3 + C(p,k+3) * u(ij,p)
    end do
    v(ij,k  ) = alpha * tmp0 + beta * v(ij,k  )
    v(ij,k+1) = alpha * tmp1 + beta * v(ij,k+1)
    v(ij,k+2) = alpha * tmp2 + beta * v(ij,k+2)
    v(ij,k+3) = alpha * tmp3 + beta * v(ij,k+3)
  end do
  end do

#endif

#if _NC2_ == _NC2_T4_ + 1

  !$acc loop vector
  !DIR$ SIMD
  do ij = 1, na_nb
    tmp0 = 0
    do p = 1, _NC1_
      tmp0 = tmp0 + C(p,_NC2_) * u(ij,p)
    end do
    v(ij,_NC2_) = alpha * tmp0 + beta * v(ij,_NC2_)
  end do

#elif _NC2_ == _NC2_T4_ + 2

  !$acc loop vector
  !DIR$ SIMD
  do ij = 1, na_nb
    tmp0 = 0
    tmp1 = 0
    do p = 1, _NC1_
      tmp0 = tmp0 + C(p,_NC2_-1) * u(ij,p)
      tmp1 = tmp1 + C(p,_NC2_  ) * u(ij,p)
    end do
    v(ij,_NC2_-1) = alpha * tmp0 + beta * v(ij,_NC2_-1)
    v(ij,_NC2_  ) = alpha * tmp1 + beta * v(ij,_NC2_  )
  end do

#elif _NC2_ == _NC2_T4_ + 3

  !$acc loop vector
  !DIR$ SIMD
  do ij = 1, na_nb
    tmp0 = 0
    tmp1 = 0
    tmp2 = 0
    do p = 1, _NC1_
      tmp0 = tmp0 + C(p,_NC2_-2) * u(ij,p)
      tmp1 = tmp1 + C(p,_NC2_-1) * u(ij,p)
      tmp2 = tmp2 + C(p,_NC2_  ) * u(ij,p)
    end do
    v(ij,_NC2_-2) = alpha * tmp0 + beta * v(ij,_NC2_-2)
    v(ij,_NC2_-1) = alpha * tmp1 + beta * v(ij,_NC2_-1)
    v(ij,_NC2_  ) = alpha * tmp2 + beta * v(ij,_NC2_  )
  end do

#endif

end subroutine PROC(CtxIxI__,_NC1_,_NC2_)

#undef _NC2_T4_
