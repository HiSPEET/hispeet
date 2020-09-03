!-------------------------------------------------------------------------------
!> Computes  v = α Cᵀ⊗I⊗I u + β v
!>
!>   * k:  block 1,  unroll 8
!>   * j:  joined with i
!>   * i:  joined with j
!>   * p:  block 1,  unroll 1
!>   * using Intel SIMD directive

#define _NC2_T8_ (_NC2_ / 8) * 8

subroutine PROC(CtxIxI__,_NC1_,_NC2_)(na_nb, C, alpha, beta, u, v)
  !$acc routine vector
  integer,   intent(in)    :: na_nb           !< size of dimensions 1+2 of u,v
  real(RNP), intent(in)    :: C(_NC1_,_NC2_)  !< rectangular matrix
  real(RNP), intent(in)    :: alpha           !< factor α
  real(RNP), intent(in)    :: beta            !< factor β
  real(RNP), intent(in)    :: u(na_nb,_NC1_)  !< operand
  real(RNP), intent(inout) :: v(na_nb,_NC2_)  !< result

  real(RNP) :: tmp0, tmp1, tmp2, tmp3
  real(RNP) :: tmp4, tmp5, tmp6, tmp7
  integer   :: ij, k, p

#if _NC2_T8_ > 0

  !$acc loop collapse(2) vector
  do k = 1, _NC2_T8_, 8
  !DIR$ SIMD
  do ij = 1, na_nb
    tmp0 = 0
    tmp1 = 0
    tmp2 = 0
    tmp3 = 0
    tmp4 = 0
    tmp5 = 0
    tmp6 = 0
    tmp7 = 0
    do p = 1, _NC1_
      tmp0 = tmp0 + C(p,k  ) * u(ij,p)
      tmp1 = tmp1 + C(p,k+1) * u(ij,p)
      tmp2 = tmp2 + C(p,k+2) * u(ij,p)
      tmp3 = tmp3 + C(p,k+3) * u(ij,p)
      tmp4 = tmp4 + C(p,k+4) * u(ij,p)
      tmp5 = tmp5 + C(p,k+5) * u(ij,p)
      tmp6 = tmp6 + C(p,k+6) * u(ij,p)
      tmp7 = tmp7 + C(p,k+7) * u(ij,p)
    end do
    v(ij,k  ) = alpha * tmp0 + beta * v(ij,k  )
    v(ij,k+1) = alpha * tmp1 + beta * v(ij,k+1)
    v(ij,k+2) = alpha * tmp2 + beta * v(ij,k+2)
    v(ij,k+3) = alpha * tmp3 + beta * v(ij,k+3)
    v(ij,k+4) = alpha * tmp4 + beta * v(ij,k+4)
    v(ij,k+5) = alpha * tmp5 + beta * v(ij,k+5)
    v(ij,k+6) = alpha * tmp6 + beta * v(ij,k+6)
    v(ij,k+7) = alpha * tmp7 + beta * v(ij,k+7)
  end do
  end do

#endif

#if _NC2_ == _NC2_T8_ + 1

  k =  _NC2_T8_ + 1

  !$acc loop vector
  !DIR$ SIMD
  do ij = 1, na_nb
    tmp0 = 0
    do p = 1, _NC1_
      tmp0 = tmp0 + C(p,k) * u(ij,p)
    end do
    v(ij,k) = alpha * tmp0 + beta * v(ij,k)
  end do

#elif _NC2_ == _NC2_T8_ + 2

  k =  _NC2_T8_ + 1

  !$acc loop vector
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

#elif _NC2_ == _NC2_T8_ + 3

  k =  _NC2_T8_ + 1

  !$acc loop vector
  !DIR$ SIMD
  do ij = 1, na_nb
    tmp0 = 0
    tmp1 = 0
    tmp2 = 0
    do p = 1, _NC1_
      tmp0 = tmp0 + C(p,k  ) * u(ij,p)
      tmp1 = tmp1 + C(p,k+1) * u(ij,p)
      tmp2 = tmp2 + C(p,k+2) * u(ij,p)
    end do
    v(ij,k  ) = alpha * tmp0 + beta * v(ij,k  )
    v(ij,k+1) = alpha * tmp1 + beta * v(ij,k+1)
    v(ij,k+2) = alpha * tmp2 + beta * v(ij,k+2)
  end do

#elif _NC2_ == _NC2_T8_ + 4

  k =  _NC2_T8_ + 1

  !$acc loop vector
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

#elif _NC2_ == _NC2_T8_ + 5

  k =  _NC2_T8_ + 1

  !$acc loop vector
  !DIR$ SIMD
  do ij = 1, na_nb
    tmp0 = 0
    tmp1 = 0
    tmp2 = 0
    tmp3 = 0
    tmp4 = 0
    do p = 1, _NC1_
      tmp0 = tmp0 + C(p,k  ) * u(ij,p)
      tmp1 = tmp1 + C(p,k+1) * u(ij,p)
      tmp2 = tmp2 + C(p,k+2) * u(ij,p)
      tmp3 = tmp3 + C(p,k+3) * u(ij,p)
      tmp4 = tmp4 + C(p,k+4) * u(ij,p)
    end do
    v(ij,k  ) = alpha * tmp0 + beta * v(ij,k  )
    v(ij,k+1) = alpha * tmp1 + beta * v(ij,k+1)
    v(ij,k+2) = alpha * tmp2 + beta * v(ij,k+2)
    v(ij,k+3) = alpha * tmp3 + beta * v(ij,k+3)
    v(ij,k+4) = alpha * tmp4 + beta * v(ij,k+4)
  end do

#elif __NA1__ == __NA1_T8__ + 6

  k =  _NC2_T8_ + 1

  !$acc loop vector
  !DIR$ SIMD
  do ij = 1, na_nb
    tmp0 = 0
    tmp1 = 0
    tmp2 = 0
    tmp3 = 0
    tmp4 = 0
    tmp5 = 0
    do p = 1, _NC1_
      tmp0 = tmp0 + C(p,k  ) * u(ij,p)
      tmp1 = tmp1 + C(p,k+1) * u(ij,p)
      tmp2 = tmp2 + C(p,k+2) * u(ij,p)
      tmp3 = tmp3 + C(p,k+3) * u(ij,p)
      tmp4 = tmp4 + C(p,k+4) * u(ij,p)
      tmp5 = tmp5 + C(p,k+5) * u(ij,p)
    end do
    v(ij,k  ) = alpha * tmp0 + beta * v(ij,k  )
    v(ij,k+1) = alpha * tmp1 + beta * v(ij,k+1)
    v(ij,k+2) = alpha * tmp2 + beta * v(ij,k+2)
    v(ij,k+3) = alpha * tmp3 + beta * v(ij,k+3)
    v(ij,k+4) = alpha * tmp4 + beta * v(ij,k+4)
    v(ij,k+5) = alpha * tmp5 + beta * v(ij,k+5)
  end do

#elif _NC2_ == _NC2_T8_ + 7

  k =  _NC2_T8_ + 1

  !$acc loop vector
  !DIR$ SIMD
  do ij = 1, na_nb
    tmp0 = 0
    tmp1 = 0
    tmp2 = 0
    tmp3 = 0
    tmp4 = 0
    tmp5 = 0
    tmp6 = 0
    do p = 1, _NC1_
      tmp0 = tmp0 + C(p,k  ) * u(ij,p)
      tmp1 = tmp1 + C(p,k+1) * u(ij,p)
      tmp2 = tmp2 + C(p,k+2) * u(ij,p)
      tmp3 = tmp3 + C(p,k+3) * u(ij,p)
      tmp4 = tmp4 + C(p,k+4) * u(ij,p)
      tmp5 = tmp5 + C(p,k+5) * u(ij,p)
      tmp6 = tmp6 + C(p,k+6) * u(ij,p)
    end do
    v(ij,k  ) = alpha * tmp0 + beta * v(ij,k  )
    v(ij,k+1) = alpha * tmp1 + beta * v(ij,k+1)
    v(ij,k+2) = alpha * tmp2 + beta * v(ij,k+2)
    v(ij,k+3) = alpha * tmp3 + beta * v(ij,k+3)
    v(ij,k+4) = alpha * tmp4 + beta * v(ij,k+4)
    v(ij,k+5) = alpha * tmp5 + beta * v(ij,k+5)
    v(ij,k+6) = alpha * tmp6 + beta * v(ij,k+6)
  end do

#endif

end subroutine PROC(CtxIxI__,_NC1_,_NC2_)

#undef _NC2_T8_
