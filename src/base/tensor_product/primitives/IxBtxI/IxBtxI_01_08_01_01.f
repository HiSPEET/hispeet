-------------------------------------------------------------------------------
!> Computes  v = α I⊗Bᵀ⊗I u + β v
!>
!>   * k:  block 1,  unroll 1
!>   * j:  block 1,  unroll 8
!>   * i:  block 1,  unroll 1
!>   * p:  block 1,  unroll 1
!>   * using Intel SIMD directive

subroutine PROC(IxBtxI__,_NB1_,_NB2_)(na, nc, B, alpha, beta, u, v)
  !$acc routine vector
  integer,   intent(in)    :: na              !< 1st dimension of u,v
  integer,   intent(in)    :: nc              !< 3rd dimension of u,v
  real(RNP), intent(in)    :: B(_NB1_,_NB2_)  !< rectangular matrix
  real(RNP), intent(in)    :: alpha           !< factor α
  real(RNP), intent(in)    :: beta            !< factor β
  real(RNP), intent(in)    :: u(na,_NB1_,nc)  !< operand
  real(RNP), intent(inout) :: v(na,_NB2_,nc)  !< result

  real(RNP) :: tmp0, tmp1, tmp2, tmp3
  real(RNP) :: tmp4, tmp5, tmp6, tmp7
  integer   :: i, j, k, p

#if _NB2_T8_ > 0

  !$acc loop collapse(3) vector
  do k = 1, nc
  do j = 1, _NB2_T8_, 8
  !DIR$ SIMD
  do i = 1, na
    tmp0 = 0
    tmp1 = 0
    tmp2 = 0
    tmp3 = 0
    tmp4 = 0
    tmp5 = 0
    tmp6 = 0
    tmp7 = 0
    do p = 1, _NB1_
      tmp0 = tmp0 + B(p,j  ) * u(i,p,k)
      tmp1 = tmp1 + B(p,j+1) * u(i,p,k)
      tmp2 = tmp2 + B(p,j+2) * u(i,p,k)
      tmp3 = tmp3 + B(p,j+3) * u(i,p,k)
      tmp4 = tmp4 + B(p,j+4) * u(i,p,k)
      tmp5 = tmp5 + B(p,j+5) * u(i,p,k)
      tmp6 = tmp6 + B(p,j+6) * u(i,p,k)
      tmp7 = tmp7 + B(p,j+7) * u(i,p,k)
    end do
    v(i,j  ,k) = alpha * tmp0 + beta * v(i,j  ,k)
    v(i,j+1,k) = alpha * tmp1 + beta * v(i,j+1,k)
    v(i,j+2,k) = alpha * tmp2 + beta * v(i,j+2,k)
    v(i,j+3,k) = alpha * tmp3 + beta * v(i,j+3,k)
    v(i,j+4,k) = alpha * tmp4 + beta * v(i,j+4,k)
    v(i,j+5,k) = alpha * tmp5 + beta * v(i,j+5,k)
    v(i,j+6,k) = alpha * tmp6 + beta * v(i,j+6,k)
    v(i,j+7,k) = alpha * tmp7 + beta * v(i,j+7,k)
  end do
  end do
  end do

#endif

#if na == _NB2_T8_ + 1

  j = _NB2_T8_ + 1

  !$acc loop collapse(2) vector
  do k = 1, nc
  !DIR$ SIMD
  do i = 1, na
    tmp0 = 0
    do p = 1, nc
      tmp0 = tmp0 + B(p,j) * u(i,p,k)
    end do
    v(i,j,k) = alpha * tmp0 + beta * v(i,j,k)
  end do
  end do

#elif na == _NB2_T8_ + 2

  j = _NB2_T8_ + 1

  !$acc loop collapse(2) vector
  do k = 1, nc
  !DIR$ SIMD
  do i = 1, na
    tmp0 = 0
    tmp1 = 0
    do p = 1, _NB1_
      tmp0 = tmp0 + B(p,j  ) * u(i,p,k)
      tmp1 = tmp1 + B(p,j+1) * u(i,p,k)
    end do
    v(i,j  ,k) = alpha * tmp0 + beta * v(i,j  ,k)
    v(i,j+1,k) = alpha * tmp1 + beta * v(i,j+1,k)
  end do
  end do

#elif na == _NB2_T8_ + 3

  j = _NB2_T8_ + 1

  !$acc loop collapse(2) vector
  do k = 1, nc
  !DIR$ SIMD
  do i = 1, na
    tmp0 = 0
    tmp1 = 0
    tmp2 = 0
    do p = 1, _NB1_
      tmp0 = tmp0 + B(p,j  ) * u(i,p,k)
      tmp1 = tmp1 + B(p,j+1) * u(i,p,k)
      tmp2 = tmp2 + B(p,j+2) * u(i,p,k)
    end do
    v(i,j  ,k) = alpha * tmp0 + beta * v(i,j  ,k)
    v(i,j+1,k) = alpha * tmp1 + beta * v(i,j+1,k)
    v(i,j+2,k) = alpha * tmp2 + beta * v(i,j+2,k)
  end do
  end do

#elif na == _NB2_T8_ + 4

  j = _NB2_T8_ + 1

  !$acc loop collapse(2) vector
  do k = 1, nc
  !DIR$ SIMD
  do i = 1, na
    tmp0 = 0
    tmp1 = 0
    tmp2 = 0
    tmp3 = 0
    do p = 1, _NB1_
      tmp0 = tmp0 + B(p,j  ) * u(i,p,k)
      tmp1 = tmp1 + B(p,j+1) * u(i,p,k)
      tmp2 = tmp2 + B(p,j+2) * u(i,p,k)
      tmp3 = tmp3 + B(p,j+3) * u(i,p,k)
    end do
    v(i,j  ,k) = alpha * tmp0 + beta * v(i,j  ,k)
    v(i,j+1,k) = alpha * tmp1 + beta * v(i,j+1,k)
    v(i,j+2,k) = alpha * tmp2 + beta * v(i,j+2,k)
    v(i,j+3,k) = alpha * tmp3 + beta * v(i,j+3,k)
  end do
  end do

#elif na == _NB2_T8_ + 5

  j = _NB2_T8_ + 1

  !$acc loop collapse(2) vector
  do k = 1, nc
  !DIR$ SIMD
  do i = 1, na
    tmp0 = 0
    tmp1 = 0
    tmp2 = 0
    tmp3 = 0
    tmp4 = 0
    do p = 1, _NB1_
      tmp0 = tmp0 + B(p,j  ) * u(i,p,k)
      tmp1 = tmp1 + B(p,j+1) * u(i,p,k)
      tmp2 = tmp2 + B(p,j+2) * u(i,p,k)
      tmp3 = tmp3 + B(p,j+3) * u(i,p,k)
      tmp4 = tmp4 + B(p,j+4) * u(i,p,k)
    end do
    v(i,j  ,k) = alpha * tmp0 + beta * v(i,j  ,k)
    v(i,j+1,k) = alpha * tmp1 + beta * v(i,j+1,k)
    v(i,j+2,k) = alpha * tmp2 + beta * v(i,j+2,k)
    v(i,j+3,k) = alpha * tmp3 + beta * v(i,j+3,k)
    v(i,j+4,k) = alpha * tmp4 + beta * v(i,j+4,k)
  end do
  end do

#elif na == _NB2_T8_ + 6

  j = _NB2_T8_ + 1

  !$acc loop collapse(2) vector
  do k = 1, nc
  !DIR$ SIMD
  do i = 1, na
    tmp0 = 0
    tmp1 = 0
    tmp2 = 0
    tmp3 = 0
    tmp4 = 0
    tmp5 = 0
    do p = 1, _NB1_
      tmp0 = tmp0 + B(p,j  ) * u(i,p,k)
      tmp1 = tmp1 + B(p,j+1) * u(i,p,k)
      tmp2 = tmp2 + B(p,j+2) * u(i,p,k)
      tmp3 = tmp3 + B(p,j+3) * u(i,p,k)
      tmp4 = tmp4 + B(p,j+4) * u(i,p,k)
      tmp5 = tmp5 + B(p,j+5) * u(i,p,k)
    end do
    v(i,j  ,k) = alpha * tmp0 + beta * v(i,j  ,k)
    v(i,j+1,k) = alpha * tmp1 + beta * v(i,j+1,k)
    v(i,j+2,k) = alpha * tmp2 + beta * v(i,j+2,k)
    v(i,j+3,k) = alpha * tmp3 + beta * v(i,j+3,k)
    v(i,j+4,k) = alpha * tmp4 + beta * v(i,j+4,k)
    v(i,j+5,k) = alpha * tmp5 + beta * v(i,j+5,k)
  end do
  end do

#elif na == _NB2_T8_ + 7

  j = _NB2_T8_ + 1

  !$acc loop collapse(2) vector
  do k = 1, nc
  !DIR$ SIMD
  do i = 1, na
    tmp0 = 0
    tmp1 = 0
    tmp2 = 0
    tmp3 = 0
    tmp4 = 0
    tmp5 = 0
    tmp6 = 0
    do p = 1, _NB1_
      tmp0 = tmp0 + B(p,j  ) * u(i,p,k)
      tmp1 = tmp1 + B(p,j+1) * u(i,p,k)
      tmp2 = tmp2 + B(p,j+2) * u(i,p,k)
      tmp3 = tmp3 + B(p,j+3) * u(i,p,k)
      tmp4 = tmp4 + B(p,j+4) * u(i,p,k)
      tmp5 = tmp5 + B(p,j+5) * u(i,p,k)
      tmp6 = tmp6 + B(p,j+6) * u(i,p,k)
    end do
    v(i,j  ,k) = alpha * tmp0 + beta * v(i,j  ,k)
    v(i,j+1,k) = alpha * tmp1 + beta * v(i,j+1,k)
    v(i,j+2,k) = alpha * tmp2 + beta * v(i,j+2,k)
    v(i,j+3,k) = alpha * tmp3 + beta * v(i,j+3,k)
    v(i,j+4,k) = alpha * tmp4 + beta * v(i,j+4,k)
    v(i,j+5,k) = alpha * tmp5 + beta * v(i,j+5,k)
    v(i,j+6,k) = alpha * tmp6 + beta * v(i,j+6,k)
  end do
  end do

#endif

end subroutine PROC(IxBtxI__,_NB1_,_NB2_)
