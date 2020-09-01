!-------------------------------------------------------------------------------
!> Computes  v = α I⊗Bᵀ⊗I u + β v
!>
!>   * k:  block 1,  unroll 1
!>   * j:  block 1,  unroll 2
!>   * i:  block 1,  unroll 1
!>   * p:  block 1,  unroll 1
!>   * using Intel SIMD directive

# define _NB1_T2_ (_NB1_ / 2) * 2

subroutine PROC(IxBtxI__,_NB1_,_NB2_)(na, nc, Bt, alpha, beta, u, v)
  !$acc routine vector
  integer,   intent(in)    :: na              !< 1st dimension of u,v
  integer,   intent(in)    :: nc              !< 3rd dimension of u,v
  real(RNP), intent(in)    :: Bt(_NB2_,_NB1_)  !< rectangular matrix
  real(RNP), intent(in)    :: alpha           !< factor α
  real(RNP), intent(in)    :: beta            !< factor β
  real(RNP), intent(in)    :: u(na,_NB2_,nc)  !< operand
  real(RNP), intent(inout) :: v(na,_NB1_,nc)  !< result

  real(RNP) :: tmp0, tmp1
  integer   :: i, j, k, p

#if _NB1_T2_ > 0

  !$acc loop collapse(3) vector
  do k = 1, nc
  do j = 1, _NB1_T2_, 2
  do i = 1, na
    tmp0 = 0
    tmp1 = 0
    do p = 1, _NB2_
      tmp0 = tmp0 + Bt(p,j  ) * u(i,p,k)
      tmp1 = tmp1 + Bt(p,j+1) * u(i,p,k)
    end do
    v(i,j  ,k) = alpha * tmp0 + beta * v(i,j  ,k)
    v(i,j+1,k) = alpha * tmp1 + beta * v(i,j+1,k)
  end do
  end do
  end do

#endif

#if _NB1_ == _NB1_T2_ + 1

  !$acc loop collapse(2) vector
  do k = 1, nc
  !DIR$ SIMD
  do i = 1, na
    tmp0 = 0
    do p = 1, _NB2_
      tmp0 = tmp0 + Bt(p,_NB1_) * u(i,p,k)
    end do
    v(i,_NB1_,k) = alpha * tmp0 + beta * v(i,_NB1_,k)
  end do
  end do

#endif

end subroutine PROC(IxBtxI__,_NB1_,_NB2_)

#undef _NB1_T2_
