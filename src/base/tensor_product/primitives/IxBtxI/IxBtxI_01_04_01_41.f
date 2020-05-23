!-------------------------------------------------------------------------------
!> Computes  v = α I⊗Bᵀ⊗I u + β v
!>
!>   * k:  block 1,  unroll 1
!>   * j:  block 1,  unroll 4
!>   * i:  block 1,  unroll 1
!>   * p:  block 4,  unroll 1
!>   * using Intel SIMD directive

subroutine PROC(IxBtxI_,_NB1_,_NB2_)(na, nc, B, alpha, beta, u, v)
  !$acc routine vector
  integer,   intent(in)    :: na              !< 1st dimension of u,v
  integer,   intent(in)    :: nc              !< 3rd dimension of u,v
  real(RNP), intent(in)    :: B(_NB1_,_NB2_)  !< rectangular matrix
  real(RNP), intent(in)    :: alpha           !< factor α
  real(RNP), intent(in)    :: beta            !< factor β
  real(RNP), intent(in)    :: u(na,_NB1_,nc)  !< operand
  real(RNP), intent(inout) :: v(na,_NB2_,nc)  !< result

  real(RNP) :: tmp0, tmp1
  integer   :: i, j, k, p, pb


#if _NB1_ > 0

  !$acc loop collapse(2) independent vector
  do k = 1, nc
  do j = 1, _NB2_T4_, 4
    !$acc loop seq
    do pb = 1, _NB1_T4_, 4
      !$acc loop independent vector
      do i = 1, na
        do p = pb, pb + 3
          v(i,j  ,k) = alpha * B(p,j  ) * u(i,p,k) + beta * v(i,j  ,k)
          v(i,j+1,k) = alpha * B(p,j+1) * u(i,p,k) + beta * v(i,j+1,k)
          v(i,j+2,k) = alpha * B(p,j+2) * u(i,p,k) + beta * v(i,j+2,k)
          v(i,j+3,k) = alpha * B(p,j+3) * u(i,p,k) + beta * v(i,j+3,k)
        end do
      end do
    end do
  end do
  end do

#if _NB1_ == _NB1_T4_ + 1

  !$acc loop collapse(3) vector
  do k = 1, nc
  do j = 1, _NB2_T4_, 4
  do i = 1, na
    v(i,j  ,k) = alpha * B(_NB1_,j  ) * u(i,nc,k) + beta * v(i,j  ,k)
    v(i,j+1,k) = alpha * B(_NB1_,j+1) * u(i,nc,k) + beta * v(i,j+1,k)
    v(i,j+2,k) = alpha * B(_NB1_,j+2) * u(i,nc,k) + beta * v(i,j+2,k)
    v(i,j+3,k) = alpha * B(_NB1_,j+3) * u(i,nc,k) + beta * v(i,j+3,k)
  end do
  end do
  end do

#elif _NB1_ == _NB1_T4_ + 2

  !$acc loop collapse(3) vector
  do k = 1, nc
  do j = 1, _NB2_T4_, 4
  do i = 1, na
    do p = _NB1_-1, _NB1_
      v(i,j  ,k) = alpha * B(p,j  ) * u(i,p,k) + beta * v(i,j  ,k)
      v(i,j+1,k) = alpha * B(p,j+1) * u(i,p,k) + beta * v(i,j+1,k)
      v(i,j+2,k) = alpha * B(p,j+2) * u(i,p,k) + beta * v(i,j+2,k)
      v(i,j+3,k) = alpha * B(p,j+3) * u(i,p,k) + beta * v(i,j+3,k)
    end do
  end do
  end do
  end do

#elif _NB1_ == _NB1_T4_ + 3

  !$acc loop collapse(3) vector
  do k = 1, nc
  do j = 1, _NB2_T4_, 4
  do i = 1, na
    do p = _NB1_-2, _NB1_
      v(i,j  ,k) = alpha * B(p,j  ) * u(i,p,k) + beta * v(i,j  ,k)
      v(i,j+1,k) = alpha * B(p,j+1) * u(i,p,k) + beta * v(i,j+1,k)
      v(i,j+2,k) = alpha * B(p,j+2) * u(i,p,k) + beta * v(i,j+2,k)
      v(i,j+3,k) = alpha * B(p,j+3) * u(i,p,k) + beta * v(i,j+3,k)
    end do
  end do
  end do
  end do

#endif

#if _NB2_ == _NB2_T4_ + 1

  !$acc loop collapse(2) vector
  do k = 1, nc
  do i = 1, na
    do p = 1, _NB1_
      v(i,_NB2_,k) = alpha * B(p,_NB2_) * u(i,p,k) + beta * v(i,_NB2_,k)
    end do
  end do
  end do

#elif _NB2_ == _NB2_T4_ + 2

  !$acc loop collapse(2) vector
  do k = 1, nc
  do i = 1, na
    do p = 1, _NB1_
      v(i,_NB2_-1,k) = alpha * B(p,_NB2_-1) * u(i,p,k) + beta * v(i,_NB2_-1,k)
      v(i,_NB2_  ,k) = alpha * B(p,_NB2_  ) * u(i,p,k) + beta * v(i,_NB2_  ,k)
    end do
  end do
  end do

#elif _NB2_ == _NB2_T4_ + 3

  !$acc loop collapse(2) vector
  do k = 1, nc
  do i = 1, na
    do p = 1, _NB1_
      v(i,_NB2_-2,k) = alpha * B(p,_NB2_-2) * u(i,p,k) + beta * v(i,_NB2_-2,k)
      v(i,_NB2_-1,k) = alpha * B(p,_NB2_-1) * u(i,p,k) + beta * v(i,_NB2_-1,k)
      v(i,_NB2_  ,k) = alpha * B(p,_NB2_  ) * u(i,p,k) + beta * v(i,_NB2_  ,k)
    end do
  end do
  end do

#endif

end subroutine PROC(IxBtxI_,_NB1_,_NB2_)
