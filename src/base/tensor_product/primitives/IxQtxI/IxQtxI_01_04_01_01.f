!-----------------------------------------------------------------------------
!> Computes  v = α I⊗Qᵀ⊗I u + β v
!>
!>   * k:  block 1,  unroll 1
!>   * j:  block 1,  unroll 4
!>   * i:  block 1,  unroll 1
!>   * p:  block 1,  unroll 1
!>   * explicit remainder handling (r)

#define _NQ_T4_ (_NQ_ / 4) * 4

subroutine PROC(IxQtxI__,_NQ_)(Q, alpha, beta, u, v)
  !$acc routine vector
  real(RWP), intent(in)    :: Q(_NQ_,_NQ_)      !< square matrix
  real(RWP), intent(in)    :: alpha             !< factor α
  real(RWP), intent(in)    :: beta              !< factor β
  real(RWP), intent(in)    :: u(_NQ_,_NQ_,_NQ_) !< operand
  real(RWP), intent(inout) :: v(_NQ_,_NQ_,_NQ_) !< result

  integer, parameter :: NQ_T4 = _NQ_T4_

  real(RWP) :: tmp0, tmp1, tmp2, tmp3
  integer   :: i, j, k, p

  !$acc loop collapse(3)
  do k = 1, _NQ_
  do j = 1, _NQ_
  do i = 1, _NQ_
    v(i,j,k) = beta * v(i,j,k)
  end do
  end do
  end do

#if _NQ_T4_ > 0


  !$acc loop collapse(3) vector
  do k = 1, _NQ_
  do j = 1, NQ_T4, 4
  do i = 1, _NQ_
    tmp0 = 0
    tmp1 = 0
    tmp2 = 0
    tmp3 = 0
    do p = 1, _NQ_
      tmp0 = tmp0 + Q(p,j  ) * u(i,p,k)
      tmp1 = tmp1 + Q(p,j+1) * u(i,p,k)
      tmp2 = tmp2 + Q(p,j+2) * u(i,p,k)
      tmp3 = tmp3 + Q(p,j+3) * u(i,p,k)
    end do
    v(i,j  ,k) = alpha * tmp0 + v(i,j  ,k)
    v(i,j+1,k) = alpha * tmp1 + v(i,j+1,k)
    v(i,j+2,k) = alpha * tmp2 + v(i,j+2,k)
    v(i,j+3,k) = alpha * tmp3 + v(i,j+3,k)
  end do
  end do
  end do

#endif

#if _NQ_ == _NQ_T4_ + 1

  ! remainder: j = nq ........................................................

  !$acc loop collapse(2) vector
  do k = 1, _NQ_
  do i = 1, _NQ_
    tmp0 = 0
    do p = 1, _NQ_
      tmp0 = tmp0 + Q(p,_NQ_) * u(i,p,k)
    end do
    v(i,_NQ_,k) = alpha * tmp0 + v(i,_NQ_,k)
  end do
  end do

#elif _NQ_ == _NQ_T4_ + 2

  ! remainder: j = nq-1:nq ...................................................

  !$acc loop collapse(2) vector
  do k = 1, _NQ_
  do i = 1, _NQ_
    tmp0 = 0
    tmp1 = 0
    do p = 1, _NQ_
      tmp0 = tmp0 + Q(p,_NQ_-1) * u(i,p,k)
      tmp1 = tmp1 + Q(p,_NQ_  ) * u(i,p,k)
    end do
    v(i,_NQ_-1,k) = alpha * tmp0 + v(i,_NQ_-1,k)
    v(i,_NQ_  ,k) = alpha * tmp1 + v(i,_NQ_  ,k)
  end do
  end do

#elif _NQ_ == _NQ_T4_ + 3

  ! remainder: j = nq-2:nq ...................................................

  !$acc loop collapse(2) vector
  do k = 1, _NQ_
  do i = 1, _NQ_
    tmp0 = 0
    tmp1 = 0
    tmp2 = 0
    do p = 1, _NQ_
      tmp0 = tmp0 + Q(p,_NQ_-2) * u(i,p,k)
      tmp1 = tmp1 + Q(p,_NQ_-1) * u(i,p,k)
      tmp2 = tmp2 + Q(p,_NQ_  ) * u(i,p,k)
    end do
    v(i,_NQ_-2,k) = alpha * tmp0 + v(i,_NQ_-2,k)
    v(i,_NQ_-1,k) = alpha * tmp1 + v(i,_NQ_-1,k)
    v(i,_NQ_  ,k) = alpha * tmp2 + v(i,_NQ_  ,k)
  end do
  end do

#endif

end subroutine PROC(IxQtxI__,_NQ_)

#undef _NQ_T4
