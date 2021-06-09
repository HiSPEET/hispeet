!-----------------------------------------------------------------------------
!> Computes  v = α Qᵀ⊗I⊗I u + β v
!>
!>   * k:  block 1,  unroll 4
!>   * j:  block 1,  unroll 1
!>   * i:  block 1,  unroll 1
!>   * p:  block 4,  unroll 1
!>   * i+j merged
!>   * explicit remainder handling

#define _NQ_T4_ (_NQ_ / 4) * 4

subroutine PROC(QtxIxI__,_NQ_)(Q, alpha, beta, u, v)
  !$acc routine vector
  real(RWP), intent(in)    :: Q(_NQ_,_NQ_)      !< square matrix
  real(RWP), intent(in)    :: alpha             !< factor α
  real(RWP), intent(in)    :: beta              !< factor β
  real(RWP), intent(in)    :: u(_NQ_*_NQ_,_NQ_) !< operand
  real(RWP), intent(inout) :: v(_NQ_*_NQ_,_NQ_) !< result

  integer, parameter :: NQ_T4 = _NQ_T4_

  real(RWP) :: tmp0, tmp1, tmp2, tmp3
  integer   :: ij, k, p   ! loop counters
  integer   :: pb         ! block counters

  !$acc loop collapse(3)
  do k = 1, _NQ_
  do ij = 1, _NQ_**2
    v(ij,k) = beta * v(ij,k)
  end do
  end do

#if _NQ_T4_ > 0

  !$acc loop independent vector
  do k = 1, NQ_T4, 4
    !$acc loop seq
    do pb = 1, NQ_T4, 4
      !$acc loop independent vector
      do ij = 1, _NQ_**2
        tmp0 = 0
        tmp1 = 0
        tmp2 = 0
        tmp3 = 0
        do p = pb, pb+3
          tmp0 = tmp0 + Q(p,k  ) * u(ij,p)
          tmp1 = tmp1 + Q(p,k+1) * u(ij,p)
          tmp2 = tmp2 + Q(p,k+2) * u(ij,p)
          tmp3 = tmp3 + Q(p,k+3) * u(ij,p)
        end do
        v(ij,k  ) = alpha * tmp0 + v(ij,k  )
        v(ij,k+1) = alpha * tmp1 + v(ij,k+1)
        v(ij,k+2) = alpha * tmp2 + v(ij,k+2)
        v(ij,k+3) = alpha * tmp3 + v(ij,k+3)
      end do
    end do
  end do

#endif

#if _NQ_ == _NQ_T4_ + 1

  ! remainder: k = 1:nq-1, p = nq ............................................

  p = _NQ_

  !$acc loop collapse(2) vector
  do k = 1, NQ_T4, 4
    do ij = 1, _NQ_**2
      v(ij,k  ) = alpha * Q(p,k  ) * u(ij,p) + v(ij,k  )
      v(ij,k+1) = alpha * Q(p,k+1) * u(ij,p) + v(ij,k+1)
      v(ij,k+2) = alpha * Q(p,k+2) * u(ij,p) + v(ij,k+2)
      v(ij,k+3) = alpha * Q(p,k+3) * u(ij,p) + v(ij,k+3)
    end do
  end do

  ! remainder: k = nq, p = 1:nq ..............................................

  !$acc loop vector
  do ij = 1, _NQ_**2
    tmp0 = 0
    do p = 1, _NQ_
      tmp0 = tmp0 + Q(p,_NQ_) * u(ij,p)
    end do
    v(ij,_NQ_) = alpha * tmp0 + v(ij,_NQ_)
  end do

#elif _NQ_ == _NQ_T4_ + 2

  ! remainder: k = 1:nq-2, p = nq-1:nq .......................................

  !$acc loop collapse(2) vector
  do k = 1, NQ_T4, 4
    do ij = 1, _NQ_**2
      do p = _NQ_-1, _NQ_
        v(ij,k  ) = alpha * Q(p,k  ) * u(ij,p) + v(ij,k  )
        v(ij,k+1) = alpha * Q(p,k+1) * u(ij,p) + v(ij,k+1)
        v(ij,k+2) = alpha * Q(p,k+2) * u(ij,p) + v(ij,k+2)
        v(ij,k+3) = alpha * Q(p,k+3) * u(ij,p) + v(ij,k+3)
      end do
    end do
  end do

  ! remainder: k = nq-1:nq, p = 1:nq .........................................

  !$acc loop vector
  do ij = 1, _NQ_**2
    tmp0 = 0
    tmp1 = 0
    do p = 1, _NQ_
      tmp0 = tmp0 + Q(p,_NQ_-1) * u(ij,p)
      tmp1 = tmp1 + Q(p,_NQ_  ) * u(ij,p)
    end do
    v(ij,_NQ_-1) = alpha * tmp0 + v(ij,_NQ_-1)
    v(ij,_NQ_  ) = alpha * tmp1 + v(ij,_NQ_  )
  end do

#elif _NQ_ == _NQ_T4_ + 3

  ! remainder: k = 1:nq-3, p = nq-2:nq .......................................

  !$acc loop collapse(2) vector
  do k = 1, NQ_T4, 4
    do ij = 1, _NQ_**2
      do p = _NQ_-2, _NQ_
        v(ij,k  ) = alpha * Q(p,k  ) * u(ij,p) + v(ij,k  )
        v(ij,k+1) = alpha * Q(p,k+1) * u(ij,p) + v(ij,k+1)
        v(ij,k+2) = alpha * Q(p,k+2) * u(ij,p) + v(ij,k+2)
        v(ij,k+3) = alpha * Q(p,k+3) * u(ij,p) + v(ij,k+3)
      end do
    end do
  end do

  ! remainder: j = nq-2:nq, p = 1:nq .........................................

  !$acc loop vector
  do ij = 1, _NQ_**2
    tmp0 = 0
    tmp1 = 0
    tmp2 = 0
    do p = 1, _NQ_
      tmp0 = tmp0 + Q(p,_NQ_-2) * u(ij,p)
      tmp1 = tmp1 + Q(p,_NQ_-1) * u(ij,p)
      tmp2 = tmp2 + Q(p,_NQ_  ) * u(ij,p)
    end do
    v(ij,_NQ_-2) = alpha * tmp0 + v(ij,_NQ_-2)
    v(ij,_NQ_-1) = alpha * tmp1 + v(ij,_NQ_-1)
    v(ij,_NQ_  ) = alpha * tmp2 + v(ij,_NQ_  )
  end do

#endif

end subroutine PROC(QtxIxI__,_NQ_)

#undef _NQ_T4
