!-------------------------------------------------------------------------------
!> Computes  v = α I⊗I⊗Qᵀ u + β v
!>
!>   * k:  block 1,  unroll 1
!>   * j:  block 1,  unroll 1
!>   * i:  block 8,  unroll 8
!>   * p:  block 1,  unroll 1
!>   * explicit remainder handling

#define _NQ_T8_ (_NQ_ / 8) * 8

subroutine PROC(IxIxQt__,_NQ_)(Q, alpha, beta, u, v)
  !$acc routine vector
  real(RWP), intent(in)    :: Q(_NQ_,_NQ_)      !< square matrix
  real(RWP), intent(in)    :: alpha             !< factor α
  real(RWP), intent(in)    :: beta              !< factor β
  real(RWP), intent(in)    :: u(_NQ_,_NQ_,_NQ_) !< operand
  real(RWP), intent(inout) :: v(_NQ_,_NQ_,_NQ_) !< result

  integer, parameter :: NQ_T8 = _NQ_T8_

  integer   :: i, j, k, p   ! loop counters
  integer   :: ib           ! block counters

  !$acc loop collapse(3)
  do k = 1, _NQ_
  do j = 1, _NQ_
  do i = 1, _NQ_
    v(i,j,k) = beta * v(i,j,k)
  end do
  end do
  end do

#if _NQ_T8_ > 0

  !$acc loop independent vector
  do k = 1, _NQ_
    !$acc loop seq
    do p = 1, _NQ_
      !$acc loop collapse(2) independent vector
      do j = 1, _NQ_
      do ib = 1, NQ_T8, 8

        v(ib  ,j,k) = alpha * Q(p,ib  ) * u(p,j,k) + v(ib  ,j,k)
        v(ib+1,j,k) = alpha * Q(p,ib+1) * u(p,j,k) + v(ib+1,j,k)
        v(ib+2,j,k) = alpha * Q(p,ib+2) * u(p,j,k) + v(ib+2,j,k)
        v(ib+3,j,k) = alpha * Q(p,ib+3) * u(p,j,k) + v(ib+3,j,k)
        v(ib+4,j,k) = alpha * Q(p,ib+4) * u(p,j,k) + v(ib+4,j,k)
        v(ib+5,j,k) = alpha * Q(p,ib+5) * u(p,j,k) + v(ib+5,j,k)
        v(ib+6,j,k) = alpha * Q(p,ib+6) * u(p,j,k) + v(ib+6,j,k)
        v(ib+7,j,k) = alpha * Q(p,ib+7) * u(p,j,k) + v(ib+7,j,k)

      end do
      end do
    end do
  end do

#endif

#if _NQ_ == _NQ_T8_ + 1

  ib = NQ_T8 + 1

  do k = 1, _NQ_
    do j = 1, _NQ_
      do p = 1, _NQ_

        v(ib  ,j,k) = alpha * Q(p,ib  ) * u(p,j,k) + v(ib  ,j,k)

      end do
    end do
  end do

#elif _NQ_ == _NQ_T8_ + 2

  ib = NQ_T8 + 1

  do k = 1, _NQ_
    do j = 1, _NQ_
      do p = 1, _NQ_

        v(ib  ,j,k) = alpha * Q(p,ib  ) * u(p,j,k) + v(ib  ,j,k)
        v(ib+1,j,k) = alpha * Q(p,ib+1) * u(p,j,k) + v(ib+1,j,k)

      end do
    end do
  end do

#elif _NQ_ == _NQ_T8_ + 3

  ib = NQ_T8 + 1

  do k = 1, _NQ_
    do j = 1, _NQ_
      do p = 1, _NQ_

        v(ib  ,j,k) = alpha * Q(p,ib  ) * u(p,j,k) + v(ib  ,j,k)
        v(ib+1,j,k) = alpha * Q(p,ib+1) * u(p,j,k) + v(ib+1,j,k)
        v(ib+2,j,k) = alpha * Q(p,ib+2) * u(p,j,k) + v(ib+2,j,k)

      end do
    end do
  end do

#elif _NQ_ == _NQ_T8_ + 4

  ib = NQ_T8 + 1

  do k = 1, _NQ_
    do j = 1, _NQ_
      do p = 1, _NQ_

        v(ib  ,j,k) = alpha * Q(p,ib  ) * u(p,j,k) + v(ib  ,j,k)
        v(ib+1,j,k) = alpha * Q(p,ib+1) * u(p,j,k) + v(ib+1,j,k)
        v(ib+2,j,k) = alpha * Q(p,ib+2) * u(p,j,k) + v(ib+2,j,k)
        v(ib+3,j,k) = alpha * Q(p,ib+3) * u(p,j,k) + v(ib+3,j,k)

      end do
    end do
  end do

#elif _NQ_ == _NQ_T8_ + 5

  ib = NQ_T8 + 1

  do k = 1, _NQ_
    do j = 1, _NQ_
      do p = 1, _NQ_

        v(ib  ,j,k) = alpha * Q(p,ib  ) * u(p,j,k) + v(ib  ,j,k)
        v(ib+1,j,k) = alpha * Q(p,ib+1) * u(p,j,k) + v(ib+1,j,k)
        v(ib+2,j,k) = alpha * Q(p,ib+2) * u(p,j,k) + v(ib+2,j,k)
        v(ib+3,j,k) = alpha * Q(p,ib+3) * u(p,j,k) + v(ib+3,j,k)
        v(ib+4,j,k) = alpha * Q(p,ib+4) * u(p,j,k) + v(ib+4,j,k)

      end do
    end do
  end do

#elif _NQ_ == _NQ_T8_ + 6

  ib = NQ_T8 + 1

  do k = 1, _NQ_
    do j = 1, _NQ_
      do p = 1, _NQ_

        v(ib  ,j,k) = alpha * Q(p,ib  ) * u(p,j,k) + v(ib  ,j,k)
        v(ib+1,j,k) = alpha * Q(p,ib+1) * u(p,j,k) + v(ib+1,j,k)
        v(ib+2,j,k) = alpha * Q(p,ib+2) * u(p,j,k) + v(ib+2,j,k)
        v(ib+3,j,k) = alpha * Q(p,ib+3) * u(p,j,k) + v(ib+3,j,k)
        v(ib+4,j,k) = alpha * Q(p,ib+4) * u(p,j,k) + v(ib+4,j,k)
        v(ib+5,j,k) = alpha * Q(p,ib+5) * u(p,j,k) + v(ib+5,j,k)

      end do
    end do
  end do

#elif _NQ_ == _NQ_T8_ + 7

  ib = NQ_T8 + 1

  do k = 1, _NQ_
    do j = 1, _NQ_
      do p = 1, _NQ_

        v(ib  ,j,k) = alpha * Q(p,ib  ) * u(p,j,k) + v(ib  ,j,k)
        v(ib+1,j,k) = alpha * Q(p,ib+1) * u(p,j,k) + v(ib+1,j,k)
        v(ib+2,j,k) = alpha * Q(p,ib+2) * u(p,j,k) + v(ib+2,j,k)
        v(ib+3,j,k) = alpha * Q(p,ib+3) * u(p,j,k) + v(ib+3,j,k)
        v(ib+4,j,k) = alpha * Q(p,ib+4) * u(p,j,k) + v(ib+4,j,k)
        v(ib+5,j,k) = alpha * Q(p,ib+5) * u(p,j,k) + v(ib+5,j,k)
        v(ib+6,j,k) = alpha * Q(p,ib+6) * u(p,j,k) + v(ib+6,j,k)

      end do
    end do
  end do

#endif

end subroutine PROC(IxIxQt__,_NQ_)

#undef _NQ_T8_
