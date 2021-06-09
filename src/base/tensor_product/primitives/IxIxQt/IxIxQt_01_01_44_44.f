!-------------------------------------------------------------------------------
!> Computes  v = α I⊗I⊗Qᵀ u + β v
!>
!>   * k:  block 1,  unroll 1
!>   * j:  block 1,  unroll 1
!>   * i:  block 4,  unroll 4
!>   * p:  block 4,  unroll 4
!>   * explicit remainder handling

#define _NQ_T4_ (_NQ_ / 4) * 4

subroutine PROC(IxIxQt__,_NQ_)(Q, alpha, beta, u, v)
  !$acc routine vector
  real(RWP), intent(in)    :: Q(_NQ_,_NQ_)      !< square matrix
  real(RWP), intent(in)    :: alpha             !< factor α
  real(RWP), intent(in)    :: beta              !< factor β
  real(RWP), intent(in)    :: u(_NQ_,_NQ_,_NQ_) !< operand
  real(RWP), intent(inout) :: v(_NQ_,_NQ_,_NQ_) !< result

  integer, parameter :: NQ_T4 = _NQ_T4_

  real(RWP) :: tmp0, tmp1, tmp2, tmp3
  integer   :: i, j, k, p   ! loop counters
  integer   :: ib, pb       ! block counters

  !$acc loop collapse(3)
  do k = 1, _NQ_
  do j = 1, _NQ_
  do i = 1, _NQ_
    v(i,j,k) = beta * v(i,j,k)
  end do
  end do
  end do

#if _NQ_T4_ > 0

  !$acc loop independent vector
  do k = 1, _NQ_
    !$acc loop seq
    do pb = 1, NQ_T4, 4
      !$acc loop collapse(2) independent vector
      do j = 1, _NQ_
      do ib = 1, NQ_T4, 4

        tmp0 =        Q(pb  ,ib  ) * u(pb  ,j,k)
        tmp1 =        Q(pb  ,ib+1) * u(pb  ,j,k)
        tmp2 =        Q(pb  ,ib+2) * u(pb  ,j,k)
        tmp3 =        Q(pb  ,ib+3) * u(pb  ,j,k)

        tmp0 = tmp0 + Q(pb+1,ib  ) * u(pb+1,j,k)
        tmp1 = tmp1 + Q(pb+1,ib+1) * u(pb+1,j,k)
        tmp2 = tmp2 + Q(pb+1,ib+2) * u(pb+1,j,k)
        tmp3 = tmp3 + Q(pb+1,ib+3) * u(pb+1,j,k)

        tmp0 = tmp0 + Q(pb+2,ib  ) * u(pb+2,j,k)
        tmp1 = tmp1 + Q(pb+2,ib+1) * u(pb+2,j,k)
        tmp2 = tmp2 + Q(pb+2,ib+2) * u(pb+2,j,k)
        tmp3 = tmp3 + Q(pb+2,ib+3) * u(pb+2,j,k)

        tmp0 = tmp0 + Q(pb+3,ib  ) * u(pb+3,j,k)
        tmp1 = tmp1 + Q(pb+3,ib+1) * u(pb+3,j,k)
        tmp2 = tmp2 + Q(pb+3,ib+2) * u(pb+3,j,k)
        tmp3 = tmp3 + Q(pb+3,ib+3) * u(pb+3,j,k)

        v(ib  ,j,k) = alpha * tmp0 + v(ib  ,j,k)
        v(ib+1,j,k) = alpha * tmp1 + v(ib+1,j,k)
        v(ib+2,j,k) = alpha * tmp2 + v(ib+2,j,k)
        v(ib+3,j,k) = alpha * tmp3 + v(ib+3,j,k)

      end do
      end do
    end do
  end do

#endif

#if _NQ_ == _NQ_T4_ + 1

  ! 1:nq/4*4, nq -----------------------------------------------------------

  pb = NQ_T4 + 1

  !$acc loop collapse(3) independent vector
  do k = 1, _NQ_
    do j = 1, _NQ_
    do ib = 1, NQ_T4, 4

      v(ib  ,j,k) = alpha * Q(pb,ib  ) * u(pb,j,k) + v(ib  ,j,k)
      v(ib+1,j,k) = alpha * Q(pb,ib+1) * u(pb,j,k) + v(ib+1,j,k)
      v(ib+2,j,k) = alpha * Q(pb,ib+2) * u(pb,j,k) + v(ib+2,j,k)
      v(ib+3,j,k) = alpha * Q(pb,ib+3) * u(pb,j,k) + v(ib+3,j,k)

    end do
    end do
  end do

  ! nq, 1:nq ---------------------------------------------------------------

  ib = NQ_T4 + 1

  !$acc loop collapse(2) independent vector
  do k = 1, _NQ_
    do j = 1, _NQ_
      do p = 1, _NQ_

        v(ib,j,k) = alpha * Q(p,ib) * u(p,j,k) + v(ib,j,k)

      end do
    end do
  end do

#elif _NQ_ == _NQ_T4_ + 2

  ! 1:nq/4*4, nq-2:nq ----------------------------------------------------

  pb = NQ_T4 + 1

  !$acc loop collapse(3) independent vector
  do k = 1, _NQ_
    do j = 1, _NQ_
    do ib = 1, NQ_T4, 4

      tmp0 =        Q(pb  ,ib  ) * u(pb  ,j,k)
      tmp1 =        Q(pb  ,ib+1) * u(pb  ,j,k)
      tmp2 =        Q(pb  ,ib+2) * u(pb  ,j,k)
      tmp3 =        Q(pb  ,ib+3) * u(pb  ,j,k)

      tmp0 = tmp0 + Q(pb+1,ib  ) * u(pb+1,j,k)
      tmp1 = tmp1 + Q(pb+1,ib+1) * u(pb+1,j,k)
      tmp2 = tmp2 + Q(pb+1,ib+2) * u(pb+1,j,k)
      tmp3 = tmp3 + Q(pb+1,ib+3) * u(pb+1,j,k)

      v(ib  ,j,k) = alpha * tmp0 + v(ib  ,j,k)
      v(ib+1,j,k) = alpha * tmp1 + v(ib+1,j,k)
      v(ib+2,j,k) = alpha * tmp2 + v(ib+2,j,k)
      v(ib+3,j,k) = alpha * tmp3 + v(ib+3,j,k)

    end do
    end do
  end do

  ! nq-2:nq, 1:nq --------------------------------------------------------

  ib = NQ_T4 + 1

  !$acc loop collapse(2) independent vector
  do k = 1, _NQ_
    do j = 1, _NQ_
    do i = ib, ib+1
      do p = 1, _NQ_

        v(i,j,k) = alpha * Q(p,i) * u(p,j,k) + v(i,j,k)

      end do
    end do
    end do
  end do

#elif _NQ_ == _NQ_T4_ + 3

  ! 1:nq/4*4, nq-2:nq ----------------------------------------------------

  pb = NQ_T4 + 1

  !$acc loop collapse(3) independent vector
  do k = 1, _NQ_
    do j = 1, _NQ_
    do ib = 1, NQ_T4, 4

      tmp0 =        Q(pb  ,ib  ) * u(pb  ,j,k)
      tmp1 =        Q(pb  ,ib+1) * u(pb  ,j,k)
      tmp2 =        Q(pb  ,ib+2) * u(pb  ,j,k)
      tmp3 =        Q(pb  ,ib+3) * u(pb  ,j,k)

      tmp0 = tmp0 + Q(pb+1,ib  ) * u(pb+1,j,k)
      tmp1 = tmp1 + Q(pb+1,ib+1) * u(pb+1,j,k)
      tmp2 = tmp2 + Q(pb+1,ib+2) * u(pb+1,j,k)
      tmp3 = tmp3 + Q(pb+1,ib+3) * u(pb+1,j,k)

      tmp0 = tmp0 + Q(pb+2,ib  ) * u(pb+2,j,k)
      tmp1 = tmp1 + Q(pb+2,ib+1) * u(pb+2,j,k)
      tmp2 = tmp2 + Q(pb+2,ib+2) * u(pb+2,j,k)
      tmp3 = tmp3 + Q(pb+2,ib+3) * u(pb+2,j,k)

      v(ib  ,j,k) = alpha * tmp0 + v(ib  ,j,k)
      v(ib+1,j,k) = alpha * tmp1 + v(ib+1,j,k)
      v(ib+2,j,k) = alpha * tmp2 + v(ib+2,j,k)
      v(ib+3,j,k) = alpha * tmp3 + v(ib+3,j,k)

    end do
    end do
  end do

  ! nq-2:nq, 1:nq --------------------------------------------------------

  ib = NQ_T4 + 1

  !$acc loop collapse(2) independent vector
  do k = 1, _NQ_
    do j = 1, _NQ_
    do i = ib, ib+2
      do p = 1, _NQ_

        v(i,j,k) = alpha * Q(p,i) * u(p,j,k) + v(i,j,k)

      end do
    end do
    end do
  end do

#endif

end subroutine PROC(IxIxQt__,_NQ_)

#undef _NQ_T4_
