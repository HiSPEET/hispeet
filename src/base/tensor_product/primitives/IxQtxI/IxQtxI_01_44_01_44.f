!-----------------------------------------------------------------------------
!> Computes  v = α I⊗Qᵀ⊗I u + β v
!>
!>   * k:  block 1,  unroll 1
!>   * j:  block 4,  unroll 4
!>   * i:  block 1,  unroll 1
!>   * p:  block 4,  unroll 4
!>   * explicit remainder handling (r)

!> loop orders:
!> k_i_pb_jb as well as
!> k_pb_i_jb where also investigated
!> and proven to be less performant
!> than k_pb_jb_i

#define _NQ_T4_ (_NQ_ / 4) * 4

subroutine PROC(IxQtxI__,_NQ_)(Q, alpha, beta, u, v)
  !$acc routine vector
  real(RNP), intent(in)    :: Q(_NQ_,_NQ_)      !< square matrix
  real(RNP), intent(in)    :: alpha             !< factor α
  real(RNP), intent(in)    :: beta              !< factor β
  real(RNP), intent(in)    :: u(_NQ_,_NQ_,_NQ_) !< operand
  real(RNP), intent(inout) :: v(_NQ_,_NQ_,_NQ_) !< result

  integer, parameter :: NQ_T4 = _NQ_T4_

  real(RNP) :: tmp0, tmp1, tmp2, tmp3
  integer   :: i, j, k, p   ! loop counters
  integer   :: jb, pb       ! block counters
  
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
  do k  = 1, _NQ_
    !$acc loop seq
    do pb = 1, NQ_T4, 4
      !$acc loop collapse(2) independent vector
      do jb = 1, NQ_T4, 4
      do i  = 1, _NQ_

        tmp0 =        Q(pb  ,jb  ) * u(i,pb  ,k)
        tmp1 =        Q(pb  ,jb+1) * u(i,pb  ,k)
        tmp2 =        Q(pb  ,jb+2) * u(i,pb  ,k)
        tmp3 =        Q(pb  ,jb+3) * u(i,pb  ,k)

        tmp0 = tmp0 + Q(pb+1,jb  ) * u(i,pb+1,k)
        tmp1 = tmp1 + Q(pb+1,jb+1) * u(i,pb+1,k)
        tmp2 = tmp2 + Q(pb+1,jb+2) * u(i,pb+1,k)
        tmp3 = tmp3 + Q(pb+1,jb+3) * u(i,pb+1,k)

        tmp0 = tmp0 + Q(pb+2,jb  ) * u(i,pb+2,k)
        tmp1 = tmp1 + Q(pb+2,jb+1) * u(i,pb+2,k)
        tmp2 = tmp2 + Q(pb+2,jb+2) * u(i,pb+2,k)
        tmp3 = tmp3 + Q(pb+2,jb+3) * u(i,pb+2,k)

        tmp0 = tmp0 + Q(pb+3,jb  ) * u(i,pb+3,k)
        tmp1 = tmp1 + Q(pb+3,jb+1) * u(i,pb+3,k)
        tmp2 = tmp2 + Q(pb+3,jb+2) * u(i,pb+3,k)
        tmp3 = tmp3 + Q(pb+3,jb+3) * u(i,pb+3,k)

        v(i,jb  ,k) = alpha * tmp0 + v(i,jb  ,k)
        v(i,jb+1,k) = alpha * tmp1 + v(i,jb+1,k)
        v(i,jb+2,k) = alpha * tmp2 + v(i,jb+2,k)
        v(i,jb+3,k) = alpha * tmp3 + v(i,jb+3,k)

      end do
      end do
    end do
  end do

#endif

#if _NQ_ == _NQ_T4_ + 1

  ! 1:nq/4*4, nq
  ! -----------------------------------------------------------

  pb = NQ_T4 + 1

  !$acc loop collapse(3) independent vector
  do k  = 1, _NQ_
    do jb = 1, NQ_T4, 4 
    do i  = 1, _NQ_

      v(i,jb  ,k) = alpha * Q(pb,jb  ) * u(i,pb,k) + v(i,jb  ,k)
      v(i,jb+1,k) = alpha * Q(pb,jb+1) * u(i,pb,k) + v(i,jb+1,k)
      v(i,jb+2,k) = alpha * Q(pb,jb+2) * u(i,pb,k) + v(i,jb+2,k)
      v(i,jb+3,k) = alpha * Q(pb,jb+3) * u(i,pb,k) + v(i,jb+3,k)

    end do
    end do
  end do

  ! nq, 1:nq
  ! ---------------------------------------------------------------

  jb = NQ_T4 + 1

  !$acc loop collapse(2) independent vector
  do k = 1, _NQ_
    do p = 1, _NQ_
      do i = 1, _NQ_

        v(i,jb,k) = alpha * Q(p,jb) * u(i,p,k) + v(i,jb,k)

      end do
    end do
  end do


#elif _NQ_ == _NQ_T4_ + 2

  ! 1:nq/4*4, nq-2:nq ----------------------------------------------------

  pb = NQ_T4 + 1

  !$acc loop collapse(3) independent vector
  do k = 1, _NQ_
    do jb = 1, NQ_T4, 4  
    do i  = 1, _NQ_

      tmp0 =        Q(pb  ,jb  ) * u(i,pb  ,k)
      tmp1 =        Q(pb  ,jb+1) * u(i,pb  ,k)
      tmp2 =        Q(pb  ,jb+2) * u(i,pb  ,k)
      tmp3 =        Q(pb  ,jb+3) * u(i,pb  ,k)

      tmp0 = tmp0 + Q(pb+1,jb  ) * u(i,pb+1,k)
      tmp1 = tmp1 + Q(pb+1,jb+1) * u(i,pb+1,k)
      tmp2 = tmp2 + Q(pb+1,jb+2) * u(i,pb+1,k)
      tmp3 = tmp3 + Q(pb+1,jb+3) * u(i,pb+1,k)

      v(i,jb  ,k) = alpha * tmp0 + v(i,jb  ,k)
      v(i,jb+1,k) = alpha * tmp1 + v(i,jb+1,k)
      v(i,jb+2,k) = alpha * tmp2 + v(i,jb+2,k)
      v(i,jb+3,k) = alpha * tmp3 + v(i,jb+3,k)

    end do
    end do
  end do

  ! nq-2:nq, 1:nq --------------------------------------------------------

  jb = NQ_T4 + 1

  !$acc loop collapse(2) independent vector
  do k = 1, _NQ_
    do p = 1, _NQ_
    do j = jb, jb+1  
      do i = 1, _NQ_

        v(i,j,k) = alpha * Q(p,j) * u(i,p,k) + v(i,j,k)

      end do
    end do
    end do
  end do

#elif _NQ_ == _NQ_T4_ + 3

 ! 1:nq/4*4, nq-2:nq ----------------------------------------------------

  pb = NQ_T4 + 1

  !$acc loop collapse(3) independent vector
  do k  = 1, _NQ_
    do jb = 1, NQ_T4, 4  
    do i  = 1, _NQ_

      tmp0 =        Q(pb  ,jb  ) * u(i,pb  ,k)
      tmp1 =        Q(pb  ,jb+1) * u(i,pb  ,k)
      tmp2 =        Q(pb  ,jb+2) * u(i,pb  ,k)
      tmp3 =        Q(pb  ,jb+3) * u(i,pb  ,k)

      tmp0 = tmp0 + Q(pb+1,jb  ) * u(i,pb+1,k)
      tmp1 = tmp1 + Q(pb+1,jb+1) * u(i,pb+1,k)
      tmp2 = tmp2 + Q(pb+1,jb+2) * u(i,pb+1,k)
      tmp3 = tmp3 + Q(pb+1,jb+3) * u(i,pb+1,k)

      tmp0 = tmp0 + Q(pb+2,jb  ) * u(i,pb+2,k)
      tmp1 = tmp1 + Q(pb+2,jb+1) * u(i,pb+2,k)
      tmp2 = tmp2 + Q(pb+2,jb+2) * u(i,pb+2,k)
      tmp3 = tmp3 + Q(pb+2,jb+3) * u(i,pb+2,k)

      v(i,jb  ,k) = alpha * tmp0 + v(i,jb  ,k)
      v(i,jb+1,k) = alpha * tmp1 + v(i,jb+1,k)
      v(i,jb+2,k) = alpha * tmp2 + v(i,jb+2,k)
      v(i,jb+3,k) = alpha * tmp3 + v(i,jb+3,k)

    end do
    end do
  end do

  ! nq-2:nq, 1:nq ---------------

  jb = NQ_T4 + 1

  !$acc loop collapse(2) independent vector
  do k = 1, _NQ_
    do p = 1, _NQ_
    do j = jb, jb+2  
      do i = 1, _NQ_

        v(i,j,k) = alpha * Q(p,j) * u(i,p,k) + v(i,j,k)

      end do
    end do
    end do
  end do

#endif

end subroutine PROC(IxQtxI__,_NQ_)

#undef _NQ_T4
