!-------------------------------------------------------------------------------
!> Performs `v += g2 IxAxI u`
!>
!>   * k: simple loop (l)
!>   * j: unrolled and jammed with length of 4 (u4)
!>   * i: simple loop (l)
!>   * p: blocked with length of 4 (b4)
!>   * explicit remainder handling (r)

subroutine SubOp_2a(g2, At, u, v)
  !$acc routine vector
  real(RNP), intent(in)    :: g2                       !< metric factor
  real(RNP), intent(in)    :: At(__NA__,__NA__)        !< transpose of A
  real(RNP), intent(in)    :: u(__NA__,__NA__,__NA__)  !< operand
  real(RNP), intent(inout) :: v(__NA__,__NA__,__NA__)  !< result

  real(RNP) :: tmp0, tmp1, tmp2, tmp3
  integer   :: i, j, k, p, q

#if __NA_T4__ > 0

  !$acc loop collapse(2) independent vector
  do k = 1, __NA__
  do j = 1, __NA_T4__, 4
    !$acc loop seq
    do q = 1, __NA_T4__, 4
      !$acc loop independent vector
      do i = 1, __NA__
        tmp0 = 0
        tmp1 = 0
        tmp2 = 0
        tmp3 = 0
        do p = q, q+3
          tmp0 = tmp0 + At(p,j  ) * u(i,p,k)
          tmp1 = tmp1 + At(p,j+1) * u(i,p,k)
          tmp2 = tmp2 + At(p,j+2) * u(i,p,k)
          tmp3 = tmp3 + At(p,j+3) * u(i,p,k)
        end do
        v(i,j  ,k) = v(i,j  ,k) + g2 * tmp0
        v(i,j+1,k) = v(i,j+1,k) + g2 * tmp1
        v(i,j+2,k) = v(i,j+2,k) + g2 * tmp2
        v(i,j+3,k) = v(i,j+3,k) + g2 * tmp3
      end do
      end do
    end do
  end do

#endif

#if __NA__ == __NA_T4__ + 1

  !$acc loop independent vector
  do k = 1, __NA__

    ! remainder: j = 1: na-1, p = na .........................................

    p = __NA__

    !$acc loop collapse(2) vector
    do j = 1, __NA_T4__, 4
    do i = 1, __NA__
      v(i,j  ,k) = v(i,j  ,k) + g2 * At(p,j  ) * u(i,p,k)
      v(i,j+1,k) = v(i,j+1,k) + g2 * At(p,j+1) * u(i,p,k)
      v(i,j+2,k) = v(i,j+2,k) + g2 * At(p,j+2) * u(i,p,k)
      v(i,j+3,k) = v(i,j+3,k) + g2 * At(p,j+3) * u(i,p,k)
    end do
    end do

    ! remainder: j = na, p = 1:na  ...........................................

    !$acc loop vector
    do i = 1, __NA__
      tmp0 = 0
      do p = 1, __NA__
        tmp0 = tmp0 + At(p,__NA__) * u(i,p,k)
      end do
      v(i,__NA__,k) = v(i,__NA__,k) + g2 * tmp0
    end do

  end do

#elif __NA__ == __NA_T4__ + 2

  !$acc loop independent vector
  do k = 1, __NA__

    ! remainder: j = 1:na-2, p = na-1:na .....................................

    !$acc loop collapse(2) vector
    do j = 1, __NA_T4__, 4
    do i = 1, __NA__
      do p = __NA__-1, __NA__
        v(i,j  ,k) = v(i,j  ,k) + g2 * At(p,j  ) * u(i,p,k)
        v(i,j+1,k) = v(i,j+1,k) + g2 * At(p,j+1) * u(i,p,k)
        v(i,j+2,k) = v(i,j+2,k) + g2 * At(p,j+2) * u(i,p,k)
        v(i,j+3,k) = v(i,j+3,k) + g2 * At(p,j+3) * u(i,p,k)
      end do
    end do
    end do

    ! remainder: j = na-1:na, p = 1:na .......................................

    !$acc loop vector
    do i = 1, __NA__
      tmp0 = 0
      tmp1 = 0
      do p = 1, __NA__
        tmp0 = tmp0 + At(p,__NA__-1) * u(i,p,k)
        tmp1 = tmp1 + At(p,__NA__  ) * u(i,p,k)
      end do
      v(i,__NA__-1,k) = v(i,__NA__-1,k) + g2 * tmp0
      v(i,__NA__  ,k) = v(i,__NA__  ,k) + g2 * tmp1
    end do

  end do

#elif __NA__ == __NA_T4__ + 3

  !$acc loop independent vector
  do k = 1, __NA__

    ! remainder: j = 1:na-3, p = na-2:na .....................................

    !$acc loop collapse(2) vector
    do j = 1, __NA_T4__, 4
    do i = 1, __NA__
      do p = __NA__-2, __NA__
        v(i,j  ,k) = v(i,j  ,k) + g2 * At(p,j  ) * u(i,p,k)
        v(i,j+1,k) = v(i,j+1,k) + g2 * At(p,j+1) * u(i,p,k)
        v(i,j+2,k) = v(i,j+2,k) + g2 * At(p,j+2) * u(i,p,k)
        v(i,j+3,k) = v(i,j+3,k) + g2 * At(p,j+3) * u(i,p,k)
      end do
    end do
    end do

    ! remainder: j = na-2:na, p = 1:na .......................................

    !$acc loop vector
    do i = 1, __NA__
      tmp0 = 0
      tmp1 = 0
      tmp2 = 0
      do p = 1, __NA__
        tmp0 = tmp0 + At(p,__NA__-2) * u(i,p,k)
        tmp1 = tmp1 + At(p,__NA__-1) * u(i,p,k)
        tmp2 = tmp2 + At(p,__NA__  ) * u(i,p,k)
      end do
      v(i,__NA__-2,k) = v(i,__NA__-2,k) + g2 * tmp0
      v(i,__NA__-1,k) = v(i,__NA__-1,k) + g2 * tmp1
      v(i,__NA__  ,k) = v(i,__NA__  ,k) + g2 * tmp2
    end do

  end do

#endif

end subroutine SubOp_2a
