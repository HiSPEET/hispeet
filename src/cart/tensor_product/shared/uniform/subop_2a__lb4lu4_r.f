!-------------------------------------------------------------------------------
!> Performs `v += g2 IxAxI u`
!>
!>   * j: blocked with length of 4 (b4)
!>   * i: simple loop (l)
!>   * p: unrolled and jammed with length of 4 (u4)
!>   * explicit remainder handling (r)

subroutine SubOp_2a(g2, At, u, v)
  !$acc routine vector
  real(RNP), intent(in)    :: g2                       !< metric factor
  real(RNP), intent(in)    :: At(__NA__,__NA__)        !< transpose of A
  real(RNP), intent(in)    :: u(__NA__,__NA__,__NA__)  !< operand
  real(RNP), intent(inout) :: v(__NA__,__NA__,__NA__)  !< result

  real(RNP) :: tmp(0:3)
  integer   :: i, j, k, p, q

#if __NA_T4__ > 0

  !$acc loop independent vector
  do k = 1, __NA__
    !$acc loop seq
    do q = 1, __NA_T4__, 4
      !$acc loop collapse(2) independent vector private(tmp)
      do j = 1, __NA_T4__, 4
      do i = 1, __NA__

        tmp =       At(q  ,j:j+3) * u(i,q  ,k)
        tmp = tmp + At(q+1,j:j+3) * u(i,q+1,k)
        tmp = tmp + At(q+2,j:j+3) * u(i,q+2,k)
        tmp = tmp + At(q+3,j:j+3) * u(i,q+3,k)

        v(i,j:j+3,k) = v(i,j:j+3,k) + g2 * tmp

      end do
      end do
    end do
  end do

#endif

#if __NA__ == __NA_T4__ + 1

  !$acc loop independent
  do k = 1, __NA__

    ! remainder: j = 1:na-1, p = na ...........................................

    !$acc loop collapse(2) vector
    do j = 1, __NA_T4__, 4
    do i = 1, __NA__

      v(i,j:j+3,k) = v(i,j:j+3,k) + g2 * At(__NA__,j:j+3) * u(i,__NA__,k)

    end do
    end do

    ! remainder: j = na, p = 1:na-1 ............................................

    !$acc loop seq
    do q = 1, __NA_T4__, 4
      !$acc loop vector private(tmp)
      do i = 1, __NA__

        tmp(0) =          At(q  ,__NA__) * u(i,q  ,k)
        tmp(0) = tmp(0) + At(q+1,__NA__) * u(i,q+1,k)
        tmp(0) = tmp(0) + At(q+2,__NA__) * u(i,q+2,k)
        tmp(0) = tmp(0) + At(q+3,__NA__) * u(i,q+3,k)

        v(i,__NA__,k) = v(i,__NA__,k) + g2 * tmp(0)

      end do
    end do

  end do

  ! remainder: j = na, p = na ................................................

  !$acc loop vector
  do i = 1, __NA__
    v(i,__NA__,k) = v(i,__NA__,k) + g2 * At(__NA__,__NA__) * u(i,__NA__,k)
  end do

#elif __NA__ == __NA_T4__ + 2

  !$acc loop independent vector
  do k = 1, __NA__

    ! remainder: j = 1:na-2, p = na-1:na .....................................

    q = __NA__ - 1

    !$acc loop collapse(2) vector private(tmp)
    do j = 1, __NA_T4__, 4
    do i = 1, __NA__

      tmp =       At(q  ,j:j+3) * u(i,q  ,k)
      tmp = tmp + At(q+1,j:j+3) * u(i,q+1,k)

      v(i,j:j+3,k) = v(i,j:j+3,k) + g2 * tmp

    end do
    end do

    ! remainder: j = na-1:na, p = 1:na-2 .....................................

    !$acc loop seq
    do q = 1, __NA_T4__, 4
      !$acc loop vector private(tmp)
      do i = 1, __NA__

        tmp(0:1) =            At(q  ,__NA__-1:__NA__) * u(i,q  ,k)
        tmp(0:1) = tmp(0:1) + At(q+1,__NA__-1:__NA__) * u(i,q+1,k)
        tmp(0:1) = tmp(0:1) + At(q+2,__NA__-1:__NA__) * u(i,q+2,k)
        tmp(0:1) = tmp(0:1) + At(q+3,__NA__-1:__NA__) * u(i,q+3,k)

        v(i,__NA__-1:__NA__,k) = v(i,__NA__-1:__NA__,k) + g2 * tmp(0:1)

      end do
    end do

    ! remainder: j = na-1:na, p = na-1:na ....................................

    !$acc loop vector
    do i = 1, __NA__
    do p = __NA__-1, __NA__

      v(i,__NA__-1:__NA__,k) = v(i,__NA__-1:__NA__,k) + g2 * At(p,__NA__-1:__NA__) * u(i,p,k)

    end do
    end do

  end do

#elif __NA__ == __NA_T4__ + 3

  !$acc loop independent vector
  do k = 1, __NA__

    ! remainder: j = 1:na-3, p = na-2:na .....................................

    q = __NA__ - 2

    !$acc loop collapse(2) vector private(tmp)
    do j = 1, __NA_T4__, 4
    do i = 1, __NA__

      tmp =       At(q  ,j:j+3) * u(i,q  ,k)
      tmp = tmp + At(q+1,j:j+3) * u(i,q+1,k)
      tmp = tmp + At(q+2,j:j+3) * u(i,q+2,k)

      v(i,j:j+3,k) = v(i,j:j+3,k) + g2 * tmp

    end do
    end do

    ! remainder: j = na-2:na, p = 1:na-3 .....................................

    !$acc loop seq
    do q = 1, __NA_T4__, 4
      !$acc loop vector private(tmp)
      do i = 1, __NA__

        tmp(0:2) =            At(q  ,__NA__-2:__NA__) * u(i,q  ,k)
        tmp(0:2) = tmp(0:2) + At(q+1,__NA__-2:__NA__) * u(i,q+1,k)
        tmp(0:2) = tmp(0:2) + At(q+2,__NA__-2:__NA__) * u(i,q+2,k)
        tmp(0:2) = tmp(0:2) + At(q+3,__NA__-2:__NA__) * u(i,q+3,k)

        v(i,__NA__-2:__NA__,k) = v(i,__NA__-2:__NA__,k) + g2 * tmp(0:2)

      end do
    end do

    ! remainder: j = na-2:na, p = na-2:na ....................................

    !$acc loop vector
    do i = 1, __NA__
    do p = __NA__-2, __NA__

      v(i,__NA__-2:__NA__,k) = v(i,__NA__-2:__NA__,k) + g2 * At(p,__NA__-2:__NA__) * u(i,p,k)

    end do
    end do

  end do

#endif

end subroutine SubOp_2a
