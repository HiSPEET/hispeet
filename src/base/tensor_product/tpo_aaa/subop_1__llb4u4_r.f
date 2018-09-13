!-------------------------------------------------------------------------------
!> Performs `z1 = IxIxA u`
!>
!>   * k: simple loop (l)
!>   * j: simple loop (l)
!>   * i: blocked with length of 4 (b4)
!>   * p: unrolled and jammed with length of 4 (u4)
!>   * explicit remainder handling (r)

subroutine SubOp_1(At, u, z1)
  !$acc routine vector
  real(RNP), intent(in)  :: At(__NA2__,__NA1__)         !< transpose of A
  real(RNP), intent(in)  :: u (__NA2__,__NA2__,__NA2__) !< u
  real(RNP), intent(out) :: z1(__NA1__,__NA2__,__NA2__) !< intermediate result

  real(RNP) :: tmp(0:3)
  integer   :: i, j, k, p, q

  !$acc loop collapse(3) vector
  do k = 1, __NA2__
  do j = 1, __NA2__
  do i = 1, __NA1__
    z1(i,j,k) = 0
  end do
  end do
  end do

  !$acc loop independent vector
  do k = 1, __NA2__
    !$acc loop seq
    do q = 1, __NA2_T4__, 4
      !$acc loop collapse(2) independent vector private(tmp)
      do j = 1, __NA2__
      do i = 1, __NA1_T4__, 4

        tmp =       At(q  ,i:i+3) * u(q  ,j,k)
        tmp = tmp + At(q+1,i:i+3) * u(q+1,j,k)
        tmp = tmp + At(q+2,i:i+3) * u(q+2,j,k)
        tmp = tmp + At(q+3,i:i+3) * u(q+3,j,k)

        z1(i:i+3,j,k) = z1(i:i+3,j,k) + tmp

      end do
      end do
    end do
  end do

#if __NA2__ == __NA2_T4__ + 1

  ! 1:na1/4*4, na2 -------------------------------------------------------------

  q = __NA2_T4__ + 1

  !$acc loop collapse(3) independent vector
  do k = 1, __NA2__
    do j = 1, __NA2__
    do i = 1, __NA1_T4__, 4

      z1(i:i+3,j,k) = z1(i:i+3,j,k) + At(q,i:i+3) * u(q,j,k)

    end do
    end do
  end do

#elif __NA2__ == __NA2_T4__ + 2

  ! 1:na1/4*4, na2-1:na2 -------------------------------------------------------

  q = __NA2_T4__ + 1

  !$acc loop collapse(3) independent vector private(tmp)
  do k = 1, __NA2__
    do j = 1, __NA2__
    do i = 1, __NA1_T4__, 4

      tmp =       At(q  ,i:i+3) * u(q  ,j,k)
      tmp = tmp + At(q+1,i:i+3) * u(q+1,j,k)

      z1(i:i+3,j,k) = z1(i:i+3,j,k) + tmp

    end do
    end do
  end do

#elif __NA2__ == __NA2_T4__ + 3

  ! 1:na1/4*4, na2-2:na2 -------------------------------------------------------

  q = __NA2_T4__ + 1

  !$acc loop collapse(3) independent vector private(tmp)
  do k = 1, __NA2__
    do j = 1, __NA2__
    do i = 1, __NA1_T4__, 4

      tmp =       At(q  ,i:i+3) * u(q  ,j,k)
      tmp = tmp + At(q+1,i:i+3) * u(q+1,j,k)
      tmp = tmp + At(q+2,i:i+3) * u(q+2,j,k)

      z1(i:i+3,j,k) = z1(i:i+3,j,k) + tmp

    end do
    end do
  end do

#endif

#if __NA1__ == __NA1_T4__ + 1

  ! na1, 1:na2 -----------------------------------------------------------------

  i = __NA1_T4__ + 1

  !$acc loop collapse(2) independent vector
  do k = 1, __NA2__
  do j = 1, __NA2__
    do p = 1, __NA2__

      z1(i,j,k) = z1(i,j,k) + At(p,i) * u(p,j,k)

    end do
  end do
  end do

#elif __NA1__ == __NA1_T4__ + 2

  ! na1-1:na1, 1:na2 -----------------------------------------------------------

  i = __NA1_T4__ + 1

  !$acc loop collapse(2) independent vector
  do k = 1, __NA2__
  do j = 1, __NA2__
    do p = 1, __NA2__

      z1(i:i+1,j,k) = z1(i:i+1,j,k) + At(p,i:i+1) * u(p,j,k)

    end do
  end do
  end do

#elif __NA1__ == __NA1_T4__ + 3

  ! na1-2:na1, 1:na2 -----------------------------------------------------------

  i = __NA1_T4__ + 1

  !$acc loop collapse(2) independent vector
  do k = 1, __NA2__
  do j = 1, __NA2__
    do p = 1, __NA2__

      z1(i:i+2,j,k) = z1(i:i+2,j,k) + At(p,i:i+2) * u(p,j,k)

    end do
  end do
  end do

#endif

end subroutine SubOp_1
