!-------------------------------------------------------------------------------
!> Performs `z2 = IxAxI z1`
!>
!>   * k: simple loop (l)
!>   * j: unrolled and jammed with length of 4 (u4)
!>   * i: simple loop (l)
!>   * p: blocked with length of 4 (b4)
!>   * explicit remainder handling (r)

subroutine SubOp_2(At, z1, z2)
  !$acc routine vector
  real(RNP), intent(in)  :: At(__NA2__,__NA1__)         !< transpose of A
  real(RNP), intent(in)  :: z1(__NA1__,__NA2__,__NA2__) !< z1 = IxIxA u
  real(RNP), intent(out) :: z2(__NA1__,__NA1__,__NA2__) !< z2

  integer :: i, j, k, p, q

  !$acc loop collapse(3) vector
  do k = 1, __NA2__
  do j = 1, __NA1__
  do i = 1, __NA1__
   z2(i,j,k) = 0
  end do
  end do
  end do

  !$acc loop collapse(2) independent vector
  do k = 1, __NA2__
  do j = 1, __NA1_T4__, 4
    !$acc loop seq
    do q = 1, __NA2_T4__, 4
      !$acc loop independent vector
      do i = 1, __NA1__
        do p = q, q + 3
          z2(i,j  ,k) = z2(i,j  ,k) + At(p,j  ) * z1(i,p,k)
          z2(i,j+1,k) = z2(i,j+1,k) + At(p,j+1) * z1(i,p,k)
          z2(i,j+2,k) = z2(i,j+2,k) + At(p,j+2) * z1(i,p,k)
          z2(i,j+3,k) = z2(i,j+3,k) + At(p,j+3) * z1(i,p,k)
        end do
      end do
    end do
  end do
  end do

#if __NA2__ == __NA2_T4__ + 1

  !$acc loop collapse(3) vector
  do k = 1, __NA2__
  do j = 1, __NA1_T4__, 4
  do i = 1, __NA1__
    z2(i,j  ,k) = z2(i,j  ,k) + At(__NA2__,j  ) * z1(i,__NA2__,k)
    z2(i,j+1,k) = z2(i,j+1,k) + At(__NA2__,j+1) * z1(i,__NA2__,k)
    z2(i,j+2,k) = z2(i,j+2,k) + At(__NA2__,j+2) * z1(i,__NA2__,k)
    z2(i,j+3,k) = z2(i,j+3,k) + At(__NA2__,j+3) * z1(i,__NA2__,k)
  end do
  end do
  end do

#elif __NA2__ == __NA2_T4__ + 2

  !$acc loop collapse(3) vector
  do k = 1, __NA2__
  do j = 1, __NA1_T4__, 4
  do i = 1, __NA1__
    do p = __NA2__-1, __NA2__
      z2(i,j  ,k) = z2(i,j  ,k) + At(p,j  ) * z1(i,p,k)
      z2(i,j+1,k) = z2(i,j+1,k) + At(p,j+1) * z1(i,p,k)
      z2(i,j+2,k) = z2(i,j+2,k) + At(p,j+2) * z1(i,p,k)
      z2(i,j+3,k) = z2(i,j+3,k) + At(p,j+3) * z1(i,p,k)
    end do
  end do
  end do
  end do

#elif __NA2__ == __NA2_T4__ + 3

  !$acc loop collapse(3) vector
  do k = 1, __NA2__
  do j = 1, __NA1_T4__, 4
  do i = 1, __NA1__
    do p = __NA2__-2, __NA2__
      z2(i,j  ,k) = z2(i,j  ,k) + At(p,j  ) * z1(i,p,k)
      z2(i,j+1,k) = z2(i,j+1,k) + At(p,j+1) * z1(i,p,k)
      z2(i,j+2,k) = z2(i,j+2,k) + At(p,j+2) * z1(i,p,k)
      z2(i,j+3,k) = z2(i,j+3,k) + At(p,j+3) * z1(i,p,k)
    end do
  end do
  end do
  end do

#endif

#if __NA1__ == __NA1_T4__ + 1

  !$acc loop collapse(2) vector
  do k = 1, __NA2__
  do i = 1, __NA1__
    do p = 1, __NA2__
      z2(i,__NA1__,k) = z2(i,__NA1__,k) + At(p,__NA1__) * z1(i,p,k)
    end do
  end do
  end do

#elif __NA1__ == __NA1_T4__ + 2

  !$acc loop collapse(2) vector
  do k = 1, __NA2__
  do i = 1, __NA1__
    do p = 1, __NA2__
      z2(i,__NA1__-1,k) = z2(i,__NA1__-1,k) + At(p,__NA1__-1) * z1(i,p,k)
      z2(i,__NA1__  ,k) = z2(i,__NA1__  ,k) + At(p,__NA1__  ) * z1(i,p,k)
    end do
  end do
  end do

#elif __NA1__ == __NA1_T4__ + 3

  !$acc loop collapse(2) vector
  do k = 1, __NA2__
  do i = 1, __NA1__
    do p = 1, __NA2__
      z2(i,__NA1__-2,k) = z2(i,__NA1__-2,k) + At(p,__NA1__-2) * z1(i,p,k)
      z2(i,__NA1__-1,k) = z2(i,__NA1__-1,k) + At(p,__NA1__-1) * z1(i,p,k)
      z2(i,__NA1__  ,k) = z2(i,__NA1__  ,k) + At(p,__NA1__  ) * z1(i,p,k)
    end do
  end do
  end do

#endif

end subroutine SubOp_2
