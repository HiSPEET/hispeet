!-------------------------------------------------------------------------------
!> Performs `z2 = IxAxI z1`
!>
!>   * k: simple loop (l)
!>   * j: simple loop (l)
!>   * i: simple loop (l)
!>   * p: simple loop (l)

subroutine SubOp_2(At, z1, z2)
  !$acc routine vector
  real(RNP), intent(in)  :: At(__NA2__,__NA1__)         !< transpose of A
  real(RNP), intent(in)  :: z1(__NA1__,__NA2__,__NA2__) !< z1 = IxIxA u
  real(RNP), intent(out) :: z2(__NA1__,__NA1__,__NA2__) !< z2

  real(RNP) :: tmp
  integer   :: i, j, k, p

  !$acc loop collapse(3) independent vector
  do k = 1, __NA2__
  do j = 1, __NA1__
    !DIR$ SIMD
    do i = 1, __NA1__
      tmp = 0
      do p = 1, __NA2__
        tmp = tmp + At(p,j) * z1(i,p,k)
      end do
      z2(i,j,k) = tmp
    end do
  end do
  end do

end subroutine SubOp_2
