!-------------------------------------------------------------------------------
!> Performs `z1 = IxIxA u`
!>
!>   * k: simple loop (l)
!>   * j: simple loop (l)
!>   * i: simple loop (l)
!>   * p: simple loop (l)

subroutine SubOp_1(At, u, z1)
  !$acc routine vector
  real(RNP), intent(in)  :: At(__NA2__,__NA1__)         !< transpose of A
  real(RNP), intent(in)  :: u (__NA2__,__NA2__,__NA2__) !< u
  real(RNP), intent(out) :: z1(__NA1__,__NA2__,__NA2__) !< intermediate result

  real(RNP) :: tmp
  integer   :: i, j, k, p

  !$acc loop collapse(3) independent vector
  do k = 1, __NA2__
  do j = 1, __NA2__
    !DIR$ SIMD
    do i = 1, __NA1__
      tmp = 0
      do p = 1, __NA2__
        tmp = tmp + At(p,i) * u(p,j,k)
      end do
      z1(i,j,k) = tmp
    end do
  end do
  end do

end subroutine SubOp_1
