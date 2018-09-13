!-------------------------------------------------------------------------------
!> Performs `z2 = IxAxI z1` for na1 = na2 = 4
!>
!>   * k: simple loop (l)
!>   * j: unrolled and jammed with length of 4 (u4)
!>   * i: simple loop (l)
!>   * p: simple loop (l)
!>   * no remainder handling

subroutine SubOp_2(At, z1, z2)
  !$acc routine vector
  real(RNP), intent(in)  :: At(4,4)   !< transpose of A
  real(RNP), intent(in)  :: z1(4,4,4) !< z1 = IxIxA u
  real(RNP), intent(out) :: z2(4,4,4) !< z2

  real(RNP) :: tmp0, tmp1, tmp2, tmp3
  integer   :: i, k, p

  !$acc loop collapse(2) vector
  do k = 1, 4
  do i = 1, 4
    tmp0 = 0
    tmp1 = 0
    tmp2 = 0
    tmp3 = 0
    do p = 1, 4
      tmp0 = tmp0 + At(p,1) * z1(i,p,k)
      tmp1 = tmp1 + At(p,2) * z1(i,p,k)
      tmp2 = tmp2 + At(p,3) * z1(i,p,k)
      tmp3 = tmp3 + At(p,4) * z1(i,p,k)
    end do
    z2(i,1,k) = tmp0
    z2(i,2,k) = tmp1
    z2(i,3,k) = tmp2
    z2(i,4,k) = tmp3
  end do
  end do

end subroutine SubOp_2
