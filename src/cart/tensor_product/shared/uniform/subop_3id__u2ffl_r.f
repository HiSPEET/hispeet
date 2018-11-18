!-------------------------------------------------------------------------------
!> Performs `v = g D AᵀxIxI u` for single element
!>
!>   * k: unrolled and jammed with length of 4 (u4)
!>   * j: joined with i and flattened (f)
!>   * i: joined with j and flattened (f)
!>   * p: simple loop (l)
!>   * explicit remainder handling (r)

subroutine SubOp_3id(g, A, D, u, v)
  !$acc routine vector
  real(RNP), intent(in)  :: g                   !< metric factor
  real(RNP), intent(in)  :: A(__NA__,__NA__)    !< A
  real(RNP), intent(in)  :: D(__NA__**2,__NA__) !< diagonal operator
  real(RNP), intent(in)  :: u(__NA__**2,__NA__) !< operand
  real(RNP), intent(out) :: v(__NA__**2,__NA__) !< result

  real(RNP) :: tmp0, tmp1
  integer   :: ij, k, p

  !$acc loop collapse(2) vector
  do k = 1, __NA_T2__, 2
  do ij = 1, __NA__**2
    tmp0 = 0
    tmp1 = 0
    do p = 1, __NA__
      tmp0 = tmp0 + A(p,k  ) * u(ij,p)
      tmp1 = tmp1 + A(p,k+1) * u(ij,p)
    end do
    v(ij,k  ) = g * D(ij,k  ) * tmp0
    v(ij,k+1) = g * D(ij,k+1) * tmp1
  end do
  end do

#if __NA__ == __NA_T2__ + 1

  ! remainder: k = na ........................................................

  !$acc loop vector
  do ij = 1, __NA__**2
    tmp0 = 0
    do p = 1, __NA__
      tmp0 = tmp0 + A(p,__NA__) * u(ij,p)
    end do
    v(ij,__NA__) = g * D(ij,__NA__) * tmp0
  end do

#endif

end subroutine SubOp_3id
