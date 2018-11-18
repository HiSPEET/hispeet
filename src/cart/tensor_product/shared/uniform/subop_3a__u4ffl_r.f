!-------------------------------------------------------------------------------
!> Performs `v += g AᵀxIxI u` for single element
!>
!>   * k: unrolled and jammed with length of 4 (u4)
!>   * j: joined with i and flattened (f)
!>   * i: joined with j and flattened (f)
!>   * p: simple loop (l)
!>   * explicit remainder handling (r)

subroutine SubOp_3a(g, A, u, v)
  !$acc routine vector
  real(RNP), intent(in)    :: g                   !< metric factor
  real(RNP), intent(in)    :: A(__NA__,__NA__)    !< A
  real(RNP), intent(in)    :: u(__NA__**2,__NA__) !< operand
  real(RNP), intent(inout) :: v(__NA__**2,__NA__) !< result

  real(RNP) :: tmp0, tmp1, tmp2, tmp3
  integer   :: ij, k, p

#if __NA_T4__ > 0

  !$acc loop collapse(2) vector
  do k = 1, __NA_T4__, 4
  do ij = 1, __NA__**2
    tmp0 = 0
    tmp1 = 0
    tmp2 = 0
    tmp3 = 0
    do p = 1, __NA__
      tmp0 = tmp0 + A(p,k  ) * u(ij,p)
      tmp1 = tmp1 + A(p,k+1) * u(ij,p)
      tmp2 = tmp2 + A(p,k+2) * u(ij,p)
      tmp3 = tmp3 + A(p,k+3) * u(ij,p)
    end do
    v(ij,k  ) = v(ij,k  ) + g * tmp0
    v(ij,k+1) = v(ij,k+1) + g * tmp1
    v(ij,k+2) = v(ij,k+2) + g * tmp2
    v(ij,k+3) = v(ij,k+3) + g * tmp3
  end do
  end do

#endif

#if __NA__ == __NA_T4__ + 1

  ! remainder: k = na, p = 1:na ..............................................

  !$acc loop vector
  do ij = 1, __NA__**2
    tmp0 = 0
    do p = 1, __NA__
      tmp0 = tmp0 + A(p,__NA__) * u(ij,p)
    end do
    v(ij,__NA__) = v(ij,__NA__) + g * tmp0
  end do

#elif __NA__ == __NA_T4__ + 2

  ! remainder: k = na-1:na, p = 1:na .........................................

  !$acc loop vector
  do ij = 1, __NA__**2
    tmp0 = 0
    tmp1 = 0
    do p = 1, __NA__
      tmp0 = tmp0 + A(p,__NA__-1) * u(ij,p)
      tmp1 = tmp1 + A(p,__NA__  ) * u(ij,p)
    end do
    v(ij,__NA__-1) = v(ij,__NA__-1) + g * tmp0
    v(ij,__NA__  ) = v(ij,__NA__  ) + g * tmp1
  end do

#elif __NA__ == __NA_T4__ + 3

  !$acc loop vector
  do ij = 1, __NA__**2
    tmp0 = 0
    tmp1 = 0
    tmp2 = 0
    do p = 1, __NA__
      tmp0 = tmp0 + A(p,__NA__-2) * u(ij,p)
      tmp1 = tmp1 + A(p,__NA__-1) * u(ij,p)
      tmp2 = tmp2 + A(p,__NA__  ) * u(ij,p)
    end do
    v(ij,__NA__-2) = v(ij,__NA__-2) + g * tmp0
    v(ij,__NA__-1) = v(ij,__NA__-1) + g * tmp1
    v(ij,__NA__  ) = v(ij,__NA__  ) + g * tmp2
  end do

#endif

end subroutine SubOp_3a
