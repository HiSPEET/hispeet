!-------------------------------------------------------------------------------
!> Performs `v = g AᵀxIxI u` for single element
!>
!>   * k: unrolled and jammed with length of 8 (u8)
!>   * j: joined with i and flattened (f)
!>   * i: joined with j and flattened (f)
!>   * p: simple loop (l)
!>   * explicit remainder handling (r)

subroutine SubOp_3i(g, A, u, v)
  !$acc routine vector
  real(RNP), intent(in)    :: g                   !< metric factor
  real(RNP), intent(in)    :: A(__NA__,__NA__)    !< A
  real(RNP), intent(in)    :: u(__NA__**2,__NA__) !< operand
  real(RNP), intent(inout) :: v(__NA__**2,__NA__) !< result

  real(RNP) :: tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7
  integer   :: ij, k, p

#if __NA_T8__ > 0

  !$acc loop collapse(2) vector
  do k = 1, __NA_T8__, 8
  do ij = 1, __NA__**2
    tmp0 = 0
    tmp1 = 0
    tmp2 = 0
    tmp3 = 0
    tmp4 = 0
    tmp5 = 0
    tmp6 = 0
    tmp7 = 0
    do p = 1, __NA__
      tmp0 = tmp0 + A(p,k  ) * u(ij,p)
      tmp1 = tmp1 + A(p,k+1) * u(ij,p)
      tmp2 = tmp2 + A(p,k+2) * u(ij,p)
      tmp3 = tmp3 + A(p,k+3) * u(ij,p)
      tmp4 = tmp4 + A(p,k+4) * u(ij,p)
      tmp5 = tmp5 + A(p,k+5) * u(ij,p)
      tmp6 = tmp6 + A(p,k+6) * u(ij,p)
      tmp7 = tmp7 + A(p,k+7) * u(ij,p)
    end do
    v(ij,k  ) = g * tmp0
    v(ij,k+1) = g * tmp1
    v(ij,k+2) = g * tmp2
    v(ij,k+3) = g * tmp3
    v(ij,k+4) = g * tmp4
    v(ij,k+5) = g * tmp5
    v(ij,k+6) = g * tmp6
    v(ij,k+7) = g * tmp7
  end do
  end do

#endif

#if __NA__ == __NA_T8__ + 1

  k = __NA_T8__ + 1

  !$acc loop vector
  do ij = 1, __NA__**2
    tmp0 = 0
    do p = 1, __NA__
      tmp0 = tmp0 + A(p,k) * u(ij,p)
    end do
    v(ij,k) = g * tmp0
  end do

#elif __NA__ == __NA_T8__ + 2

  k = __NA_T8__ + 1

  !$acc loop vector
  do ij = 1, __NA__**2
    tmp0 = 0
    tmp1 = 0
    do p = 1, __NA__
      tmp0 = tmp0 + A(p,k  ) * u(ij,p)
      tmp1 = tmp1 + A(p,k+1) * u(ij,p)
    end do
    v(ij,k  ) = g * tmp0
    v(ij,k+1) = g * tmp1
  end do

#elif __NA__ == __NA_T8__ + 3

  k = __NA_T8__ + 1

  !$acc loop vector
  do ij = 1, __NA__**2
    tmp0 = 0
    tmp1 = 0
    tmp2 = 0
    do p = 1, __NA__
      tmp0 = tmp0 + A(p,k  ) * u(ij,p)
      tmp1 = tmp1 + A(p,k+1) * u(ij,p)
      tmp2 = tmp2 + A(p,k+2) * u(ij,p)
    end do
    v(ij,k  ) = g * tmp0
    v(ij,k+1) = g * tmp1
    v(ij,k+2) = g * tmp2
  end do

#elif __NA__ == __NA_T8__ + 4

  k = __NA_T8__ + 1

  !$acc loop vector
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
    v(ij,k  ) = g * tmp0
    v(ij,k+1) = g * tmp1
    v(ij,k+2) = g * tmp2
    v(ij,k+3) = g * tmp3
  end do

#elif __NA__ == __NA_T8__ + 5

  k = __NA_T8__ + 1

  !$acc loop vector
  do ij = 1, __NA__**2
    tmp0 = 0
    tmp1 = 0
    tmp2 = 0
    tmp3 = 0
    tmp4 = 0
    do p = 1, __NA__
      tmp0 = tmp0 + A(p,k  ) * u(ij,p)
      tmp1 = tmp1 + A(p,k+1) * u(ij,p)
      tmp2 = tmp2 + A(p,k+2) * u(ij,p)
      tmp3 = tmp3 + A(p,k+3) * u(ij,p)
      tmp4 = tmp4 + A(p,k+4) * u(ij,p)
    end do
    v(ij,k  ) = g * tmp0
    v(ij,k+1) = g * tmp1
    v(ij,k+2) = g * tmp2
    v(ij,k+3) = g * tmp3
    v(ij,k+4) = g * tmp4
  end do

#elif __NA__ == __NA_T8__ + 6

  k = __NA_T8__ + 1

  !$acc loop vector
  do ij = 1, __NA__**2
    tmp0 = 0
    tmp1 = 0
    tmp2 = 0
    tmp3 = 0
    tmp4 = 0
    tmp5 = 0
    do p = 1, __NA__
      tmp0 = tmp0 + A(p,k  ) * u(ij,p)
      tmp1 = tmp1 + A(p,k+1) * u(ij,p)
      tmp2 = tmp2 + A(p,k+2) * u(ij,p)
      tmp3 = tmp3 + A(p,k+3) * u(ij,p)
      tmp4 = tmp4 + A(p,k+4) * u(ij,p)
      tmp5 = tmp5 + A(p,k+5) * u(ij,p)
    end do
    v(ij,k  ) = g * tmp0
    v(ij,k+1) = g * tmp1
    v(ij,k+2) = g * tmp2
    v(ij,k+3) = g * tmp3
    v(ij,k+4) = g * tmp4
    v(ij,k+5) = g * tmp5
  end do

#elif __NA__ == __NA_T8__ + 7

  k = __NA_T8__ + 1

  !$acc loop vector
  do ij = 1, __NA__**2
    tmp0 = 0
    tmp1 = 0
    tmp2 = 0
    tmp3 = 0
    tmp4 = 0
    tmp5 = 0
    tmp6 = 0
    do p = 1, __NA__
      tmp0 = tmp0 + A(p,k  ) * u(ij,p)
      tmp1 = tmp1 + A(p,k+1) * u(ij,p)
      tmp2 = tmp2 + A(p,k+2) * u(ij,p)
      tmp3 = tmp3 + A(p,k+3) * u(ij,p)
      tmp4 = tmp4 + A(p,k+4) * u(ij,p)
      tmp5 = tmp5 + A(p,k+5) * u(ij,p)
      tmp6 = tmp6 + A(p,k+6) * u(ij,p)
    end do
    v(ij,k  ) = g * tmp0
    v(ij,k+1) = g * tmp1
    v(ij,k+2) = g * tmp2
    v(ij,k+3) = g * tmp3
    v(ij,k+4) = g * tmp4
    v(ij,k+5) = g * tmp5
    v(ij,k+6) = g * tmp6
  end do

#endif

end subroutine SubOp_3i
