!-------------------------------------------------------------------------------
!> Performs `v =(d₀IxIxI + d₁IxIxV₁ + d₂IxIxV₂ + d₃IxIxV₃)⁻¹ (IxIxA') u`
!>
!>   * k: simple loop (l)
!>   * j: unrolled and jammed with length of 4 (u4)
!>   * i: unrolled and jammed with length of 2 (u2)
!>   * p: simple loop (l)
!>   * explicit remainder handling (r)

subroutine SubOp_1d(d0, d1, d2, d3, V1, V2, V3, A, u, v)
  !$acc routine vector
  real(RNP), intent(in)  :: d0                       !< λ Δx₁ Δx₂ Δx₃
  real(RNP), intent(in)  :: d1                       !< ν Δx₂ Δx₃ / Δx₁
  real(RNP), intent(in)  :: d2                       !< ν Δx₃ Δx₁ / Δx₂
  real(RNP), intent(in)  :: d3                       !< ν Δx₁ Δx₂ / Δx₃
  real(RNP), intent(in)  :: V1(__NA__)               !< 1D eigenvalues in dir 1
  real(RNP), intent(in)  :: V2(__NA__)               !< 1D eigenvalues in dir 2
  real(RNP), intent(in)  :: V3(__NA__)               !< 1D eigenvalues in dir 3
  real(RNP), intent(in)  :: A(__NA__,__NA__)         !< A
  real(RNP), intent(in)  :: u(__NA__,__NA__,__NA__)  !< operand
  real(RNP), intent(out) :: v(__NA__,__NA__,__NA__)  !< result

  real(RNP) :: tmp00, tmp10
  real(RNP) :: tmp01, tmp11
  real(RNP) :: tmp02, tmp12
  real(RNP) :: tmp03, tmp13
  integer   :: i, j, k, p

#if __NA_T4__ > 0

  !$acc loop collapse(3) vector
  do k = 1, __NA__
  do j = 1, __NA_T4__, 4
  do i = 1, __NA_T2__, 2
    tmp00 = 0
    tmp10 = 0
    tmp01 = 0
    tmp11 = 0
    tmp02 = 0
    tmp12 = 0
    tmp03 = 0
    tmp13 = 0
    do p = 1, __NA__
      tmp00 = tmp00 + A(p,i  ) * u(p,j  ,k)
      tmp10 = tmp10 + A(p,i+1) * u(p,j  ,k)
      tmp01 = tmp01 + A(p,i  ) * u(p,j+1,k)
      tmp11 = tmp11 + A(p,i+1) * u(p,j+1,k)
      tmp02 = tmp02 + A(p,i  ) * u(p,j+2,k)
      tmp12 = tmp12 + A(p,i+1) * u(p,j+2,k)
      tmp03 = tmp03 + A(p,i  ) * u(p,j+3,k)
      tmp13 = tmp13 + A(p,i+1) * u(p,j+3,k)
    end do
    v(i  ,j  ,k) = tmp00
    v(i+1,j  ,k) = tmp10
    v(i  ,j+1,k) = tmp01
    v(i+1,j+1,k) = tmp11
    v(i  ,j+2,k) = tmp02
    v(i+1,j+2,k) = tmp12
    v(i  ,j+3,k) = tmp03
    v(i+1,j+3,k) = tmp13
  end do
  end do
  end do

#endif

#if __NA__ == __NA_T4__ + 1

  ! remainder: i = 1:na-1, j = na ............................................

  j = __NA__

  !$acc loop collapse(2) vector
  do k = 1, __NA__
  do i = 1, __NA_T2__, 2
    tmp00 = 0
    tmp10 = 0
    do p = 1, __NA__
      tmp00 = tmp00 + A(p,i  ) * u(p,j,k)
      tmp10 = tmp10 + A(p,i+1) * u(p,j,k)
    end do
    v(i  ,j,k) = tmp00
    v(i+1,j,k) = tmp10
  end do
  end do

  ! remainder: i = na, j = 1:na ..............................................

  i = __NA__

  !$acc loop collapse(2) vector
  do k = 1, __NA__
  do j = 1, __NA__
    tmp00 = 0
    do p = 1, __NA__
      tmp00 = tmp00 + A(p,i) * u(p,j,k)
    end do
    v(i,j,k) = tmp00
  end do
  end do

#elif __NA__ == __NA_T4__ + 2

  ! remainder: i = 1:na, j = na-1:na .........................................

  j = __NA__ - 1

  !$acc loop collapse(2) vector
  do k = 1, __NA__
  do i = 1, __NA_T2__, 2
    tmp00 = 0
    tmp10 = 0
    tmp01 = 0
    tmp11 = 0
    do p = 1, __NA__
      tmp00 = tmp00 + A(p,i  ) * u(p,j  ,k)
      tmp10 = tmp10 + A(p,i+1) * u(p,j  ,k)
      tmp01 = tmp01 + A(p,i  ) * u(p,j+1,k)
      tmp11 = tmp11 + A(p,i+1) * u(p,j+1,k)
    end do
    v(i  ,j  ,k) = tmp00
    v(i+1,j  ,k) = tmp10
    v(i  ,j+1,k) = tmp01
    v(i+1,j+1,k) = tmp11
  end do
  end do

#elif __NA__ == __NA_T4__ + 3

  ! remainder: i = 1:na-1, j = na-2:na .......................................

  j = __NA__ - 2

  !$acc loop collapse(2) vector
  do k = 1, __NA__
  do i = 1, __NA_T2__, 2
     tmp00 = 0
     tmp10 = 0
     tmp01 = 0
     tmp11 = 0
     tmp02 = 0
     tmp12 = 0
     do p = 1, __NA__
       tmp00 = tmp00 + A(p,i  ) * u(p,j  ,k)
       tmp10 = tmp10 + A(p,i+1) * u(p,j  ,k)
       tmp01 = tmp01 + A(p,i  ) * u(p,j+1,k)
       tmp11 = tmp11 + A(p,i+1) * u(p,j+1,k)
       tmp02 = tmp02 + A(p,i  ) * u(p,j+2,k)
       tmp12 = tmp12 + A(p,i+1) * u(p,j+2,k)
     end do
     v(i  ,j  ,k) = tmp00
     v(i+1,j  ,k) = tmp10
     v(i  ,j+1,k) = tmp01
     v(i+1,j+1,k) = tmp11
     v(i  ,j+2,k) = tmp02
     v(i+1,j+2,k) = tmp12
  end do
  end do

  ! remainder: i = na, j = 1:na ..............................................

  i = __NA__

  !$acc loop collapse(2) vector
  do k = 1, __NA__
  do j = 1, __NA__
    tmp00 = 0
    do p = 1, __NA__
      tmp00 = tmp00 + A(p,i) * u(p,j,k)
    end do
    v(i,j,k) = tmp00
  end do
  end do

#endif

  !$acc loop collapse(3) vector
  do k = 1, __NA__
  do j = 1, __NA__
  do i = 1, __NA__
    v(i,j,k) = v(i,j,k) / (d0 + d1*V1(i) + d2*V2(j) + d3*V3(k))
  end do
  end do
  end do

end subroutine SubOp_1d
