!-------------------------------------------------------------------------------
!> Performs `v =(d₀IxIxI + d₁IxIxV₁ + d₂IxIxV₂ + d₃IxIxV₃)⁻¹ (IxIxA') u`
!>
!>   * k: simple loop (l)
!>   * j: simple loop (l)
!>   * i: blocked with length of 4 (b4)
!>   * p: unrolled and jammed with length of 4 (u4)
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

  real(RNP) :: tmp(0:3)
  integer   :: i, j, k, p, q

  !$acc loop collapse(3) vector
  do k = 1, __NA__
  do j = 1, __NA__
  do i = 1, __NA__
    v(i,j,k) = 0
  end do
  end do
  end do

#if __NA_T4__ > 0

  !$acc loop independent vector
  do k = 1, __NA__
    !$acc loop seq
    do q = 1, __NA_T4__, 4
      !$acc loop collapse(2) vector private(tmp)
      do j = 1, __NA__
      do i = 1, __NA_T4__, 4

        tmp =       A(q  ,i:i+3) * u(q  ,j,k)
        tmp = tmp + A(q+1,i:i+3) * u(q+1,j,k)
        tmp = tmp + A(q+2,i:i+3) * u(q+2,j,k)
        tmp = tmp + A(q+3,i:i+3) * u(q+3,j,k)

        v(i:i+3,j,k) = v(i:i+3,j,k) + tmp

      end do
      end do
    end do
  end do

#endif

#if __NA__ == __NA_T4__ + 1

  ! 1:na/4*4, na-2:na ----------------------------------------------------------

  q = __NA_T4__ + 1

  !$acc loop collapse(3) independent vector private(tmp)
  do k = 1, __NA__
    do j = 1, __NA__
    do i = 1, __NA_T4__, 4

      v(i:i+3,j,k) = v(i:i+3,j,k) + A(q,i:i+3) * u(q,j,k)

    end do
    end do
  end do

  ! na, 1:na -------------------------------------------------------------------

  i = __NA_T4__ + 1

  !$acc loop collapse(2) independent vector
  do k = 1, __NA__
  do j = 1, __NA__
    do p = 1, __NA__

      v(i,j,k) = v(i,j,k) + A(p,i) * u(p,j,k)

    end do
  end do
  end do

#elif __NA__ == __NA_T4__ + 2

  ! 1:na/4*4, na-2:na ----------------------------------------------------------

  q = __NA_T4__ + 1

  !$acc loop collapse(3) independent vector private(tmp)
  do k = 1, __NA__
    do j = 1, __NA__
    do i = 1, __NA_T4__, 4

      tmp =       A(q  ,i:i+3) * u(q  ,j,k)
      tmp = tmp + A(q+1,i:i+3) * u(q+1,j,k)

      v(i:i+3,j,k) = v(i:i+3,j,k) + tmp

    end do
    end do
  end do

  ! na-2:na, 1:na --------------------------------------------------------------

  i = __NA_T4__ + 1

  !$acc loop collapse(2) independent vector
  do k = 1, __NA__
  do j = 1, __NA__
    do p = 1, __NA__

      v(i:i+1,j,k) = v(i:i+1,j,k) + A(p,i:i+1) * u(p,j,k)

    end do
  end do
  end do

#elif __NA__ == __NA_T4__ + 3

  ! 1:na/4*4, na-2:na ----------------------------------------------------------

  q = __NA_T4__ + 1

  !$acc loop collapse(3) independent vector private(tmp)
  do k = 1, __NA__
    do j = 1, __NA__
    do i = 1, __NA_T4__, 4

      tmp =       A(q  ,i:i+3) * u(q  ,j,k)
      tmp = tmp + A(q+1,i:i+3) * u(q+1,j,k)
      tmp = tmp + A(q+2,i:i+3) * u(q+2,j,k)

      v(i:i+3,j,k) = v(i:i+3,j,k) + tmp

    end do
    end do
  end do

  ! na-2:na, 1:na --------------------------------------------------------------

  i = __NA_T4__ + 1

  !$acc loop collapse(2) independent vector
  do k = 1, __NA__
  do j = 1, __NA__
    do p = 1, __NA__

      v(i:i+2,j,k) = v(i:i+2,j,k) + A(p,i:i+2) * u(p,j,k)

    end do
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
