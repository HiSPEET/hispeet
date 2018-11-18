!-------------------------------------------------------------------------------
!> Performs `v = g D IxIxAᵀ u`
!>
!>   * k: simple loop (l)
!>   * j: simple loop (l)
!>   * i: blocked with length of 2 (b2)
!>   * p: unrolled and jammed with length of 2 (u2)
!>   * explicit remainder handling (r)

subroutine SubOp_1id(g, A, D, u, v)
  !$acc routine vector
  real(RNP), intent(in)  :: g                       !< metric factor
  real(RNP), intent(in)  :: A(__NA__,__NA__)        !< A
  real(RNP), intent(in)  :: D(__NA__,__NA__,__NA__) !< diagonal operator
  real(RNP), intent(in)  :: u(__NA__,__NA__,__NA__) !< operand
  real(RNP), intent(out) :: v(__NA__,__NA__,__NA__) !< result

  real(RNP) :: tmp(0:1)
  integer   :: i, j, k, q

  !$acc loop collapse(3) vector
  do k = 1, __NA__
  do j = 1, __NA__
  do i = 1, __NA__
    v(i,j,k) = 0
  end do
  end do
  end do

  !$acc loop independent vector
  do k = 1, __NA__
    !$acc loop seq
    do q = 1, __NA_T2__, 2
      !$acc loop collapse(2) independent vector private(tmp)
      do j = 1, __NA__
      do i = 1, __NA_T2__, 2

        tmp =       A(q  ,i:i+1) * u(q  ,j,k)
        tmp = tmp + A(q+1,i:i+1) * u(q+1,j,k)

        v(i:i+1,j,k) = v(i:i+1,j,k) + g * tmp

      end do
      end do
    end do
  end do

#if __NA__ == __NA_T2__ + 1

  ! remainder: i = 1:np-1, p = np ..........................................

  !$acc loop collapse(3) vector
  do k = 1, __NA__
  do j = 1, __NA__
  do i = 1, __NA_T2__, 2

    v(i:i+1,j,k) = v(i:i+1,j,k) + g * A(__NA__,i:i+1) * u(__NA__,j,k)

  end do
  end do
  end do

  ! remainder: i = np, p = 1:np-1 ..........................................

  !$acc loop independent vector
  do k = 1, __NA__
    !$acc loop seq
    do q = 1, __NA_T2__, 2
      !$acc loop private(tmp)
      do j = 1, __NA__

        tmp(0) =          A(q  ,__NA__) * u(q  ,j,k)
        tmp(0) = tmp(0) + A(q+1,__NA__) * u(q+1,j,k)

        v(__NA__,j,k) = v(__NA__,j,k) + g *tmp(0)

      end do
    end do
  end do

  ! remainder: i = np, p = np ..............................................

  !$acc loop collapse(2) independent vector
  do k = 1, __NA__
  do j = 1, __NA__
    v(__NA__,j,k) = v(__NA__,j,k) + g * A(__NA__,__NA__) * u(__NA__,j,k)
  end do
  end do

#endif

  ! v = D v --------------------------------------------------------------------

  !$acc loop collapse(3) vector
  do k = 1, __NA__
  do j = 1, __NA__
  do i = 1, __NA__
    v(i,j,k) = D(i,j,k) * v(i,j,k)
  end do
  end do
  end do

end subroutine SubOp_1id
