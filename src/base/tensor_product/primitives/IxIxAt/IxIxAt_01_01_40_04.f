!-------------------------------------------------------------------------------
!> Computes  v = α I⊗I⊗Aᵀ u + β v
!>
!>   * k:  block 1,  unroll 1
!>   * j:  block 1,  unroll 1
!>   * i:  block 4,  unroll 1
!>   * p:  block 1,  unroll 4
!>   * using Intel SIMD directive

subroutine PROC(IxIxAt__,_NA1_,_NA2_)(nb, nc, A, alpha, beta, u, v)
  !$acc routine vector
  integer,   intent(in)    :: nb              !< 2nd dimension of u,v
  integer,   intent(in)    :: nc              !< 3rd dimension of u,v
  real(RNP), intent(in)    :: A(_NA1_,_NA2_)  !< rectangular matrix
  real(RNP), intent(in)    :: alpha           !< factor α
  real(RNP), intent(in)    :: beta            !< factor β
  real(RNP), intent(in)    :: u(_NA1_,nb,nc)  !< operand
  real(RNP), intent(inout) :: v(_NA2_,nb,nc)  !< result

  real(RNP) :: tmp(0:3)
  integer   :: i, j, k, p

  !$acc loop independent vector
  do k = 1, nc
    !$acc loop seq
    do p = 1, _NA1_T4_, 4
      !$acc loop collapse(2) independent vector private(tmp)
      do j = 1, nb
      do i = 1, _NA2_T4_, 4

        tmp =       A(p  ,i:i+3) * u(p  ,j,k)
        tmp = tmp + A(p+1,i:i+3) * u(p+1,j,k)
        tmp = tmp + A(p+2,i:i+3) * u(p+2,j,k)
        tmp = tmp + A(p+3,i:i+3) * u(p+3,j,k)

        v(i:i+3,j,k) = alpha * tmp + beta * v(i:i+3,j,k)

      end do
      end do
    end do
  end do

#if _NA1_ == _NA1_T4__ + 1

  ! 1:na1/4*4, na2 -------------------------------------------------------------

  p = _NA1_T4_ + 1

  !$acc loop collapse(3) independent vector
  do k = 1, nc
    do j = 1, nb
    do i = 1, _NA2_T4_, 4

      v(i:i+3,j,k) = v(i:i+3,j,k) + A(p,i:i+3) * u(p,j,k)

    end do
    end do
  end do

#elif _NA1_ == _NA1_T4__ + 2

  ! 1:na1/4*4, na2-1:na2 -------------------------------------------------------

  p = _NA1_T4_ + 1

  !$acc loop collapse(3) independent vector private(tmp)
  do k = 1, nc
    do j = 1, nb
    do i = 1, _NA2_T4_, 4

      tmp =       A(p  ,i:i+3) * u(p  ,j,k)
      tmp = tmp + A(p+1,i:i+3) * u(p+1,j,k)

      v(i:i+3,j,k) = alpha * tmp + beta * v(i:i+3,j,k)

    end do
    end do
  end do

#elif _NA1_ == _NA1_T4__ + 3

  ! 1:na1/4*4, na2-2:na2 -------------------------------------------------------

  p = _NA1_T4_ + 1

  !$acc loop collapse(3) independent vector private(tmp)
  do k = 1, nc
    do j = 1, nb
    do i = 1, _NA2_T4_, 4

      tmp =       A(p  ,i:i+3) * u(p  ,j,k)
      tmp = tmp + A(p+1,i:i+3) * u(p+1,j,k)
      tmp = tmp + A(p+2,i:i+3) * u(p+2,j,k)

      v(i:i+3,j,k) = alpha * tmp + beta * v(i:i+3,j,k)

    end do
    end do
  end do

#endif

#if _NA2_ == _NA2_T4_ + 1

  ! na1, 1:na2 -----------------------------------------------------------------

  i = _NA2_T4_ + 1

  !$acc loop collapse(2) independent vector
  do k = 1, nc
  do j = 1, nb
    do p = 1, _NA1_

      v(i,j,k) = v(i,j,k) + A(p,i) * u(p,j,k)

    end do
  end do
  end do

#elif _NA2_ == _NA2_T4_ + 2

  ! na1-1:na1, 1:na2 -----------------------------------------------------------

  i = _NA2_T4_ + 1

  !$acc loop collapse(2) independent vector
  do k = 1, nc
  do j = 1, nb
    do p = 1, _NA1_

      v(i:i+1,j,k) = v(i:i+1,j,k) + A(p,i:i+1) * u(p,j,k)

    end do
  end do
  end do

#elif _NA2_ == _NA2_T4_ + 3

  ! na1-2:na1, 1:na2 -----------------------------------------------------------

  i = _NA1_T4_ + 1

  !$acc loop collapse(2) independent vector
  do k = 1, nc
  do j = 1, nb
    do p = 1, _NA1_

      v(i:i+2,j,k) = v(i:i+2,j,k) + A(p,i:i+2) * u(p,j,k)

    end do
  end do
  end do

#endif

end subroutine PROC(IxIxAt__,_NA1_,_NA2_)
