!-------------------------------------------------------------------------------
!> Computes  v = α Cᵀ⊗I⊗I u + β v
!>
!>   * k:  block 1,  unroll 4
!>   * j:  joined with i
!>   * i:  joined with j
!>   * p:  block 4,  unroll 1
!>   * using Intel SIMD directive

subroutine PROC(CtxIxI__,_NC1_,_NC2_)(na_nb, C, alpha, beta, u, v)
  !$acc routine vector
  integer,   intent(in)    :: na_nb           !< size of dimensions 1+2 of u,v
  real(RNP), intent(in)    :: C(_NC1_,_NC2_)  !< transpose of C
  real(RNP), intent(in)    :: alpha           !< factor α
  real(RNP), intent(in)    :: beta            !< factor β
  real(RNP), intent(in)    :: u(na_nb,_NC1_)  !< operand
  real(RNP), intent(inout) :: v(na_nb,_NC2_)  !< result

  integer   :: ij, k, p, pb


  !!! ???
  !$acc loop collapse(2) vector
  do k = 1, _NC2_
  do ij = 1, na_nb
   v(ij,k) = 0
  end do
  end do

  !$acc loop independent vector
  do k = 1, _NC2_T4_, 4
    !$acc loop seq
    do pb = 1, _NC1_T4_, 4
      !$acc loop independent vector
      !DIR$ SIMD
      do ij = 1, na_nb
        do p = pb, pb+3
          v(ij,k  ) = alpha * C(p,k  ) * u(ij,p) + beta * v(ij,k  )
          v(ij,k+1) = alpha * C(p,k+1) * u(ij,p) + beta * v(ij,k+1)
          v(ij,k+2) = alpha * C(p,k+2) * u(ij,p) + beta * v(ij,k+2)
          v(ij,k+3) = alpha * C(p,k+3) * u(ij,p) + beta * v(ij,k+3)
        end do
      end do
    end do
  end do

#if _NC1_ == _NC1_T4_ + 1

  !$acc loop collapse(2) independent vector
  do k = 1, _NC2_T4_, 4
  !DIR$ SIMD
  do ij = 1, na_nb
    v(ij,k  ) = alpha * C(_NC1_,k  ) * u(ij,_NC1_) + beta * v(ij,k  )
    v(ij,k+1) = alpha * C(_NC1_,k+1) * u(ij,_NC1_) + beta * v(ij,k+1)
    v(ij,k+2) = alpha * C(_NC1_,k+2) * u(ij,_NC1_) + beta * v(ij,k+2)
    v(ij,k+3) = alpha * C(_NC1_,k+3) * u(ij,_NC1_) + beta * v(ij,k+3)
  end do
  end do

#elif _NC1_ == _NC1_T4_ + 2

  !$acc loop collapse(2) independent vector
  do k = 1, _NC2_T4_, 4
  !DIR$ SIMD
  do ij = 1, na_nb
    do p = _NC1_-1, _NC1_
      v(ij,k  ) = alpha * C(p,k  ) * u(ij,p) + beta * v(ij,k  )
      v(ij,k+1) = alpha * C(p,k+1) * u(ij,p) + beta * v(ij,k+1)
      v(ij,k+2) = alpha * C(p,k+2) * u(ij,p) + beta * v(ij,k+2)
      v(ij,k+3) = alpha * C(p,k+3) * u(ij,p) + beta * v(ij,k+3)
    end do
  end do
  end do

#elif _NC1_ == _NC1_T4_ + 3

  !$acc loop collapse(2) independent vector
  do k = 1, _NC2_T4_, 4
  !DIR$ SIMD
  do ij = 1, na_nb
    do p = _NC1_-2, _NC1_
      v(ij,k  ) = alpha * C(p,k  ) * u(ij,p) + beta * v(ij,k  )
      v(ij,k+1) = alpha * C(p,k+1) * u(ij,p) + beta * v(ij,k+1)
      v(ij,k+2) = alpha * C(p,k+2) * u(ij,p) + beta * v(ij,k+2)
      v(ij,k+3) = alpha * C(p,k+3) * u(ij,p) + beta * v(ij,k+3)
    end do
  end do
  end do

#endif

#if _NC2_ == _NC2_T4_ + 1

  !$acc loop vector
  !DIR$ SIMD
  do ij = 1, na_nb
    do p = 1, _NC1_
      v(ij,_NC2_) = alpha * C(p,_NC2_) * u(ij,p) + beta * v(ij,_NC2_)
    end do
  end do

#elif _NC2_ == _NC2_T4_ + 2

  !$acc loop vector
  !DIR$ SIMD
  do ij = 1, na_nb
    do p = 1, _NC1_
      v(ij,_NC2_-1) = alpha * C(p,_NC2_-1) * u(ij,p) + beta * v(ij,_NC2_-1)
      v(ij,_NC2_  ) = alpha * C(p,_NC2_  ) * u(ij,p) + beta * v(ij,_NC2_  )
    end do
  end do

#elif _NC2_ == _NC2_T4_ + 3

  !$acc loop vector
  !DIR$ SIMD
  do ij = 1, na_nb
    do p = 1, _NC1_
      v(ij,_NC2_-2) = alpha * C(p,_NC2_-2) * u(ij,p) + beta * v(ij,_NC2_-2)
      v(ij,_NC2_-1) = alpha * C(p,_NC2_-1) * u(ij,p) + beta * v(ij,_NC2_-1)
      v(ij,_NC2_  ) = alpha * C(p,_NC2_  ) * u(ij,p) + beta * v(ij,_NC2_  )
    end do
  end do

#endif

end subroutine PROC(CtxIxI__,_NC1_,_NC2_)
