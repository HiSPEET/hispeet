!-------------------------------------------------------------------------------
!> Computes  v = α Cᵀ⊗I⊗I u + β v
!>
!>   * k:  block 1,  unroll 4
!>   * j:  joined with i
!>   * i:  joined with j
!>   * p:  block 4,  unroll 1
!>   * using Intel SIMD directive

#define _NC2_T4_ (_NC2_ / 4) * 4
#define _NC1_T4_ (_NC1_ / 4) * 4

subroutine PROC(CtxIxI__,_NC1_,_NC2_)(na_nb, Ct, alpha, beta, u, v)
  !$acc routine vector
  integer,   intent(in)    :: na_nb           !< size of dimensions 1+2 of u,v
  real(RNP), intent(in)    :: Ct(_NC2_,_NC1_)  !< transpose of C
  real(RNP), intent(in)    :: alpha           !< factor α
  real(RNP), intent(in)    :: beta            !< factor β
  real(RNP), intent(in)    :: u(na_nb,_NC2_)  !< operand
  real(RNP), intent(inout) :: v(na_nb,_NC1_)  !< result

  integer   :: ij, k, p, pb

  do k = 1, _NC2_
  do ij = 1, na_nb
   v(ij,k) = beta * v(ij,k)
  end do
  end do

  !$acc loop independent vector
  do k = 1, _NC1_T4_, 4
    !$acc loop seq
    do pb = 1, _NC2_T4_, 4
      !$acc loop independent vector
      do ij = 1, na_nb
        do p = pb, pb+3
          v(ij,k  ) = alpha * Ct(p,k  ) * u(ij,p) + beta * v(ij,k  )
          v(ij,k+1) = alpha * Ct(p,k+1) * u(ij,p) + beta * v(ij,k+1)
          v(ij,k+2) = alpha * Ct(p,k+2) * u(ij,p) + beta * v(ij,k+2)
          v(ij,k+3) = alpha * Ct(p,k+3) * u(ij,p) + beta * v(ij,k+3)
        end do
      end do
    end do
  end do

#if _NC2_ == _NC2_T4_ + 1

  !$acc loop collapse(2) independent vector
  do k = 1, _NC1_T4_, 4
  do ij = 1, na_nb
    v(ij,k  ) = alpha * Ct(_NC2_,k  ) * u(ij,_NC2_) + beta * v(ij,k  )
    v(ij,k+1) = alpha * Ct(_NC2_,k+1) * u(ij,_NC2_) + beta * v(ij,k+1)
    v(ij,k+2) = alpha * Ct(_NC2_,k+2) * u(ij,_NC2_) + beta * v(ij,k+2)
    v(ij,k+3) = alpha * Ct(_NC2_,k+3) * u(ij,_NC2_) + beta * v(ij,k+3)
  end do
  end do

#elif _NC2_ == _NC2_T4_ + 2

  !$acc loop collapse(2) independent vector
  do k = 1, _NC1_T4_, 4
  !DIR$ SIMD
  do ij = 1, na_nb
    do p = _NC2_-1, _NC2_
      v(ij,k  ) = alpha * Ct(p,k  ) * u(ij,p) + beta * v(ij,k  )
      v(ij,k+1) = alpha * Ct(p,k+1) * u(ij,p) + beta * v(ij,k+1)
      v(ij,k+2) = alpha * Ct(p,k+2) * u(ij,p) + beta * v(ij,k+2)
      v(ij,k+3) = alpha * Ct(p,k+3) * u(ij,p) + beta * v(ij,k+3)
    end do
  end do
  end do

#elif _NC2_ == _NC2_T4_ + 3

  !$acc loop collapse(2) independent vector
  do k = 1, _NC1_T4_, 4
  !DIR$ SIMD
  do ij = 1, na_nb
    do p = _NC2_-2, _NC2_
      v(ij,k  ) = alpha * Ct(p,k  ) * u(ij,p) + beta * v(ij,k  )
      v(ij,k+1) = alpha * Ct(p,k+1) * u(ij,p) + beta * v(ij,k+1)
      v(ij,k+2) = alpha * Ct(p,k+2) * u(ij,p) + beta * v(ij,k+2)
      v(ij,k+3) = alpha * Ct(p,k+3) * u(ij,p) + beta * v(ij,k+3)
    end do
  end do
  end do

#endif

#if _NC1_ == _NC1_T4_ + 1

  !$acc loop vector
  !DIR$ SIMD
  do ij = 1, na_nb
    do p = 1, _NC2_
      v(ij,_NC1_) = alpha * Ct(p,_NC1_) * u(ij,p) + beta * v(ij,_NC1_)
    end do
  end do

#elif _NC1_ == _NC1_T4_ + 2

  !$acc loop vector
  !DIR$ SIMD
  do ij = 1, na_nb
    do p = 1, _NC2_
      v(ij,_NC1_-1) = alpha * Ct(p,_NC1_-1) * u(ij,p) + beta * v(ij,_NC1_-1)
      v(ij,_NC1_  ) = alpha * Ct(p,_NC1_  ) * u(ij,p) + beta * v(ij,_NC1_  )
    end do
  end do

#elif _NC1_ == _NC1_T4_ + 3

  !$acc loop vector
  !DIR$ SIMD
  do ij = 1, na_nb
    do p = 1, _NC2_
      v(ij,_NC1_-2) = alpha * Ct(p,_NC1_-2) * u(ij,p) + beta * v(ij,_NC1_-2)
      v(ij,_NC1_-1) = alpha * Ct(p,_NC1_-1) * u(ij,p) + beta * v(ij,_NC1_-1)
      v(ij,_NC1_  ) = alpha * Ct(p,_NC1_  ) * u(ij,p) + beta * v(ij,_NC1_  )
    end do
  end do

#endif

end subroutine PROC(CtxIxI__,_NC1_,_NC2_)

#undef _NC2_T4_
#undef _NC1_T4_
