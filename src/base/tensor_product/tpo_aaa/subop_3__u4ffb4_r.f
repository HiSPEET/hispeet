!-------------------------------------------------------------------------------
!> Performs `v = AxIxI z2` for single element
!>
!>   * k: unrolled and jammed with length of 4 (u4)
!>   * j: joined with i and flattened (f)
!>   * i: joined with j and flattened (f)
!>   * p: blocked with length of 4 (b4)
!>   * explicit remainder handling (r)

subroutine SubOp_3(At, z2, v)
  !$acc routine vector
  real(RNP), intent(in)  :: At(__NA2__,__NA1__)    !< transpose of A
  real(RNP), intent(in)  :: z2(__NA1__**2,__NA2__) !< z2 = IxIxA u
  real(RNP), intent(out) :: v (__NA1__**2,__NA1__) !< result

  integer   :: ij, k, p, q

  !$acc loop collapse(2) vector
  do k = 1, __NA1__
  do ij = 1, __NA1__**2
   v(ij,k) = 0
  end do
  end do

  !$acc loop independent vector
  do k = 1, __NA1_T4__, 4
    !$acc loop seq
    do q = 1, __NA2_T4__, 4
      !$acc loop independent vector
      !DIR$ SIMD
      do ij = 1, __NA1__**2
        do p = q, q+3
          v(ij,k  ) = v(ij,k  ) + At(p,k  ) * z2(ij,p)
          v(ij,k+1) = v(ij,k+1) + At(p,k+1) * z2(ij,p)
          v(ij,k+2) = v(ij,k+2) + At(p,k+2) * z2(ij,p)
          v(ij,k+3) = v(ij,k+3) + At(p,k+3) * z2(ij,p)
        end do
      end do
    end do
  end do

#if __NA2__ == __NA2_T4__ + 1

  !$acc loop collapse(2) independent vector
  do k = 1, __NA1_T4__, 4
  !DIR$ SIMD
  do ij = 1, __NA1__**2
    v(ij,k  ) = v(ij,k  ) + At(__NA2__,k  ) * z2(ij,__NA2__)
    v(ij,k+1) = v(ij,k+1) + At(__NA2__,k+1) * z2(ij,__NA2__)
    v(ij,k+2) = v(ij,k+2) + At(__NA2__,k+2) * z2(ij,__NA2__)
    v(ij,k+3) = v(ij,k+3) + At(__NA2__,k+3) * z2(ij,__NA2__)
  end do
  end do

#elif __NA2__ == __NA2_T4__ + 2

  !$acc loop collapse(2) independent vector
  do k = 1, __NA1_T4__, 4
  !DIR$ SIMD
  do ij = 1, __NA1__**2
    do p = __NA2__-1, __NA2__
      v(ij,k  ) = v(ij,k  ) + At(p,k  ) * z2(ij,p)
      v(ij,k+1) = v(ij,k+1) + At(p,k+1) * z2(ij,p)
      v(ij,k+2) = v(ij,k+2) + At(p,k+2) * z2(ij,p)
      v(ij,k+3) = v(ij,k+3) + At(p,k+3) * z2(ij,p)
    end do
  end do
  end do

#elif __NA2__ == __NA2_T4__ + 3

  !$acc loop collapse(2) independent vector
  do k = 1, __NA1_T4__, 4
  !DIR$ SIMD
  do ij = 1, __NA1__**2
    do p = __NA2__-2, __NA2__
      v(ij,k  ) = v(ij,k  ) + At(p,k  ) * z2(ij,p)
      v(ij,k+1) = v(ij,k+1) + At(p,k+1) * z2(ij,p)
      v(ij,k+2) = v(ij,k+2) + At(p,k+2) * z2(ij,p)
      v(ij,k+3) = v(ij,k+3) + At(p,k+3) * z2(ij,p)
    end do
  end do
  end do

#endif

#if __NA1__ == __NA1_T4__ + 1

  !$acc loop vector
  !DIR$ SIMD
  do ij = 1, __NA1__**2
    do p = 1, __NA2__
      v(ij,__NA1__) = v(ij,__NA1__) + At(p,__NA1__) * z2(ij,p)
    end do
  end do

#elif __NA1__ == __NA1_T4__ + 2

  !$acc loop vector
  !DIR$ SIMD
  do ij = 1, __NA1__**2
    do p = 1, __NA2__
      v(ij,__NA1__-1) = v(ij,__NA1__-1) + At(p,__NA1__-1) * z2(ij,p)
      v(ij,__NA1__  ) = v(ij,__NA1__  ) + At(p,__NA1__  ) * z2(ij,p)
    end do
  end do

#elif __NA1__ == __NA1_T4__ + 3

  !$acc loop vector
  !DIR$ SIMD
  do ij = 1, __NA1__**2
    do p = 1, __NA2__
      v(ij,__NA1__-2) = v(ij,__NA1__-2) + At(p,__NA1__-2) * z2(ij,p)
      v(ij,__NA1__-1) = v(ij,__NA1__-1) + At(p,__NA1__-1) * z2(ij,p)
      v(ij,__NA1__  ) = v(ij,__NA1__  ) + At(p,__NA1__  ) * z2(ij,p)
    end do
  end do

#endif

end subroutine SubOp_3
