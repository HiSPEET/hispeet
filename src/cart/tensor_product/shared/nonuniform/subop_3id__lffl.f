!-------------------------------------------------------------------------------
!> Performs `v = g D CᵀxIxI u` for single element
!>
!>   * k: simple loop (l)
!>   * j: joined with i and flattened (f)
!>   * i: joined with j and flattened (f)
!>   * p: simple loop (l)

subroutine SubOp_3id(g, C, D, u, v)
  !$acc routine vector
  real(RNP), intent(in)  :: g                       !< metric factor
  real(RNP), intent(in)  :: C(__NC__,__NC__)        !< C
  real(RNP), intent(in)  :: D(__NA__*__NB__,__NC__) !< diagonal operator
  real(RNP), intent(in)  :: u(__NA__*__NB__,__NC__) !< operand
  real(RNP), intent(out) :: v(__NA__*__NB__,__NC__) !< result

  real(RNP) :: tmp
  integer   :: ij, k, p

  !$acc loop collapse(2) vector
  do k = 1, __NC__
  !DIR$ SIMD VECREMAINDER
  do ij = 1, __NA__*__NB__
    tmp = 0
    do p = 1, __NC__
      tmp = tmp + C(p,k) * u(ij,p)
    end do
    v(ij,k) = g * D(ij,k) * tmp
  end do
  end do

end subroutine SubOp_3id
