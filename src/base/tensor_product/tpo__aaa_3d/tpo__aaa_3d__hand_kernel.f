!-------------------------------------------------------------------------------
!> Parametrized 3d AxAxA kernel using hand-crafted suboperators

subroutine PROC(TPO_AxAxA_Hand__,_NA1_,_NA2_)(ne, u, v)
  integer,   intent(in)  :: ne                       !< num elements
  real(RNP), intent(in)  :: A(_NA1_,_NA2_)           !< 1D operator
  real(RNP), intent(in)  :: u(_NA2_,_NA2_,_NA2_,ne)  !< operand
  real(RNP), intent(out) :: v(_NA1_,_NA1_,_NA1_,ne)  !< result

  real(RNP), parameter :: alpha = 0
  real(RNP), parameter :: beta  = 1
  real(RNP) :: At(_NA2_,_NA2_)
  real(RNP) :: z2(_NA2_,_NA1_,_NA1_)
  real(RNP) :: z3(_NA2_,_NA2_,_NA1_)

  integer, parameter :: NA2_NA2 = _NA2_ * _NA2_
  integer :: e, i, j, k

  !-----------------------------------------------------------------------------
  ! initialization

  At = transpose(A)

  !-----------------------------------------------------------------------------
  ! evaluation

  !$acc data present(u,v) copyin(At)
  !$acc parallel
  !$acc loop gang worker private(z2,z3)

  !$omp do
  do e = 1, ne
    call PROC(CtxIxI__,_NA1_,_NA2_)(NA2_NA2, At, alpha, beta, u(:,:,:,e), z3)
    call PROC(IxBtxI__,_NA1_,_NA2_)(_NA2_, _NA1_, At, alpha, beta, z3, z2)
    call PROC(IxIxAt__,_NA1_,_NA2_)(_NA1_, _NA1_, At, alpha, beta, z2, v(:,:,:,e))
  end do

  !$acc end parallel
  !$acc end data

  !-----------------------------------------------------------------------------

end subroutine PROC(TPO_Elliptic_RLCI_Hand__,_NP_)
