!-------------------------------------------------------------------------------
!> Parametrized 3d gradient kernel using hand-crafted suboperators (RCLI)

subroutine PROC(TPO_Grad_R_Hand__,_NP_)(ne, Ds, dx, u, v)

  use Constants, only: ZERO

  integer,   intent(in)  :: ne                     !< num elements
  real(RNP), intent(in)  :: Ds(_NP_,_NP_)          !< standard diff matrix
  real(RNP), intent(in)  :: u(_NP_,_NP_,_NP_,ne)   !< operand
  real(RNP), intent(out) :: v(_NP_,_NP_,_NP_,ne,3) !< result
  real(RNP), intent(in)  :: dx(3)                  !< element extensions

  real(RNP) :: A(_NP_,_NP_)
  real(RNP) :: g(3)

  integer :: e

  !-----------------------------------------------------------------------------
  ! initialization

  A = transpose(Ds)

  ! metric coefficients
  g = 2 / dx

  !---------------------------------------------------------------------------
  ! evaluation

  !$acc data present(u,v) copyin(c,M,Lm)
  !$acc parallel
  !$acc loop gang worker private(M_u)

  !$omp do
  do e = 1, ne

    call PROC(IxIxQt__,_NP_)(A, g(1), ZERO, u(:,:,:,e), v(:,:,:,e,1))
    call PROC(IxQtxI__,_NP_)(A, g(2), ZERO, u(:,:,:,e), v(:,:,:,e,2))
    call PROC(QtxIxI__,_NP_)(A, g(3), ZERO, u(:,:,:,e), v(:,:,:,e,3))

  end do

  !$acc end parallel
  !$acc end data

  !-----------------------------------------------------------------------------

end subroutine PROC(TPO_Grad_R_Hand__,_NP_)
