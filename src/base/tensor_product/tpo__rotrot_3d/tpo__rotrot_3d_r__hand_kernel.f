!-------------------------------------------------------------------------------
!> Parametrized 3d rotrotrotrot kernel using hand-crafted suboperators (RCLI)

subroutine PROC(TPO_RotRot_R_Hand__,_NP_)(ne, Ds, dx, u, v)
 
  use Constants, only: ZERO, ONE

  integer,   intent(in)  :: ne                     !< num elements
  real(RNP), intent(in)  :: Ds(_NP_,_NP_)          !< standard diff matrix
  real(RNP), intent(in)  :: u(_NP_,_NP_,_NP_,ne,3) !< operand
  real(RNP), intent(out) :: v(_NP_,_NP_,_NP_,ne,3) !< result
  real(RNP), intent(in)  :: dx(3)                  !< element extensions

  real(RNP) :: A(_NP_,_NP_), DA(_NP_,_NP_)
  real(RNP) :: div_u(_NP_,_NP_,_NP_)
  real(RNP) :: g(3)

  integer :: e, i, j

  !-----------------------------------------------------------------------------
  ! initialization

  A  = transpose(Ds)

  do j = 1, _NP_
  do i = 1, _NP_
    DA(j,i) = sum(Ds(i,:) * Ds(:,j))
  end do
  end do

  ! metric coefficients
  g = 2 / dx

  !---------------------------------------------------------------------------
  ! evaluation

  !$acc data present(u,v) copyin(c,M,Lm)
  !$acc parallel
  !$acc loop gang worker private(M_u)

  !$omp do
  do e = 1, ne

    call PROC(IxIxQt__,_NP_)( A, g(1)     , ZERO, u(:,:,:,e,1), div_u(:,:,:))
    call PROC(IxQtxI__,_NP_)( A, g(2)     , ONE , u(:,:,:,e,2), div_u(:,:,:))
    call PROC(QtxIxI__,_NP_)( A, g(3)     , ONE , u(:,:,:,e,3), div_u(:,:,:))
    call PROC(IxIxQt__,_NP_)( A, g(1)     , ZERO, div_u(:,:,:), v(:,:,:,e,1))
    call PROC(IxIxQt__,_NP_)(DA,-g(1)*g(1), ONE , u(:,:,:,e,1), v(:,:,:,e,1))
    call PROC(IxQtxI__,_NP_)(DA,-g(2)*g(2), ONE , u(:,:,:,e,1), v(:,:,:,e,1))
    call PROC(QtxIxI__,_NP_)(DA,-g(3)*g(3), ONE , u(:,:,:,e,1), v(:,:,:,e,1))
    call PROC(IxQtxI__,_NP_)( A, g(2)     , ZERO, div_u(:,:,:), v(:,:,:,e,2))
    call PROC(IxIxQt__,_NP_)(DA,-g(1)*g(1), ONE , u(:,:,:,e,2), v(:,:,:,e,2))
    call PROC(IxQtxI__,_NP_)(DA,-g(2)*g(2), ONE , u(:,:,:,e,2), v(:,:,:,e,2))
    call PROC(QtxIxI__,_NP_)(DA,-g(3)*g(3), ONE , u(:,:,:,e,2), v(:,:,:,e,2)) 
    call PROC(QtxIxI__,_NP_)( A, g(3)     , ZERO, div_u(:,:,:), v(:,:,:,e,3))
    call PROC(IxIxQt__,_NP_)(DA,-g(1)*g(1), ONE , u(:,:,:,e,3), v(:,:,:,e,3))
    call PROC(IxQtxI__,_NP_)(DA,-g(2)*g(2), ONE , u(:,:,:,e,3), v(:,:,:,e,3))
    call PROC(QtxIxI__,_NP_)(DA,-g(3)*g(3), ONE , u(:,:,:,e,3), v(:,:,:,e,3))

  end do

  !$acc end parallel
  !$acc end data

  !-----------------------------------------------------------------------------

end subroutine PROC(TPO_RotRot_R_Hand__,_NP_)
