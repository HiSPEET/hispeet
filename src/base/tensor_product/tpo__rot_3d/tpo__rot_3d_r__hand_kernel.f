!-------------------------------------------------------------------------------
!> Parametrized 3d rotation kernel using hand-crafted suboperators (RCLI)

subroutine PROC(TPO_Rot_R_Hand__,_NP_)(ne, Ds, dx, u, v)
  
  use Constants, only: ZERO, ONE
  
  integer,   intent(in)  :: ne                     !< num elements
  real(RNP), intent(in)  :: Ds(_NP_,_NP_)          !< standard diff matrix
  real(RNP), intent(in)  :: dx(3)                  !< element extensions
  real(RNP), intent(in)  :: u(_NP_,_NP_,_NP_,ne,3) !< 3D vector field
  real(RNP), intent(out) :: v(_NP_,_NP_,_NP_,ne,3) !< element-wise rotation of u

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
    
    ! v2 = du1/dx3
    call PROC(QtxIxI__,_NP_)(A, g(3), ZERO, u(:,:,:,e,1), v(:,:,:,e,2))
    ! v3 =-du1/dx2
    call PROC(IxQtxI__,_NP_)(A,-g(2), ZERO, u(:,:,:,e,1), v(:,:,:,e,3))
    ! v3+= du2/dx1
    call PROC(IxIxQt__,_NP_)(A, g(1), ONE , u(:,:,:,e,2), v(:,:,:,e,3))
    ! v1 =-du2/dx3
    call PROC(QtxIxI__,_NP_)(A,-g(3), ZERO, u(:,:,:,e,2), v(:,:,:,e,1))
    ! v1+= du3/dx2
    call PROC(IxQtxI__,_NP_)(A, g(2), ONE , u(:,:,:,e,3), v(:,:,:,e,1))
    ! v2+=-du3/dx1
    call PROC(IxIxQt__,_NP_)(A,-g(1), ONE , u(:,:,:,e,3), v(:,:,:,e,2))

  end do

  !$acc end parallel
  !$acc end data

  !-----------------------------------------------------------------------------

end subroutine PROC(TPO_Rot_R_Hand__,_NP_)
