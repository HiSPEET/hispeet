!-------------------------------------------------------------------------------
!> Parametrized 3d divergence kernel using hand-crafted suboperators (RCLI)

subroutine PROC(TPO_Div_R_Hand__,_NP_)(ne, Ds, dx, u, v)

  use Constants, only: ZERO, ONE

  integer,   intent(in)  :: ne                    !< num elements
  real(RNP), intent(in)  :: Ds(_NP_,_NP_)         !< standard diff matrix
  real(RNP), intent(in)  :: u(_NP_,_NP_,_NP_,ne,3)!< 3D vector field
  real(RNP), intent(out) :: v(_NP_,_NP_,_NP_,ne) !< element-wise divergence of u
  real(RNP), intent(in)  :: dx(3)                 !< element extensions

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
    
    v(:,:,:,e) = 0  ! avoid trouble with NaNs

    ! v is initialized   v = du1/dx1 (beta = 0)
    call PROC(IxIxQt__,_NP_)(A, g(1), ZERO , u(:,:,:,e,1), v(:,:,:,e))
    ! v is added up with v+= du2/dx2 (beta = 1) 
    call PROC(IxQtxI__,_NP_)(A, g(2), ONE  , u(:,:,:,e,2), v(:,:,:,e))
    ! v is added up with v+= du3/dx3 (beta = 1)
    call PROC(QtxIxI__,_NP_)(A, g(3), ONE  , u(:,:,:,e,3), v(:,:,:,e))

  end do

  !$acc end parallel
  !$acc end data

  !-----------------------------------------------------------------------------

end subroutine PROC(TPO_Div_R_Hand__,_NP_)
