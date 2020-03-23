!-------------------------------------------------------------------------------
!> Parametrized 3d rotation kernel using hand-crafted suboperators (RCLI)

subroutine PROC(TPO_Rot_R_Hand__,_NP_)(ne, Ds, dx, u, v)
  integer,   intent(in)  :: ne                     !< num elements
  real(RNP), intent(in)  :: Ds(_NP_,_NP_)          !< standard diff matrix
  real(RNP), intent(in)  :: u(_NP_,_NP_,_NP_,ne,3) !< operand
  real(RNP), intent(out) :: v(_NP_,_NP_,_NP_,ne,3) !< result
  real(RNP), intent(in)  :: dx(3)                  !< element extensions

  real(RNP) :: A(_NP_,_NP_)
  real(RNP), parameter :: beta_1 = 0
  real(RNP), parameter :: beta   = 1
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

    call PROC(QtxIxI__,_NP_)(A, g(3), beta_1, u(:,:,:,e,1), v(:,:,:,e,2))
    call PROC(IxQtxI__,_NP_)(A,-g(2), beta_1, u(:,:,:,e,1), v(:,:,:,e,3))
    call PROC(IxIxQt__,_NP_)(A, g(1), beta  , u(:,:,:,e,2), v(:,:,:,e,3))
    call PROC(QtxIxI__,_NP_)(A,-g(3), beta_1, u(:,:,:,e,2), v(:,:,:,e,1))
    call PROC(IxQtxI__,_NP_)(A, g(2), beta  , u(:,:,:,e,3), v(:,:,:,e,1))
    call PROC(IxIxQt__,_NP_)(A,-g(1), beta  , u(:,:,:,e,3), v(:,:,:,e,2))

  end do

  !$acc end parallel
  !$acc end data

  !-----------------------------------------------------------------------------

end subroutine PROC(TPO_Rot_R_Hand__,_NP_)
