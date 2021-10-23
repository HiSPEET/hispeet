!-------------------------------------------------------------------------------
!> Parametrized 3d gradient kernel using hand-crafted suboperators (D)

subroutine PROC(TPO_Grad_D_Hand__,_NP_)(ne, Ds, Ji, u, v)

  use Constants, only: ZERO, ONE

  integer,   intent(in)  :: ne                        !< num elements
  real(RWP), intent(in)  :: Ds(_NP_,_NP_)             !< standard diff matrix
  real(RWP), intent(in)  :: Ji(_NP_,_NP_,_NP_,ne,3,3) !< element extensions
  real(RWP), intent(in)  :: u(_NP_,_NP_,_NP_,ne)   !< 3D scalar field
  real(RWP), intent(out) :: v(_NP_,_NP_,_NP_,ne,3) !< element-wise gradient of u

  real(RWP) :: r(_NP_,_NP_,_NP_)
  real(RWP) :: s(_NP_,_NP_,_NP_)
  real(RWP) :: t(_NP_,_NP_,_NP_)

  real(RWP) :: A(_NP_,_NP_)

  integer :: e


  !-----------------------------------------------------------------------------
  ! initialization

  A = transpose(Ds)

  r = 0   ! avoid trouble with NaNs and memory warnings
  s = 0   ! avoid trouble with NaNs and memory warnings
  t = 0   ! avoid trouble with NaNs and memory warnings

  !---------------------------------------------------------------------------
  ! evaluation

  !$omp do
  do e = 1, ne

    call PROC(IxIxQt__,_NP_)(A, ONE, ZERO, u(:,:,:,e), r)

    call PROC(IxQtxI__,_NP_)(A, ONE, ZERO, u(:,:,:,e), s)

    call PROC(QtxIxI__,_NP_)(A, ONE, ZERO, u(:,:,:,e), t)


    v(:,:,:,e,1) = r * Ji(:,:,:,e,1,1) + &
                   s * Ji(:,:,:,e,2,1) + &
                   t * Ji(:,:,:,e,3,1)

    v(:,:,:,e,2) = r * Ji(:,:,:,e,1,2) + &
                   s * Ji(:,:,:,e,2,2) + &
                   t * Ji(:,:,:,e,3,2)

    v(:,:,:,e,3) = r * Ji(:,:,:,e,1,3) + &
                   s * Ji(:,:,:,e,2,3) + &
                   t * Ji(:,:,:,e,3,3)

  end do


  !-----------------------------------------------------------------------------

end subroutine PROC(TPO_Grad_D_Hand__,_NP_)
