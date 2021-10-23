!-------------------------------------------------------------------------------
!> Parametrized 3d divergence kernel using hand-crafted suboperators (D)

subroutine PROC(TPO_Div_D_Hand__,_NP_)(ne, Ds, Ji, u, v)

  use Constants, only: ZERO, ONE

  integer,   intent(in)  :: ne                        !< num elements
  real(RWP), intent(in)  :: Ds(_NP_,_NP_)             !< standard diff matrix
  real(RWP), intent(in)  :: Ji(_NP_,_NP_,_NP_,ne,3,3) !< inverse Jacobi matrix
  real(RWP), intent(in)  :: u(_NP_,_NP_,_NP_,ne,3)    !< 3D vector field
  real(RWP), intent(out) :: v(_NP_,_NP_,_NP_,ne) !< element-wise divergence of u

  real(RWP) :: r(_NP_,_NP_,_NP_)
  real(RWP) :: s(_NP_,_NP_,_NP_)
  real(RWP) :: t(_NP_,_NP_,_NP_)
  real(RWP) :: w(_NP_,_NP_,_NP_)
  real(RWP) :: A(_NP_,_NP_)

  integer :: e

  !-----------------------------------------------------------------------------
  ! initialization

  A = transpose(Ds)

  r = 0  ! avoid trouble with NaNs
  s = 0  ! avoid trouble with NaNs
  t = 0  ! avoid trouble with NaNs
  w = 0  ! avoid trouble with NaNs

  !-----------------------------------------------------------------------------
  ! evaluation

  !$omp do
  do e = 1, ne

    ! direction 1 ..............................................................

    ! r = ∂u₁/∂ξ
    call PROC(IxIxQt__,_NP_)(A, ONE, ZERO , u(:,:,:,e,1), r)
    ! s = ∂u₁/∂η
    call PROC(IxQtxI__,_NP_)(A, ONE, ZERO , u(:,:,:,e,1), s)
    ! t = ∂u₁/∂ζ
    call PROC(QtxIxI__,_NP_)(A, ONE, ZERO , u(:,:,:,e,1), t)

    w = r * Ji(:,:,:,e,1,1) &
      + s * Ji(:,:,:,e,2,1) &
      + t * Ji(:,:,:,e,3,1)

    ! direction 2 ..............................................................

    ! r = ∂u₂/∂ξ
    call PROC(IxIxQt__,_NP_)(A, ONE, ZERO , u(:,:,:,e,2), r)
    ! s = ∂u₂/∂η
    call PROC(IxQtxI__,_NP_)(A, ONE, ZERO , u(:,:,:,e,2), s)
    ! t = ∂u₂/∂ζ
    call PROC(QtxIxI__,_NP_)(A, ONE, ZERO , u(:,:,:,e,2), t)

    w = w + r * Ji(:,:,:,e,1,2) &
          + s * Ji(:,:,:,e,2,2) &
          + t * Ji(:,:,:,e,3,2)

    ! direction 3 ..............................................................

    ! r = ∂u₃/∂ξ
    call PROC(IxIxQt__,_NP_)(A, ONE, ZERO , u(:,:,:,e,3), r)
    ! s = ∂u₃/∂η
    call PROC(IxQtxI__,_NP_)(A, ONE, ZERO , u(:,:,:,e,3), s)
    ! t = ∂u₃/∂ζ
    call PROC(QtxIxI__,_NP_)(A, ONE, ZERO , u(:,:,:,e,3), t)


    w = w + r * Ji(:,:,:,e,1,3) &
          + s * Ji(:,:,:,e,2,3) &
          + t * Ji(:,:,:,e,3,3)

    v(:,:,:,e) = w

  end do

  !-----------------------------------------------------------------------------

end subroutine PROC(TPO_Div_D_Hand__,_NP_)
