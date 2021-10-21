!-------------------------------------------------------------------------------
!> Parametrized 3d diffusion kernel using hand-crafted suboperators (RLCI)

subroutine PROC(TPO_Diffusion_DLCI_Hand__,_NP_) &
  (ne, Ms, Ds, Jd, G, lambda, nu, u, v)
  integer,   intent(in)  :: ne                     !< num elements
  real(RWP), intent(in)  :: Ms(_NP_)               !< standard mass matrix
  real(RWP), intent(in)  :: Ds(_NP_,_NP_)          !< 1D standard diff matrix
  real(RWP), intent(in)  :: Jd(_NP_,_NP_,_NP_,ne)  !< Jacobian determinant
  real(RWP), intent(in)  :: G(_NP_,_NP_,_NP_,ne,6) !< Laplacian metrics
  real(RWP), intent(in)  :: lambda                 !< Helmholtz parameter λ
  real(RWP), intent(in)  :: nu                     !< diffusivity nu
  real(RWP), intent(in)  :: u(_NP_,_NP_,_NP_,ne)   !< operand
  real(RWP), intent(out) :: v(_NP_,_NP_,_NP_,ne)   !< result

  real(RWP), parameter :: ONE = 1
  real(RWP), parameter :: ZERO = 0

  real(RWP) :: M(_NP_,_NP_,_NP_), z(_NP_,_NP_,_NP_)
  real(RWP) :: r(_NP_,_NP_,_NP_), s(_NP_,_NP_,_NP_), t(_NP_,_NP_,_NP_)

  real(RWP) :: Ds_t(_NP_,_NP_)

  integer, parameter :: NP_E = _NP_**3
  integer :: e, i, j, k

  !-----------------------------------------------------------------------------
  ! initialization

  r = 0   ! avoid trouble with NaNs
  s = 0   ! avoid trouble with NaNs
  t = 0   ! avoid trouble with NaNs

  Ds_t = transpose(Ds)

    ! element mass matrix
    do k = 1, _NP_
    do j = 1, _NP_
    do i = 1, _NP_
      M(i,j,k) =  Ms(k) * Ms(j) * Ms(i)
    end do
    end do
    end do

  !---------------------------------------------------------------------------
  ! evaluation

  !$omp do
  do e = 1, ne

    ! projection to element mass matrix and computation of the Helmholtz term
    do k = 1, _NP_
    do j = 1, _NP_
    do i = 1, _NP_
      v(i,j,k,e) = lambda * M(i,j,k) * Jd(i,j,k,e) * u(i,j,k,e)
    end do
    end do
    end do

    ! first derivatives

    call PROC(IxIxQt__,_NP_)(Ds_t, ONE, ZERO, u(:,:,:,e), r)
    call PROC(IxQtxI__,_NP_)(Ds_t, ONE, ZERO, u(:,:,:,e), s)
    call PROC(QtxIxI__,_NP_)(Ds_t, ONE, ZERO, u(:,:,:,e), t)


    ! second derivative direction 1 and application of metric terms

    do k = 1, _NP_
    do j = 1, _NP_
    do i = 1, _NP_
      z(i,j,k) = nu * M(i,j,k) * ( G(i,j,k,e,1) * r(i,j,k) &
                                 + G(i,j,k,e,2) * s(i,j,k) &
                                 + G(i,j,k,e,3) * t(i,j,k) )
    end do
    end do
    end do

    call PROC(IxIxQt__,_NP_)(Ds, ONE, ONE, z(:,:,:), v(:,:,:,e))


    ! second derivative direction 2 and application of metric terms

    do k = 1, _NP_
    do j = 1, _NP_
    do i = 1, _NP_
      z(i,j,k) = nu * M(i,j,k) * ( G(i,j,k,e,2) * r(i,j,k) &
                                 + G(i,j,k,e,4) * s(i,j,k) &
                                 + G(i,j,k,e,5) * t(i,j,k) )
    end do
    end do
    end do

    call PROC(IxQtxI__,_NP_)(Ds, ONE, ONE, z(:,:,:), v(:,:,:,e))


    ! second derivative direction 3 and application of metric terms

    do k = 1, _NP_
    do j = 1, _NP_
    do i = 1, _NP_
      z(i,j,k) = nu * M(i,j,k) * ( G(i,j,k,e,3) * r(i,j,k) &
                                 + G(i,j,k,e,5) * s(i,j,k) &
                                 + G(i,j,k,e,6) * t(i,j,k) )
    end do
    end do
    end do

    call PROC(QtxIxI__,_NP_)(Ds, ONE, ONE, z(:,:,:), v(:,:,:,e))


  end do


  !-----------------------------------------------------------------------------

end subroutine PROC(TPO_Diffusion_DLCI_Hand__,_NP_)
