!-------------------------------------------------------------------------------
!> Parametrized 3d isotropic spectral kernel using hand-crafted suboperators

subroutine PROC(TPO_Spectral_CI_Hand__,_NP_)(ne, S, Lambda, u, v)
  integer,   intent(in)  :: ne                     !< num elements
  real(RNP), intent(in)  :: S(_NP_,_NP_)           !< 1D eigenvectors
  real(RNP), intent(in)  :: Lambda(_NP_,_NP_,_NP_) !< 3D eigenvalues
  real(RNP), intent(in)  :: u(_NP_,_NP_,_NP_,ne)   !< operand
  real(RNP), intent(out) :: v(_NP_,_NP_,_NP_,ne)   !< result

  real(RNP), parameter :: alpha = 1
  real(RNP), parameter :: beta  = 0
  real(RNP) :: St(_NP_,_NP_)
  real(RNP) :: z (_NP_,_NP_,_NP_)
  real(RNP) :: tmp
  integer   :: e, i, j, k, p

  !---------------------------------------------------------------------------
  ! initialization

  St = transpose(S)

  !---------------------------------------------------------------------------
  ! evaluation

  !$acc data present(u,v) copyin(S,St,Lambda)
  !$acc parallel
  !$acc loop gang worker private(z)

  !$omp do private(e)
  do e = 1, ne
    call PROC(QtxIxI__,_NP_)(S, alpha, beta, u(:,:,:,e), z)  ! z  = SᵀxIxI uᵉ
    call PROC(IxQtxI__,_NP_)(S, alpha, beta, z, v(:,:,:,e))  ! vᵉ = IxSᵀxI z
    call PROC(IxIxQt__,_NP_)(S, alpha, beta, v(:,:,:,e), z)  ! z  = IxIxSᵀ vᵉ
    z = Lambda * z
    call PROC(IxIxQt__,_NP_)(St, alpha, beta, z, v(:,:,:,e)) ! vᵉ = IxIxS z
    call PROC(IxQtxI__,_NP_)(St, alpha, beta, v(:,:,:,e), z) ! z  = IxSxI vᵉ
    call PROC(QtxIxI__,_NP_)(St, alpha, beta, z, v(:,:,:,e)) ! vᵉ = SxIxI z
  end do
  !$omp end do

  !$acc end parallel
  !$acc end data

end subroutine PROC(TPO_Spectral_CI_Hand__,_NP_)
