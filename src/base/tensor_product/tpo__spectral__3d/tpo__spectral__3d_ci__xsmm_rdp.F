!> summary:  3d generic isotropic spectral operator using LIBXSMM
!> author:   Erik Pfister, Joerg Stiller
!> date:     2020/06/01
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Spectral__3D_CI__XSMM_RDP

#ifdef __LIBXSMM__

  use Kind_Parameters, only: RDP
  use LIBXSMM
  implicit none
  private

  public :: TPO_Spectral_CI_XSMM

  interface TPO_Spectral_CI_XSMM
    module procedure TPO_Spectral_CI_XSMM_RDP
  end interface

contains

  !-----------------------------------------------------------------------------
  !> Generic procedure for evaluating  v = [(SxSxS) Λ (SxSxS)ᵀ] u

  subroutine TPO_Spectral_CI_XSMM_RDP(S, Lambda, u, v)
    real(RDP), intent(in)  :: S(:,:)        !< 1D eigenvectors
    real(RDP), intent(in)  :: Lambda(:,:,:) !< 3D eigenvalues
    real(RDP), intent(in)  :: u(:,:,:,:)    !< operand
    real(RDP), intent(out) :: v(:,:,:,:)    !< result
    real(RDP), parameter :: alpha  = 1
    real(RDP), parameter :: beta   = 0

    !-----------------------------------------------------------------------------
    ! local variables

    type(LIBXSMM_DMMFunction), save :: xmm_1, xmm_2, xmm_3

    real(RDP) :: St( size(S,1), size(S,1) )
    real(RDP) :: z ( size(S,1), size(S,1), size(S,1) )
    integer   :: np, ne
    integer   :: e, k

    !---------------------------------------------------------------------------
    ! initialization

    np = size(S,1)
    ne = size(u,4)

    St = transpose(S)

    !-----------------------------------------------------------------------------
    ! XSMM dispatch

    !$omp master
    call LIBXSMM_Dispatch(xmm_1, np   , np**2, np, alpha=alpha, beta=beta)
    call LIBXSMM_Dispatch(xmm_2, np   , np   , np, alpha=alpha, beta=beta)
    call LIBXSMM_Dispatch(xmm_3, np**2, np   , np, alpha=alpha, beta=beta)

    if (.not. ( libxsmm_available(xmm_1) .and.  &
                libxsmm_available(xmm_2) .and.  &
                libxsmm_available(xmm_3)        ) ) then

      stop "TPO_Div_R_XSMM: LIBXSMM_Dispatch failed"

    end if
    !$omp end master
    !$omp barrier

    !---------------------------------------------------------------------------
    ! evaluation

    !$omp do
    do e = 1, ne

      ! z = SᵀxIxI uᵉ
      call LIBXSMM_DMMCall(xmm_3, u(:,:,:,e), S, z(:,:,:))

      ! vᵉ = IxSᵀxI z
      do k = 1, np
        call LIBXSMM_DMMCall(xmm_2, z(:,:,k), S, v(:,:,k,e))
      end do

      ! z  = IxIxSᵀ vᵉ
      call LIBXSMM_DMMCall(xmm_1, St, v(:,:,:,e), z(:,:,:))

      ! apply eigenvalues
      z = Lambda * z

      ! vᵉ = IxIxS z
      call LIBXSMM_DMMCall(xmm_1, S, z(:,:,:), v(:,:,:,e))

      ! z = IxSxI vᵉ
      do k = 1, np
        call LIBXSMM_DMMCall(xmm_2, v(:,:,k,e), St, z(:,:,k))
      end do

      ! vᵉ = SxIxI z
      call LIBXSMM_DMMCall(xmm_3, z(:,:,:), St, v(:,:,:,e))

    end do

  end subroutine TPO_Spectral_CI_XSMM_RDP

  !=============================================================================

#endif

end module TPO__Spectral__3D_CI__XSMM_RDP

