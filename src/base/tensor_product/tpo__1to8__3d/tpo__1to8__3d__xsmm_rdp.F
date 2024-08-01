!> summary:   3d regular 1:8 interpolation operator using LIBXSMM for real(RDP)
!> author:    Joerg Stiller
!> date:      2024/07/29
!> license:   Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__1To8__3D__XSMM_RDP

#ifdef __LIBXSMM__

  use Kind_Parameters, only: RDP
  use LIBXSMM
  implicit none
  private

  public :: TPO_1To8_XSMM

  interface TPO_1To8_XSMM
    module procedure TPO_1To8_XSMM_RDP
  end interface

contains

  !-----------------------------------------------------------------------------
  !> XSMM AxAxA operator

  subroutine TPO_1To8_XSMM_RDP(A, u, v)
    real(RDP), intent(in)  :: A(:,:,:)   !< 1D operator
    real(RDP), intent(in)  :: u(:,:,:,:) !< operand
    real(RDP), intent(out) :: v(:,:,:,:) !< result
    real(RDP), parameter :: alpha  = 1
    real(RDP), parameter :: beta   = 0

    !---------------------------------------------------------------------------
    ! local variables

    type(LIBXSMM_DMMFunction), save :: xmm_1, xmm_2, xmm_3

    real(RDP) :: y ( size(A,2), size(A,1), size(A,1), 2, 2 )
    real(RDP) :: z ( size(A,2), size(A,2), size(A,1), 2 )
    real(RDP) :: At( size(A,2), size(A,1), 2)
    integer   :: na1, na2
    integer   :: e, k, o

    !---------------------------------------------------------------------------
    ! initialization

    na1 = size(A,1)
    na2 = size(A,2)

    At(:,:,1) = transpose(A(:,:,1))
    At(:,:,2) = transpose(A(:,:,2))

    !-----------------------------------------------------------------------------
    ! XSMM dispatch
    !$omp master
    call LIBXSMM_Dispatch(xmm_3, na2**2, na1 *2, na2, alpha=alpha, beta=beta)
    call LIBXSMM_Dispatch(xmm_2, na2   , na1   , na2, alpha=alpha, beta=beta)
    call LIBXSMM_Dispatch(xmm_1, na1   , na1**2, na2, alpha=alpha, beta=beta)

    if (.not. ( libxsmm_available(xmm_1) .and.  &
                libxsmm_available(xmm_2) .and.  &
                libxsmm_available(xmm_3)        ) ) then

      stop "TPO_1To8_XSMM_RDP: LIBXSMM_Dispatch failed"

    end if
    !$omp end master
    !$omp barrier

    !---------------------------------------------------------------------------
    ! evaluation

    !$omp do
    do e = 1, size(u,4)

      ! z = AxIxI uᵉ
      call LIBXSMM_DMMCall(xmm_3, u(:,:,:,e), At , z(:,:,:,:))

      ! y = IxAxI z
      do k = 1, na1
        call LIBXSMM_DMMCall(xmm_2, z(:,:,k,1), At(:,:,1), y(:,:,k,1,1))
        call LIBXSMM_DMMCall(xmm_2, z(:,:,k,1), At(:,:,2), y(:,:,k,2,1))
        call LIBXSMM_DMMCall(xmm_2, z(:,:,k,2), At(:,:,1), y(:,:,k,1,2))
        call LIBXSMM_DMMCall(xmm_2, z(:,:,k,2), At(:,:,2), y(:,:,k,2,2))
      end do

      ! offset of fine element index
      o = 8 * (e-1)

      ! vᵉ = IxIxA y
      call LIBXSMM_DMMCall(xmm_1, A(:,:,1), y(:,:,:,1,1), v(:,:,:,o+1))
      call LIBXSMM_DMMCall(xmm_1, A(:,:,2), y(:,:,:,1,1), v(:,:,:,o+2))
      call LIBXSMM_DMMCall(xmm_1, A(:,:,1), y(:,:,:,2,1), v(:,:,:,o+3))
      call LIBXSMM_DMMCall(xmm_1, A(:,:,2), y(:,:,:,2,1), v(:,:,:,o+4))
      call LIBXSMM_DMMCall(xmm_1, A(:,:,1), y(:,:,:,1,2), v(:,:,:,o+5))
      call LIBXSMM_DMMCall(xmm_1, A(:,:,2), y(:,:,:,1,2), v(:,:,:,o+6))
      call LIBXSMM_DMMCall(xmm_1, A(:,:,1), y(:,:,:,2,2), v(:,:,:,o+7))
      call LIBXSMM_DMMCall(xmm_1, A(:,:,2), y(:,:,:,2,2), v(:,:,:,o+8))

    end do

  end subroutine TPO_1To8_XSMM_RDP

  !=============================================================================

#endif

end module TPO__1To8__3D__XSMM_RDP
