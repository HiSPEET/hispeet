!> summary:   3d generic AxAxA operator using LIBXSMM
!> author:    Erik Pfister, Joerg Stiller
!> date:      2020/06/01
!> license:   Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> @todo      Replace this dummy by true LIBXSMM version
!===============================================================================

module TPO__AAA_3d__XSMM

#ifdef __LIBXSMM__

  use Kind_Parameters, only: RNP
  use LIBXSMM
  implicit none
  private

  public :: TPO_AAA_XSMM

contains

  !-----------------------------------------------------------------------------
  !> XSMM AxAxA operator

  subroutine TPO_AAA_XSMM(A, u, v)
    real(RNP), intent(in)  :: A(:,:)     !< 1D operator
    real(RNP), intent(in)  :: u(:,:,:,:) !< operand
    real(RNP), intent(out) :: v(:,:,:,:) !< result
    real(RNP), parameter :: alpha  = 1
    real(RNP), parameter :: beta   = 0

    !---------------------------------------------------------------------------
    ! local variables

    type(LIBXSMM_DMMFunction), save :: xmm_1, xmm_2, xmm_3

    real(RNP) :: At( size(A,2), size(A,1) )
    real(RNP) :: z2( size(A,2), size(A,1), size(A,1))
    real(RNP) :: z3( size(A,2), size(A,2), size(A,1))
    real(RNP) :: tmp
    integer   :: na1, na2, ne
    integer   :: e, k

    !---------------------------------------------------------------------------
    ! initialization

    na1 = size(A,1)
    na2 = size(A,2)
    ne  = size(u,4)

    At = transpose(A)

    !-----------------------------------------------------------------------------
    ! XSMM dispatch
    !$omp master
    call LIBXSMM_Dispatch(xmm_1, na2**2, na1   , na2, alpha=alpha, beta=beta)
    call LIBXSMM_Dispatch(xmm_2, na2   , na1   , na2, alpha=alpha, beta=beta)
    call LIBXSMM_Dispatch(xmm_3, na2   , na1**2, na2, alpha=alpha, beta=beta)

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

      ! z₃ = AxIxI uᵉ
      call LIBXSMM_DMMCall(xmm_3, u(:,:,:,e), At , z3(:,:,:))

      ! z₂ = IxAxI z₃
      do k = 1, na1
        call LIBXSMM_DMMCall(xmm_2, z3(:,:,k), At, z2(:,:,k))
      end do

      ! vᵉ = IxIxA z₂
      call LIBXSMM_DMMCall(xmm_1, A, z2(:,:,:), v(:,:,:,e))

    end do

  end subroutine TPO_AAA_XSMM

  !=============================================================================

#endif

end module TPO__AAA_3d__XSMM
