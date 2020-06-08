!> summary:  3d generic anisotropic spectral operator using LIBXSMM
!> author:   Erik Pfister, Joerg Stiller
!> date:     2020/06/02
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Spectral_3d_CA__XSMM

#ifdef __LIBXSMM__

  use LIBXSMM
  use Kind_Parameters, only: RNP
  implicit none
  private

  public :: TPO_Spectral_CA_XSMM

contains

  !-----------------------------------------------------------------------------
  !> Generic procedure for evaluating  v = [(S₃ x S₂ x S₁) Λ (S₃ x S₂ x S₁)ᵀ] u

  subroutine TPO_Spectral_CA_XSMM(S1, S2, S3, Lambda, u, v)
    real(RNP), intent(in)  :: S1(:,:)       !< eigenvectors for dir 1: S₁(n1,n1)
    real(RNP), intent(in)  :: S2(:,:)       !< eigenvectors for dir 2: S₂(n2,n2)
    real(RNP), intent(in)  :: S3(:,:)       !< eigenvectors for dir 3: S₃(n3,n3)
    real(RNP), intent(in)  :: Lambda(:,:,:) !< 3D eigenvalues: Λ(n1,n2,n3)
    real(RNP), intent(in)  :: u(:,:,:,:)    !< operand: u(n1.n2,n3,ne)
    real(RNP), intent(out) :: v(:,:,:,:)    !< result:  v(n1,n2,n3,ne)
    real(RNP), parameter :: alpha  = 1
    real(RNP), parameter :: beta   = 0

    !-----------------------------------------------------------------------------
    ! local variables

    type(LIBXSMM_DMMFunction), save :: xmm_1, xmm_2, xmm_3

    real(RNP) :: z( size(S1,1), size(S2,1), size(S3,1) )
    real(RNP) :: tmp
    integer   :: n1, n2, n3, ne
    integer   :: e, i, j, k, p    
    real(RNP) :: S1_t (size(S1,1), size(S1,1) )
    real(RNP) :: S2_t (size(S2,1), size(S2,1) )
    real(RNP) :: S3_t (size(S3,1), size(S3,1) )

    !---------------------------------------------------------------------------
    ! initialization

    n1 = size(S1,1)
    n2 = size(S2,1)
    n3 = size(S3,1)
    ne = size(u ,4)

    S1_t = transpose(S1)
    S2_t = transpose(S2)
    S3_t = transpose(S3)

    !-----------------------------------------------------------------------------
    ! XSMM dispatch

    !$omp master
    call LIBXSMM_Dispatch(xmm_1, n1   , n2*n3, n1, alpha=alpha, beta=beta)
    call LIBXSMM_Dispatch(xmm_2, n1   , n2   , n2, alpha=alpha, beta=beta)
    call LIBXSMM_Dispatch(xmm_3, n1*n2, n3   , n3, alpha=alpha, beta=beta)

    if (.not. ( libxsmm_available(xmm_1) .and.  &
                libxsmm_available(xmm_2) .and.  &
                libxsmm_available(xmm_3)        ) ) then

      stop "TPO_Div_R_XSMM: LIBXSMM_Dispatch failed"

    end if
    !$omp end master
    !$omp barrier

    !---------------------------------------------------------------------------
    ! evaluation

    !$omp do private(e)
    do e = 1, ne

      ! z = S₃ᵀ x I x I u ......................................................
      call LIBXSMM_DMMCall(xmm_3, u(:,:,:,e), S3(:,:), z(:,:,:))

      ! vᵉ = I x S₂ᵀ x I z .....................................................
      do k = 1, n3
        call LIBXSMM_DMMCall(xmm_2, z(:,:,k), S2(:,:), v(:,:,k,e))
      end do      

      ! z = I x I x S₁ᵀ vᵉ .....................................................
      call LIBXSMM_DMMCall(xmm_1, S1_t(:,:), v(:,:,:,e), z(:,:,:))

      ! apply eigenvalues
      z = Lambda * z

      ! vᵉ = I x I x S₁ z ......................................................
      call LIBXSMM_DMMCall(xmm_1, S1(:,:), z(:,:,:), v(:,:,:,e))

      ! z = I x S₂ x I vᵉ ......................................................
      do k = 1, n3
        call LIBXSMM_DMMCall(xmm_2, v(:,:,k,e), S2_t(:,:), z(:,:,k))
      end do     
      
      ! vᵉ = S₃ x I x I z ......................................................
      call LIBXSMM_DMMCall(xmm_3, z(:,:,:), S3_t(:,:), v(:,:,:,e))      

    end do

  end subroutine TPO_Spectral_CA_XSMM

  !=============================================================================

#endif

end module TPO__Spectral_3d_CA__XSMM
