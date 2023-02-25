!> summary:  3d generic anisotropic Schwarz operator using LIBXSMM
!> author:   Joerg Stiller, Erik Pfister
!> date:     2020/06/02
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Schwarz__3D_A__XSMM_RDP

#ifdef __LIBXSMM__

  use LIBXSMM
  use Kind_Parameters, only: RDP
  implicit none
  private

  public :: TPO_Schwarz_A_XSMM

  interface TPO_Schwarz_A_XSMM
    module procedure TPO_Schwarz_A_XSMM_RDP
  end interface

contains

  !-----------------------------------------------------------------------------
  !> 3d generic anisotropic Schwarz operator

  subroutine TPO_Schwarz_A_XSMM_RDP( S1, S2, S3, V1, V2, V3, W1, W2, W3, g, &
                                     cfg, lambda, nu, f, u                  )

    real(RDP), intent(in)  :: S1(:,:,:)  !< eigenvecs, dir 1   (n1,n1,nc)
    real(RDP), intent(in)  :: S2(:,:,:)  !< eigenvecs, dir 2   (n2,n2,nc)
    real(RDP), intent(in)  :: S3(:,:,:)  !< eigenvecs, dir 3   (n3,n3,nc)
    real(RDP), intent(in)  :: V1(:,:)    !< eigenvals, dir 1   (n1,nc)
    real(RDP), intent(in)  :: V2(:,:)    !< eigenvals, dir 2   (n2,nc)
    real(RDP), intent(in)  :: V3(:,:)    !< eigenvals, dir 3   (n3,nc)
    real(RDP), intent(in)  :: W1(:,:)    !< weights,   dir 1   (n1,nc)
    real(RDP), intent(in)  :: W2(:,:)    !< weights,   dir 2   (n2,nc)
    real(RDP), intent(in)  :: W3(:,:)    !< weights,   dir 3   (n3,nc)
    real(RDP), intent(in)  :: g(:,:)     !< metrics            (4,nd)
    integer,   intent(in)  :: cfg(:,:)   !< subdomain configs  (3,nd)
    real(RDP), intent(in)  :: lambda     !< λ
    real(RDP), intent(in)  :: nu(:)      !< ν + νˢ             (nd)
    real(RDP), intent(in)  :: f(:,:,:,:) !< RHS                (n1,n2,n3,nd)
    real(RDP), intent(out) :: u(:,:,:,:) !< solution           (n1,n2,n3,nd)

    real(RDP), parameter :: alpha = 1
    real(RDP), parameter :: beta  = 0
    real(RDP), parameter :: eps   = epsilon(lambda)

    !---------------------------------------------------------------------------
    ! local variables

    type(LIBXSMM_DMMFunction), save :: xmm_1, xmm_2, xmm_3

    real(RDP) :: S1_t  ( size(S1,2), size(S1,1), size(S1,3) )
    real(RDP) :: WS1   ( size(S1,1), size(S1,2), size(S1,3) )
    real(RDP) :: WS2_t ( size(S2,1), size(S2,2), size(S2,3) )
    real(RDP) :: WS3_t ( size(S3,1), size(S3,2), size(S3,3) )
    real(RDP) :: y     ( size(f ,1), size(f ,2), size(f ,3) )
    real(RDP) :: z     ( size(f ,1), size(f ,2), size(f ,3) )
    real(RDP) :: a1, a2, a3, a4, d

    integer :: n1, n2, n3, nc, nd
    integer :: i, j, k, l
    integer :: c, c1, c2, c3

    !-----------------------------------------------------------------------------
    ! initialization

    n1 = size(S1, 1)
    n2 = size(S2, 1)
    n3 = size(S3, 1)
    nc = size(S3, 3)
    nd = size(f , 4)

    ! WS1(:,:,c) = (W₁S₁)(:,:,c)
    do c = 1, nc
    do j = 1, n1
    do i = 1, n1
      WS1(i,j,c) = S1(i,j,c) * W1(i,c)
    end do
    end do
    end do

    ! WS2_t(:,:,c) = (W₂S₂)ᵀ(:,:,c)
    do c = 1, nc
    do j = 1, n2
    do i = 1, n2
      WS2_t(j,i,c) = W2(i,c) * S2(i,j,c)
    end do
    end do
    end do

    ! WS3_t(:,:,c) = (W₃S₃)ᵀ(:,:,c)
    do c = 1, nc
    do j = 1, n3
    do i = 1, n3
      WS3_t(j,i,c) = W3(i,c) * S3(i,j,c)
    end do
    end do
    end do

    do c = 1, nc
    do j = 1, n1
    do i = 1, n1
      S1_t(j,i,c) = S1(i,j,c)
    end do
    end do
    end do

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

    !$omp do private(l)
    subdomains: do l = 1, nd

      ! element-boundary configurations
      c1 = cfg(1,l)
      c2 = cfg(2,l)
      c3 = cfg(3,l)

      ! coefficients
      a1 = g(1,l) * nu(l)
      a2 = g(2,l) * nu(l)
      a3 = g(3,l) * nu(l)
      a4 = g(4,l) * lambda

      ! y = S₃ᵀ x I x I f ......................................................
      call LIBXSMM_DMMCall(xmm_3, f(:,:,:,l), S3(:,:,c3) , y(:,:,:))

      ! z = I x S₂ᵀ x I y ......................................................
      do k = 1, n3
        call LIBXSMM_DMMCall(xmm_2, y(:,:,k), S2(:,:,c2), z(:,:,k))
      end do

      ! y = I x I x S₁ᵀ z ......................................................
      call LIBXSMM_DMMCall(xmm_1, S1_t(:,:,c1), z(:,:,:) , y(:,:,:))

      ! y = D⁻¹ y
      do k = 1, n3
      do j = 1, n2
      do i = 1, n1
        d = a1 * V1(i,c1) + a2 * V2(j,c2) + a3 * V3(k,c3) + a4
        y(i,j,k) = y(i,j,k) / sign(max(abs(d),eps), d)
      end do
      end do
      end do

      ! z = I x I x W₁S₁ y .....................................................
      call LIBXSMM_DMMCall(xmm_1, WS1(:,:,c1), y(:,:,:) , z(:,:,:))

      ! y = I x W₂S₂ x I z .....................................................
      do k = 1, n3
        call LIBXSMM_DMMCall(xmm_2, z(:,:,k), WS2_t(:,:,c2), y(:,:,k))
      end do

      ! u = W₃S₃ x I x I y .....................................................
      call LIBXSMM_DMMCall(xmm_3, y(:,:,:), WS3_t(:,:,c3) , u(:,:,:,l))

    end do subdomains
    !$omp end do

  end subroutine TPO_Schwarz_A_XSMM_RDP

  !=============================================================================

#endif

end module TPO__Schwarz__3D_A__XSMM_RDP
