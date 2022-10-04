!> summary:  3d constant anisotropic Schwarz TPO
!> author:   Joerg Stiller
!> date:     2020/06/02
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Schwarz__3D_A
  use Kind_Parameters, only: RSP, RDP
  use TPO__Schwarz__3D_A__Gen
  use TPO__Schwarz__3D_A__XSMM
  implicit none
  private

  public :: TPO_Schwarz

  interface TPO_Schwarz
    module procedure TPO_Schwarz_A_RSP
    module procedure TPO_Schwarz_A_RDP
  end interface

contains

  !-----------------------------------------------------------------------------
  !> Evaluates  u = [(W₃S₃ x W₂S₂ x W₁S₁) D⁻¹ (S₃ x S₂ x S₁)ᵀ] f,
  !> where Wᵢ and Sᵢ depend on the subdomain configurations, for real(RSP)

  subroutine TPO_Schwarz_A_RSP( S1, S2, S3, V1, V2, V3, W1, W2, W3, g, cfg, &
                                lambda, nu, f, u                            )

    real(RSP), intent(in)  :: S1(:,:,:)  !< eigenvecs, dir 1   (n1,n1,nc)
    real(RSP), intent(in)  :: S2(:,:,:)  !< eigenvecs, dir 2   (n2,n2,nc)
    real(RSP), intent(in)  :: S3(:,:,:)  !< eigenvecs, dir 3   (n3,n3,nc)
    real(RSP), intent(in)  :: V1(:,:)    !< eigenvals, dir 1   (n1,nc)
    real(RSP), intent(in)  :: V2(:,:)    !< eigenvals, dir 2   (n2,nc)
    real(RSP), intent(in)  :: V3(:,:)    !< eigenvals, dir 3   (n3,nc)
    real(RSP), intent(in)  :: W1(:,:)    !< weights,   dir 1   (n1,nc)
    real(RSP), intent(in)  :: W2(:,:)    !< weights,   dir 2   (n2,nc)
    real(RSP), intent(in)  :: W3(:,:)    !< weights,   dir 3   (n3,nc)
    real(RSP), intent(in)  :: g(:,:)     !< metrics            (4,nd)
    integer,   intent(in)  :: cfg(:,:)   !< subdomain configs  (3,nd)
    real(RSP), intent(in)  :: lambda     !< λ
    real(RSP), intent(in)  :: nu(:)      !< ν + νˢ             (nd)
    real(RSP), intent(in)  :: f(:,:,:,:) !< RHS                (n1,n2,n3,nd)
    real(RSP), intent(out) :: u(:,:,:,:) !< solution           (n1,n2,n3,nd)

    select case(min( size(S1,1), size(S2,1), size(S3,1) ))

#ifdef __LIBXSMM__
    case(6:)
      call TPO_Schwarz_A_XSMM( S1, S2, S3, V1, V2, V3, W1, W2, W3, g, cfg, &
                               lambda, nu, f, u                            )
#endif

    case default
      call TPO_Schwarz_A_Gen( S1, S2, S3, V1, V2, V3, W1, W2, W3, g, cfg, &
                              lambda, nu, f, u                            )

    end select

  end subroutine TPO_Schwarz_A_RSP

  !-----------------------------------------------------------------------------
  !> Evaluates  u = [(W₃S₃ x W₂S₂ x W₁S₁) D⁻¹ (S₃ x S₂ x S₁)ᵀ] f,
  !> where Wᵢ and Sᵢ depend on the subdomain configurations, for real(RDP)

  subroutine TPO_Schwarz_A_RDP( S1, S2, S3, V1, V2, V3, W1, W2, W3, g, cfg, &
                                lambda, nu, f, u                            )

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

    select case(min( size(S1,1), size(S2,1), size(S3,1) ))

#ifdef __LIBXSMM__
    case(6:)
      call TPO_Schwarz_A_XSMM( S1, S2, S3, V1, V2, V3, W1, W2, W3, g, cfg, &
                               lambda, nu, f, u                            )
#endif

    case default
      call TPO_Schwarz_A_Gen( S1, S2, S3, V1, V2, V3, W1, W2, W3, g, cfg, &
                              lambda, nu, f, u                            )

    end select

  end subroutine TPO_Schwarz_A_RDP

  !=============================================================================

end module TPO__Schwarz__3D_A
