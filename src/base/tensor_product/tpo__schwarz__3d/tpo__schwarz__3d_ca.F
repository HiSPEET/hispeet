!> summary:  3d constant anisotropic Schwarz TPO
!> author:   Joerg Stiller
!> date:     2020/06/02
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Schwarz__3D_CA
  use Kind_Parameters, only: RDP
  use TPO__Schwarz__3D_CA__Gen
  use TPO__Schwarz__3D_CA__XSMM
  implicit none
  private

  public :: TPO_Schwarz

  interface TPO_Schwarz
    module procedure TPO_Schwarz_CA_RDP
  end interface

contains

  !-----------------------------------------------------------------------------
  !> Evaluates  u = [(W₃S₃ x W₂S₂ x W₁S₁) D⁻¹ (S₃ x S₂ x S₁)ᵀ] f
  !> where Wᵢ and Sᵢ depend on the subdomain configurations

  subroutine TPO_Schwarz_CA_RDP(S1, S2, S3, W1, W2, W3, cfg, D_inv, f, u)
    real(RDP), intent(in)  :: S1(:,:,:)      !< eigenvecs, dir 1   (n1,n1,nc)
    real(RDP), intent(in)  :: S2(:,:,:)      !< eigenvecs, dir 2   (n2,n2,nc)
    real(RDP), intent(in)  :: S3(:,:,:)      !< eigenvecs, dir 3   (n3,n3,nc)
    real(RDP), intent(in)  :: W1(:,:)        !< weights,   dir 1   (n1,nc)
    real(RDP), intent(in)  :: W2(:,:)        !< weights,   dir 2   (n2,nc)
    real(RDP), intent(in)  :: W3(:,:)        !< weights,   dir 3   (n3,nc)
    integer,   intent(in)  :: cfg(:,:)       !< subdomain configs  (3,nd)
    real(RDP), intent(in)  :: D_inv(:,:,:,:) !< inv 3D eigenvals   (n1,n2,n3,nd)
    real(RDP), intent(in)  :: f(:,:,:,:)     !< RHS                (n1,n2,n3,nd)
    real(RDP), intent(out) :: u(:,:,:,:)     !< solution           (n1,n2,n3,nd)

    select case(min( size(S1,1), size(S2,1), size(S3,1) ))

#ifdef __LIBXSMM__
    case(6:)
      call TPO_Schwarz_CA_XSMM(S1, S2, S3, W1, W2, W3, cfg, D_inv, f, u)
#endif

    case default
      call TPO_Schwarz_CA_Gen(S1, S2, S3, W1, W2, W3, cfg, D_inv, f, u)

    end select

  end subroutine TPO_Schwarz_CA_RDP

  !=============================================================================

end module TPO__Schwarz__3D_CA
