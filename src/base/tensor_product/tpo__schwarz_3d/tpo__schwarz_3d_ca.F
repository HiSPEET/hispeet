!> summary:  3d constant anisotropic Schwarz TPO
!> author:   Joerg Stiller
!> date:     2020/06/02
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Schwarz_3d_CA
  use Kind_Parameters, only: RNP
  use TPO__Schwarz_3d_CA__Gen
  use TPO__Schwarz_3d_CA__XSMM
  implicit none
  private

  public :: TPO_Schwarz_CA

contains

  !-----------------------------------------------------------------------------
  !> Evaluates  u = [(W₃S₃ x W₂S₂ x W₁S₁) D⁻¹ (S₃ x S₂ x S₁)ᵀ] f
  !> where Wᵢ and Sᵢ depend on the subdomain configurations

  subroutine TPO_Schwarz_CA(S1, S2, S3, W1, W2, W3, cfg, D_inv, f, u)
    real(RNP), intent(in)  :: S1(:,:,:)      !< eigenvecs, dir 1   (n1,n1,nc)
    real(RNP), intent(in)  :: S2(:,:,:)      !< eigenvecs, dir 2   (n2,n2,nc)
    real(RNP), intent(in)  :: S3(:,:,:)      !< eigenvecs, dir 3   (n3,n3,nc)
    real(RNP), intent(in)  :: W1(:,:)        !< weights,   dir 1   (n1,nc)
    real(RNP), intent(in)  :: W2(:,:)        !< weights,   dir 2   (n2,nc)
    real(RNP), intent(in)  :: W3(:,:)        !< weights,   dir 3   (n3,nc)
    integer,   intent(in)  :: cfg(:,:)       !< subdomain configs  (3,nd)
    real(RNP), intent(in)  :: D_inv(:,:,:,:) !< inv 3D eigenvals   (n1,n2,n3,nd)
    real(RNP), intent(in)  :: f(:,:,:,:)     !< RHS                (n1,n2,n3,nd)
    real(RNP), intent(out) :: u(:,:,:,:)     !< solution           (n1,n2,n3,nd)

    select case(min( size(S1,1), size(S2,1), size(S3,1) ))

#ifdef __LIBXSMM__
    case(6:)
      call TPO_Schwarz_CA_XSMM(S1, S2, S3, W1, W2, W3, cfg, D_inv, f, u)
#endif

    case default
      call TPO_Schwarz_CA_Gen(S1, S2, S3, W1, W2, W3, cfg, D_inv, f, u)

    end select

  end subroutine TPO_Schwarz_CA

  !=============================================================================

end module TPO__Schwarz_3d_CA
