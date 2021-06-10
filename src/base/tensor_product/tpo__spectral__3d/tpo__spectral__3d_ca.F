!> summary:  Evaluation of  v = [(S₃xS₂xS₁) Λ (S₃xS₂xS₁)ᵀ] u
!>           with 3D mesh variables
!> author:   Joerg Stiller
!> date:     2020/06/02
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Spectral__3D_CA
  use Kind_Parameters, only: RDP
  use TPO__Spectral__3D_CA__Gen
  use TPO__Spectral__3D_CA__XSMM
  implicit none
  private

  public :: TPO_Spectral_CA

  interface TPO_Spectral_CA
    module procedure TPO_Spectral_CA_RDP
  end interface

contains

  !-----------------------------------------------------------------------------
  !> 3D (S₃ x S₂ x S₁) Λ (S₃ x S₂ x S₁)ᵀ operator

  subroutine TPO_Spectral_CA_RDP(S1, S2, S3, Lambda, u, v)
    real(RDP), intent(in)  :: S1(:,:)       !< eigenvectors for dir 1: S₁(n1,n1)
    real(RDP), intent(in)  :: S2(:,:)       !< eigenvectors for dir 2: S₂(n2,n2)
    real(RDP), intent(in)  :: S3(:,:)       !< eigenvectors for dir 3: S₃(n3,n3)
    real(RDP), intent(in)  :: Lambda(:,:,:) !< 3D eigenvalues: Λ(n1,n2,n3)
    real(RDP), intent(in)  :: u(:,:,:,:)    !< operand: u(n1.n2,n3,ne)
    real(RDP), intent(out) :: v(:,:,:,:)    !< result:  v(n1,n2,n3,ne)

    select case(min( size(S1,1), size(S2,1), size(S3,1) ))

#ifdef __LIBXSMM__
    case(6:)
      call TPO_Spectral_CA_XSMM(S1, S2, S3, Lambda, u, v)
#endif

    case default
      call TPO_Spectral_CA_Gen(S1, S2, S3, Lambda, u, v)

    end select

  end subroutine TPO_Spectral_CA_RDP

  !=============================================================================

end module TPO__Spectral__3D_CA
