!> summary:  Evaluation of 3D constant isotropic Schwarz TPO
!> author:   Joerg Stiller
!> date:     2020/06/02
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Schwarz__3D_CI
  use Kind_Parameters, only: RSP, RDP
  use TPO__Schwarz__3D_CI__Hand
  use TPO__Schwarz__3D_CI__XSMM
  implicit none
  private

  public :: TPO_Schwarz

  interface TPO_Schwarz
    module procedure TPO_Schwarz_CI_RSP
    module procedure TPO_Schwarz_CI_RDP
  end interface

contains

  !-----------------------------------------------------------------------------
  !> 3D constant isotropic operator for real(RSP)

  subroutine TPO_Schwarz_CI_RSP(S, W, cfg, D_inv, f, u)
    real(RSP), intent(in)  :: S(:,:,:)       !< 1D eigenvectors
    real(RSP), intent(in)  :: W(:,:)         !< 1D weights
    integer,   intent(in)  :: cfg(:,:)       !< subdomain configurations
    real(RSP), intent(in)  :: D_inv(:,:,:,:) !< inv. 3D eigenvalues
    real(RSP), intent(in)  :: f(:,:,:,:)     !< RHS
    real(RSP), intent(out) :: u(:,:,:,:)     !< solution

    select case(size(S,1))

#ifdef __LIBXSMM__
    case(8:)
      call TPO_Schwarz_CI_XSMM(S, W, cfg, D_inv, f, u)
#endif

    case default
      call TPO_Schwarz_CI_Hand(S, W, cfg, D_inv, f, u)

    end select

  end subroutine TPO_Schwarz_CI_RSP

  !-----------------------------------------------------------------------------
  !> 3D constant isotropic operator for real(RDP)

  subroutine TPO_Schwarz_CI_RDP(S, W, cfg, D_inv, f, u)
    real(RDP), intent(in)  :: S(:,:,:)       !< 1D eigenvectors
    real(RDP), intent(in)  :: W(:,:)         !< 1D weights
    integer,   intent(in)  :: cfg(:,:)       !< subdomain configurations
    real(RDP), intent(in)  :: D_inv(:,:,:,:) !< inv. 3D eigenvalues
    real(RDP), intent(in)  :: f(:,:,:,:)     !< RHS
    real(RDP), intent(out) :: u(:,:,:,:)     !< solution

    select case(size(S,1))

#ifdef __LIBXSMM__
    case(8:)
      call TPO_Schwarz_CI_XSMM(S, W, cfg, D_inv, f, u)
#endif

    case default
      call TPO_Schwarz_CI_Hand(S, W, cfg, D_inv, f, u)

    end select

  end subroutine TPO_Schwarz_CI_RDP

  !=============================================================================

end module TPO__Schwarz__3D_CI
