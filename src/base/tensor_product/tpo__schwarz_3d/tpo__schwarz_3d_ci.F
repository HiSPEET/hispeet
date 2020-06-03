!> summary:  Evaluation of 3D constant isotropic Schwarz TPO
!> author:   Joerg Stiller
!> date:     2020/06/02
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Schwarz_3d_CI
  use Kind_Parameters, only: RNP
  use TPO__Schwarz_3d_CI__Hand
  use TPO__Schwarz_3d_CI__XSMM
  implicit none
  private

  public :: TPO_Schwarz_CI

contains

  !-----------------------------------------------------------------------------
  !> 3D constant isotropic operator

  subroutine TPO_Schwarz_CI(S, W, cfg, D_inv, f, u)
    real(RNP), intent(in)  :: S(:,:,:)       !< 1D eigenvectors
    real(RNP), intent(in)  :: W(:,:)         !< 1D weights
    integer,   intent(in)  :: cfg(:,:)       !< subdomain configurations
    real(RNP), intent(in)  :: D_inv(:,:,:,:) !< inv. 3D eigenvalues
    real(RNP), intent(in)  :: f(:,:,:,:)     !< RHS
    real(RNP), intent(out) :: u(:,:,:,:)     !< solution

    select case(size(S,1))

#ifdef __LIBXSMM__
    case(8:)
      call TPO_Schwarz_CI_XSMM(S, W, cfg, D_inv, f, u)
#endif

    case default
      call TPO_Schwarz_CI_Hand(S, W, cfg, D_inv, f, u)

    end select

  end subroutine TPO_Schwarz_CI

  !=============================================================================

end module TPO__Schwarz_3d_CI
