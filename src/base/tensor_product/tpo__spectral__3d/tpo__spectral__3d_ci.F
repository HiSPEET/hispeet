!> summary:  Evaluation of  v = [(SxSxS) Λ (SxSxS)ᵀ] u  with 3D mesh variables
!> author:   Joerg Stiller
!> date:     2020/05/29
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Spectral__3D_CI
  use Kind_Parameters, only: RDP
  use TPO__Spectral__3D_CI__Hand
  use TPO__Spectral__3D_CI__XSMM
  implicit none
  private

  public :: TPO_Spectral_CI

  interface TPO_Spectral_CI
    module procedure TPO_Spectral_CI_RDP
  end interface

contains

  !-----------------------------------------------------------------------------
  !> 3D (SxSxS) Λ (SxSxS)ᵀ operator

  subroutine TPO_Spectral_CI_RDP(S, Lambda, u, v)
    real(RDP), intent(in)  :: S(:,:)        !< 1D eigenvectors
    real(RDP), intent(in)  :: Lambda(:,:,:) !< 3D eigenvalues
    real(RDP), intent(in)  :: u(:,:,:,:)    !< operand
    real(RDP), intent(out) :: v(:,:,:,:)    !< result

    select case(size(S,1))

#ifdef __LIBXSMM__
    case(8:)
      call TPO_Spectral_CI_XSMM(S, Lambda, u, v)
#endif

    case default
      call TPO_Spectral_CI_Hand(S, Lambda, u, v)

    end select

  end subroutine TPO_Spectral_CI_RDP

  !=============================================================================

end module TPO__Spectral__3D_CI
