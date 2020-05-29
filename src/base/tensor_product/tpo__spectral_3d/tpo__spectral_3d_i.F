!> summary:  Evaluation of  v = [(S₃xS₂xS₁) Λ (S₃xS₂xS₁)ᵀ] u
!>           with 3D mesh variables
!> author:   Joerg Stiller
!> date:     2020/05/29
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Spectral_3d_I
  use Kind_Parameters, only: RNP
  use TPO__Spectral_3d_I__Hand
  use TPO__Spectral_3d_I__XSMM
  implicit none
  private

  public :: TPO_Spectral_I

contains

  !-----------------------------------------------------------------------------
  !> 3D (SxSxS) Λ (SxSxS)ᵀ operator

  subroutine TPO_Spectral_I(S, Lambda, u, v)
    use Kind_Parameters,  only: RNP
    real(RNP), intent(in)  :: S(:,:)        !< 1D eigenvectors
    real(RNP), intent(in)  :: Lambda(:,:,:) !< 3D eigenvalues
    real(RNP), intent(in)  :: u(:,:,:,:)    !< operand
    real(RNP), intent(out) :: v(:,:,:,:)    !< result

    select case(size(S,1))

#ifdef __LIBXSMM__
    case(8:)
      call TPO_Spectral_I_XSMM(S, Lambda, u, v)
#endif

    case default
      call TPO_Spectral_I_Hand(S, Lambda, u, v)

    end select

  end subroutine TPO_Spectral_I

  !=============================================================================

end module TPO__Spectral_3d_I
