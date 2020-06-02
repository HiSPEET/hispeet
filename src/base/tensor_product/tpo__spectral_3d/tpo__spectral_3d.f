!> summary:  Evaluation of  v = [(S₃xS₂xS₁) Λ (S₃xS₂xS₁)ᵀ] u
!> author:   Joerg Stiller
!> date:     2020/05/29
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Spectral_3d
  use TPO__Spectral_3d_CI
  use TPO__Spectral_3d_CA
  implicit none
  private

  public :: TPO_Spectral

  interface TPO_Spectral
    module procedure TPO_Spectral_CI ! (S  x S  x S ) Λ (S  x S  x S )ᵀ
    module procedure TPO_Spectral_CA ! (S₃ x S₂ x S₁) Λ (S₃ x S₂ x S₁)ᵀ
  end interface

end module TPO__Spectral_3d
