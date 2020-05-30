!> summary:  Evaluation of  v = [(S₃xS₂xS₁) Λ (S₃xS₂xS₁)ᵀ] u
!> author:   Joerg Stiller
!> date:     2020/05/29
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Spectral_3d
  use TPO__Spectral_3d_I
  implicit none
  private

  public :: TPO_Spectral

  interface TPO_Spectral
    module procedure TPO_Spectral_I ! (SxSxS) Λ (SxSxS)ᵀ
   !module procedure TPO_Spectral_A ! (S₃xS₂xS₁) Λ (S₃xS₂xS₁)ᵀ
  end interface

end module TPO__Spectral_3d
