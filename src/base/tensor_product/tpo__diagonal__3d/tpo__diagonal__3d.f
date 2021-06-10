!> summary:  3D scaled diagonal operator
!> author:   Joerg Stiller
!> date:     2020/06/04
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Diagonal__3D
  use TPO__Diagonal__3D_CI
  implicit none
  private

  public :: TPO_Diagonal

  interface TPO_Diagonal
    module procedure TPO_Diagonal_CI_RDP
  end interface

end module TPO__Diagonal__3D
