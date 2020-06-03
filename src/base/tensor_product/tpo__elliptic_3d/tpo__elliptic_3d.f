!> summary:  3D Cartesian elliptic element operators
!> author:   Joerg Stiller
!> date:     2020/05/29
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Elliptic_3d
  use TPO__Elliptic_3d_RLCI
  use TPO__Elliptic_3d_RLVI
  implicit none
  private

  public :: TPO_Elliptic

  interface TPO_Elliptic
    module procedure TPO_Elliptic_RLCI
    module procedure TPO_Elliptic_RLVI
  end interface

end module TPO__Elliptic_3d
