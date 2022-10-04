!> summary:  3D scaled averaging operator
!> author:   Joerg Stiller
!> date:     2022/10/02
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Average__3D
  use TPO__Average__3D_I
  implicit none
  private

  public :: TPO_Average

  interface TPO_Average
    module procedure TPO_Average_I_RDP
  end interface

end module TPO__Average__3D
