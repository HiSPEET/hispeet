!> summary:  Defines a class providing the density of a quantity or interest
!> author:   Joerg Stiller
!> date:     2024/12/22
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module QOI__Distribution__3D
  use Kind_Parameters, only: RNP
  implicit none
  private

  public :: QOI_Distribution_3D

  !-----------------------------------------------------------------------------
  !> Class providing the density of a quantity or interest

  type, abstract :: QOI_Distribution_3D
  contains
    procedure(Density), deferred :: Density
  end type QOI_Distribution_3D

  abstract interface

    !---------------------------------------------------------------------------
    !> Density distribution of the quantity of interest

    elemental real(RNP) function Density(this, x, y, z) result(f)
      import :: QOI_Distribution_3D, RNP
      class(QOI_Distribution_3D), intent(in) :: this
      real(RNP), intent(in) :: x !< x-coordinate
      real(RNP), intent(in) :: y !< y-coordinate
      real(RNP), intent(in) :: z !< z-coordinate
    end function Density

  end interface

  !=============================================================================

end module QOI__Distribution__3D

