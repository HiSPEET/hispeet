module QOI__Point__3D
  use Kind_Parameters, only: RNP
  use Constants, only: ONE
  use QOI__Distribution__3D
  implicit none
  private

  public :: QOI_Point_3D

  type, extends(QOI_Distribution_3D) :: QOI_Point_3D
    real(RNP) :: xp = 0 !< x-coordinate
    real(RNP) :: yp = 0 !< y-coordinate
    real(RNP) :: zp = 0 !< z-coordinate
  contains
    procedure :: Density
  end type QOI_Point_3D

contains

  elemental real(RNP) function Density(this, x, y, z) result(f)
    class(QOI_Point_3D), intent(in) :: this
    real(RNP), intent(in) :: x !< x-coordinate
    real(RNP), intent(in) :: y !< y-coordinate
    real(RNP), intent(in) :: z !< z-coordinate

    real(RNP) :: r

    r = sqrt((x - this%xp)**2 + (y - this%yp)**2 + (z - this%zp)**2)

    if (r < epsilon(r)) then
      f = 1
    else
      f = 0
    end if

  end function Density

  !=============================================================================

end module QOI__Point__3D
