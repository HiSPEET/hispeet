module QOI__Sphere__3D
  use Kind_Parameters, only: RNP
  use Constants, only: ONE
  use QOI__Distribution__3D
  implicit none
  private

  public :: QOI_Sphere_3D

  type, extends(QOI_Distribution_3D) :: QOI_Sphere_3D
    real(RNP) :: cs = 100 !< radial scaling factor
    real(RNP) :: rs = 0.5 !< sphere radius
  contains
    procedure :: Density
  end type QOI_Sphere_3D

contains

  elemental real(RNP) function Density(this, x, y, z) result(f)
    class(QOI_Sphere_3D), intent(in) :: this
    real(RNP), intent(in) :: x !< x-coordinate
    real(RNP), intent(in) :: y !< y-coordinate
    real(RNP), intent(in) :: z !< z-coordinate

    real(RNP) :: r

    r = sqrt(x**2 + y**2 + z**2)
    f = ONE / ((this%cs * (r - this%rs))**2 + 1)

  end function Density

  !=============================================================================

end module QOI__Sphere__3D
