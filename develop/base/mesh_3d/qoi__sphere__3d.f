!------------------------------------------------------------------------------!
! This file is part of HiSPEET: High-order Spectral Element Techniques         !
!                                                                              !
! Copyright (C) 2026 by the HiSPEET authors and the                            !
! Chair of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany                 !
!                                                                              !
! HiSPEET is free software: you can redistribute it and/or modify              !
! it under the terms of the GNU General Public License as published by         !
! the Free Software Foundation, either version 3 of the License, or            !
! (at your option) any later version.                                          !
!                                                                              !
! HiSPEET is distributed in the hope that it will be useful,                   !
! but WITHOUT ANY WARRANTY; without even the implied warranty of               !
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                         !
! See the GNU General Public License for more details.                         !
!                                                                              !
! You should have received a copy of the GNU General Public License            !
! along with HiSPEET. If not, see <http://www.gnu.org/licenses/>.              !
!------------------------------------------------------------------------------!

!> summary:  QOI distribution based on a spherical shell
!> author:   Joerg Stiller
!> date:     2024/12/27
!===============================================================================

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
