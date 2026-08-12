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

!> summary:  Defines a class providing the density of a quantity or interest
!> author:   Joerg Stiller
!> date:     2024/12/22
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

