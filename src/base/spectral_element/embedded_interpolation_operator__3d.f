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

!> summary:  Embedded interpolation of 3D mesh variables
!> author:   Joerg Stiller
!> date:     2018/11/02
!===============================================================================

module Embedded_Interpolation_Operator__3D
  use Kind_Parameters, only: RNP
  use Standard_Element_Operators__1D
  use Embedded_Interpolation_Operator__1D
  use TPO__AAA__3D
  implicit none
  private

  !-----------------------------------------------------------------------------
  !> Embedded interpolation of 3D mesh variables

  type, extends(EmbeddedInterpolationOperator_1D), &
      public :: EmbeddedInterpolationOperator_3D

  contains

    generic :: Apply  =>  Interpolate_S, Interpolate_A
    procedure, private :: Interpolate_S
    procedure, private :: Interpolate_A

  end type EmbeddedInterpolationOperator_3D

  ! constructor interface
  interface EmbeddedInterpolationOperator_3D
    module procedure New_SX
  end interface

contains

  !=============================================================================
  ! Constructor

  !-----------------------------------------------------------------------------
  !> New embedded 3D interpolation operator

  type(EmbeddedInterpolationOperator_3D) function New_SX(eop, xi) result(this)
    class(StandardElementOperators_1D), intent(in) :: eop !< standard operators
    real(RNP), intent(in) :: xi(:) !< points in [-1,1]

    call this % Init_EmbeddedInterpolationOperator_1D(eop, xi)

  end function New_SX

  !=============================================================================
  ! Type-bound procedures

  !-----------------------------------------------------------------------------
  !> Interpolation of scalar variables

  subroutine Interpolate_S(this, uo, ui)
    class(EmbeddedInterpolationOperator_3D), intent(in) :: this
    real(RNP), intent(in)  :: uo(:,:,:,:) !< original mesh variable
    real(RNP), intent(out) :: ui(:,:,:,:) !< interpolated mesh variable

    call TPO_AAA(this%A, uo, ui)

  end subroutine Interpolate_S

  !-----------------------------------------------------------------------------
  !> Interpolation of array variables

  subroutine Interpolate_A(this, uo, ui)
    class(EmbeddedInterpolationOperator_3D), intent(in) :: this
    real(RNP), intent(in)  :: uo(:,:,:,:,:) !< original mesh variable
    real(RNP), intent(out) :: ui(:,:,:,:,:) !< interpolated mesh variable

    integer :: c

    do c = 1, size(uo,5)
      call TPO_AAA(this%A, uo(:,:,:,:,c), ui(:,:,:,:,c))
    end do

  end subroutine Interpolate_A

  !=============================================================================

end module Embedded_Interpolation_Operator__3D
