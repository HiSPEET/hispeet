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

!> summary:  Element operators for discontinuous cuboidal elements
!> author:   Joerg Stiller
!> date:     2016/03/18
!>
!>### Element operators for discontinuous cuboidal elements
!>
!> @note
!> Obsolete. Avoid using this type.
!> @endnote
!===============================================================================

module CART__DG_Element_Operators
  use Kind_Parameters,   only: RNP
  use Execution_Control, only: Error
  use DG__Element_Operators__1D
  implicit none
  private

  public :: DG_ElementOperators3D

  !-----------------------------------------------------------------------------
  !> Element operators for discontinuous cuboidal elements

  type, extends(DG_ElementOperators_1D) :: DG_ElementOperators3D

    real(RNP) :: dx(3) !< element extensions
    real(RNP) :: mu(3) !< penalty coefficients

  contains

    procedure :: Init_DG_ElementOperators3D
    procedure :: Get_1D_StiffnessMatrix

  end type DG_ElementOperators3D

contains

!-------------------------------------------------------------------------------
!> Initialize a new DG element operators

subroutine Init_DG_ElementOperators3D(this, po, dx, penalty)
  class(DG_ElementOperators3D), intent(inout) :: this
  integer,             intent(in) :: po       !< polynomial order
  real(RNP),           intent(in) :: dx(3)    !< element extensions
  real(RNP), optional, intent(in) :: penalty  !< penalty parameter > 1 [2]

  type(DG_ElementOptions_1D) :: opt

  opt % po = po
  if (present(penalty)) then
    opt % penalty = penalty
  end if

  call this % Init_DG_ElementOperators_1D(opt)

  this % dx    = dx
  this % mu(1) = this % PenaltyFactor([dx(1), dx(1)])
  this % mu(2) = this % PenaltyFactor([dx(2), dx(2)])
  this % mu(3) = this % PenaltyFactor([dx(3), dx(3)])

end subroutine Init_DG_ElementOperators3D

!-------------------------------------------------------------------------------
!> 1D stiffness matrix for given direction and boundary conditions
!>
!> The element stiffness matrix `Le` represents the nontrivial row entries
!> of the global 1D stiffness matrix corresponding to the given element. It
!> must be dimensioned as `Le(0:po,0:po,-1:1)`, where `po = this%po` is the
!> polynomial order. The third index refers to the preceding (-1), current (0)
!> and succeeding (1) element, respectively.

subroutine Get_1D_StiffnessMatrix(this, direction, bc, Le)
  class(DG_ElementOperators3D), intent(in) :: this
  integer,   intent(in)  :: direction     !< coordinate direction {1,2,3}
  character, intent(in)  :: bc(2)         !< boundary conditions {'','D','N'}
  real(RNP), intent(out) :: Le(0:,0:,-1:) !< 1D element stiffness matrix

  real(RNP) :: dx(-1:1)

  select case(direction)
  case(1:3)
    dx = this % dx(direction)
    call this % Get_StiffnessMatrix(dx, bc, Le)
  case default
    call Error('Get_1D_StiffnessMatrix',    &
               'direction must 1,2 or 3',   &
               'CART__DG_Element_Operators' )
  end select

end subroutine Get_1D_StiffnessMatrix

!===============================================================================

end module CART__DG_Element_Operators
