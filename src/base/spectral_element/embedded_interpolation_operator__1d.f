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

!> summary:  Embedded interpolation operator
!> author:   Joerg Stiller
!> date:     2018/03/09
!>
!> Provides one-dimensional interpolation operators for the following nodal
!> bases
!>
!>   - `'G'`   Lagrangian with Gauss points
!>   - `'L'`   Lagrangian with on Lobatto points
!>   - `'RL'`  Lagrangian with left-sided Radau points
!>   - `'RR'`  Lagrangian with right-sided Radau points
!>   - `'N'`   Lagrangian with arbitrary nodes
!>
!===============================================================================

module Embedded_Interpolation_Operator__1D
  use Kind_Parameters, only: RNP
  use Execution_Control, only: Error
  use Gauss_Jacobi
  use Lagrange_Interpolation
  use Standard_Element_Operators__1D
  implicit none
  private

  !-----------------------------------------------------------------------------
  !> Embedded interpolation operator

  type, public :: EmbeddedInterpolationOperator_1D

    integer :: no = -1 !< number of original points per direction
    integer :: ni = -1 !< number of interpolated points per direction
    real(RNP), allocatable :: A(:,:)  !< 1D interpolation operator

  contains

    generic :: Init_EmbeddedInterpolationOperator_1D => Init_Points, Init_StdOps
    procedure, private :: Init_Points
    procedure, private :: Init_StdOps

  end type EmbeddedInterpolationOperator_1D

  ! constructor interface
  interface EmbeddedInterpolationOperator_1D
    module procedure New_Points
    module procedure New_StdOps
  end interface

contains

  !=============================================================================
  ! Constructors

  !-----------------------------------------------------------------------------
  !> New interpolation operator from given nodes and interpolation points

  type(EmbeddedInterpolationOperator_1D) function New_Points(nodes, xo, xi) &
      result(this)

    character(*), intent(in) :: nodes  !< nodal basis type
    real(RNP),    intent(in) :: xo(0:) !< node coordinates
    real(RNP),    intent(in) :: xi(1:) !< interpolation points in [-1,1]

    call Init_Points(this, nodes, xo, xi)

  end function New_Points

  !-----------------------------------------------------------------------------
  !> New interpolation operator from standard operators and interpolation points

  type(EmbeddedInterpolationOperator_1D) function New_StdOps(eop, xi) &
      result(this)

    class(StandardElementOperators_1D), intent(in) :: eop !< standard operators
    real(RNP), intent(in) :: xi(1:) !< interpolation points in [-1,1]

    call Init_StdOps(this, eop, xi)

  end function New_StdOps

  !=============================================================================
  ! Type-bound procedures

  !-----------------------------------------------------------------------------
  !> Initialization for given nodes and interpolation points

  subroutine Init_Points(this, nodes, xo, xi)
    class(EmbeddedInterpolationOperator_1D), intent(inout) :: this
    character(*), intent(in) :: nodes  !< nodal basis type
    real(RNP),    intent(in) :: xo(0:) !< node coordinates
    real(RNP),    intent(in) :: xi(1:) !< interpolation points in [-1,1]

    integer :: j, k, po

    if (allocated(this % A)) deallocate(this % A)

    associate(no => this % no, ni => this % ni)

      po = ubound(xo, 1)
      no = size(xo)
      ni = size(xi)

      allocate(this % A(ni,0:po))

      select case(nodes)
      case('G') ! Gauss
        do k = 0, po
        do j = 1, ni
          this % A(j,k) = GaussPolynomial(k, xo, xi(j))
        end do
        end do
      case('RL','RR') ! Radau left or right
        do k = 0, po
        do j = 1, ni
          this % A(j,k) = RadauPolynomial(k, xo, xi(j))
        end do
        end do
      case('L') ! Lobatto
        do k = 0, po
        do j = 1, ni
          this % A(j,k) = LobattoPolynomial(k, xo, xi(j))
        end do
        end do
      case('N') ! Nodal with arbitrary spacing
        do k = 0, po
        do j = 1, ni
          this % A(j,k) = LagrangePolynomial(k, xo, xi(j))
        end do
        end do
      case default
        call Error('Init_Points', 'Invalid nodal basis', &
                   'Embedded_Interpolation_Operator__1D')
      end select

    end associate

  end subroutine Init_Points

  !-----------------------------------------------------------------------------
  !> Initialization for given 1D standard operators and interpolation points

  subroutine Init_StdOps(this, eop, xi)
    class(EmbeddedInterpolationOperator_1D), intent(inout) :: this
    class(StandardElementOperators_1D), intent(in) :: eop  !< standard operators
    real(RNP), intent(in) :: xi(1:) !< interpolation points in [-1,1]

    call Init_Points(this, eop%nodes, eop%x, xi)

  end subroutine Init_StdOps

  !=============================================================================

end module Embedded_Interpolation_Operator__1D
