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

!> summary:  Projection operator
!> author:   Joerg Stiller
!> date:     2019/02/22, 2024/12/09, 2026/02/20
!===============================================================================

module Projection_Operator__1D
  use Kind_Parameters, only: RNP
  use Constants,       only: ZERO, HALF
  use Gauss_Jacobi
  use Lagrange_Interpolation
  use Standard_Element_Operators__1D
  implicit none
  private

  public :: ProjectionOperator_1D

  !-----------------------------------------------------------------------------
  !> Element-based projection operators
  !>
  !> Provides the one-dimensional operators `A` and `MA` for projecting a
  !> function `f` given at the collocation points `x(0:po_f) ∈ [-1,1]` to
  !> the polynomial basis functions `{𝓁ᵢ}` of degree `po_b` as defined in
  !> the standard operators `eop`.
  !>
  !> The operator `MA` yields the mass-weighted projection
  !>
  !>      (MA f)ᵢ ≈ ∫ 𝓁ᵢ(ξ) f(ξ) dξ
  !>
  !> and `A` the L² projection
  !>
  !>      A f = M⁻¹ MA f.
  !>
  !> Both projections refer to the standard element [-1,1] and must be
  !> scaled properly when applied to physical elements.

  type ProjectionOperator_1D
    integer :: po_b = -1 !< polynomial order of the basis {𝓁ᵢ}
    integer :: po_f = -1 !< polynomial order of the projected function
    real(RNP), allocatable :: A (:,:) !< L2 projection operator (0:po_b,0:po_f)
    real(RNP), allocatable :: MA(:,:) !< mass-weighted operator (0:po_b,0:po_f)
  end type ProjectionOperator_1D

  ! constructor interface
  interface ProjectionOperator_1D
    module procedure New_ProjectionOperator_1D
  end interface

contains

  !-----------------------------------------------------------------------------
  !> New ProjectionOperator_1D from 1D standard operators and interpolation

  function New_ProjectionOperator_1D(eop, x, nodes) result(this)
    class(StandardElementOperators_1D), intent(in) :: eop !< standard operators
    real(RNP),              intent(in) :: x(0:) !< node coordinates
    character(*), optional, intent(in) :: nodes !< node type
    type(ProjectionOperator_1D) :: this

    call Init_ProjectionOperator_1D(this, eop, x, nodes)

  end function New_ProjectionOperator_1D

  !-----------------------------------------------------------------------------
  !> Initialization

  subroutine Init_ProjectionOperator_1D(this, eop, x, nodes)
    class(ProjectionOperator_1D), intent(inout) :: this
    class(StandardElementOperators_1D), intent(in) :: eop !< standard operators
    real(RNP),              intent(in) :: x(0:) !< node coordinates
    character(*), optional, intent(in) :: nodes !< node type

    real(RNP), allocatable :: x_q(:), w_q(:)
    real(RNP), allocatable :: l_b(:,:), l_f(:,:)
    real(RNP), allocatable :: VL(:,:), ML_inv(:), M_inv(:,:)

    character(2) :: nodes_f
    integer :: po_b, po_f, po_q
    integer :: i, j, k, q

    ! initialization ...........................................................

    po_b = eop % po
    po_f = ubound(x, 1)
    po_q = max(po_b, po_f)

    if (present(nodes)) then
      nodes_f = nodes
    else
      nodes_f = ''
    end if

    this % po_b = po_b
    this % po_f = po_f

    allocate(this % A  (0:po_b, 0:po_f), source = ZERO)
    allocate(this % MA (0:po_b, 0:po_f), source = ZERO)

    ! Gauss points and weights for exact integration
    allocate(x_q(0:po_q), source = GaussPoints(po_q))
    allocate(w_q(0:po_q), source = GaussWeights(x_q))

    ! operator basis functions at quadrature points
    allocate(l_b(0:po_q, 0:po_b))
    do i = 0, po_b
    do q = 0, po_q
      select case(eop % nodes)
      case('G')
        l_b(q,i) = GaussPolynomial(i, eop%x, x_q(q))
      case('L')
        l_b(q,i) = LobattoPolynomial(i, eop%x, x_q(q))
      case('RL','RR')
        l_b(q,i) = RadauPolynomial(i, eop%x, x_q(q))
      end select
    end do
    end do

    ! Lagrange cardinal functions to `x` at quadrature points
    allocate(l_f(0:po_q, 0:po_f))
    do j = 0, po_f
    do q = 0, po_q
      select case(nodes_f)
      case('G')
        l_f(q,j) = GaussPolynomial(j, x, x_q(q))
      case('L')
        l_f(q,j) = LobattoPolynomial(j, x, x_q(q))
      case('RL','RR')
        l_f(q,j) = RadauPolynomial(j, x, x_q(q))
      case default
        l_f(q,j) = LagrangePolynomial(j, x, x_q(q))
      end select
    end do
    end do


    associate(A => this % A, MA => this % MA)

      ! mass weighted projection operator ......................................

      do j = 0, po_f
      do i = 0, po_b
        MA(i,j) = sum(w_q * l_b(:,i) * l_f(:,j))
      end do
      end do

      ! L² projection operator .................................................

      ! Vandermode matrix to operator basis functions
      allocate(VL(0:po_b,0:po_b))
      call eop % Get_Legendre_VDM(VL)

      ! inverse diagonal Legendre mass matrix
      allocate(ML_inv(0:po_b))
      do k = 0, po_b
        ML_inv(k) = HALF + k
      end do

      ! inverse mass matrix
      allocate(M_inv(0:po_b,0:po_b), source = ZERO)
      do j = 0, po_b
      do i = 0, po_b
        M_inv(i,j) = sum(VL(i,:) * ML_inv * VL(j,:))
      end do
      end do

      ! projection operator
      A = matmul(M_inv, MA)

    end associate

  end subroutine Init_ProjectionOperator_1D

  !=============================================================================

end module Projection_Operator__1D
