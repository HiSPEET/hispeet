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

!> summary:  Fine-to-coarse hp-projection operators
!> author:   Joerg Stiller
!> date:     2023/08/18
!===============================================================================

module HP__Coarsening_Operator__1D
  use Kind_Parameters
  use Constants
  use Gauss_Jacobi
  use Lagrange_Interpolation
  use Execution_Control
  use Standard_Element_Operators__1D

  implicit none
  private

  !-----------------------------------------------------------------------------
  !> hp-coarsening operator
  !>
  !> Provides 1D operators to project data from a fine mesh to a coarser one.
  !> Pure p-coarsening as well as hp-coarsening are supported. In the latter
  !> case, a regular 1:2 refinement of the coarse element is assumed.
  !>
  !> For p-coarsening specify `mode = 1` and for hp-coarsening `mode = 2.
  !> To obtain the identity operator pass `mode = 0`.
  !>
  !> The fine and coarse variables are given as nodal values according to one
  !> of the following Lagrange bases:
  !>   - `E`  equidistant,
  !>   - `G`  Gauss,
  !>   - `L`  Lobatto (default),
  !>   - `RL` Radau left,
  !>   - `RR` Radau right.
  !>
  !> Two different projection methods are provided:
  !>   - `P`  L² projection,
  !>   - `I`  interpolation using the fine basis.
  !>
  !> In the case of hp-coarsening, the treatment of discontinuities between the
  !> two fine elements is controlled by the `smooth` option. The following
  !> choices are available:
  !>   - `0`  no discontinuity handling (default)
  !>   - `1`  jump removal using antisymmetric linear correction
  !>   - `2`  jump removal by averaging interface coefficients
  !>
  !> Option `2` requires a boundary-interior decomposition and is therefore
  !> restricted to equidistant (`E`) or Lobatto (`L`) bases.
  !>
  !> If desired, use the option `filter` to enable exponential filtering.

  type, public :: HP_CoarseningOperator_1D
    character(len=2) :: nodes          !< nodal basis type
    integer          :: po_f           !< polynomial order of fine mesh
    integer          :: po_c           !< polynomial order of coarse mesh
    integer          :: mode           !< coarsening mode
    character        :: method         !< projection method
    real(RNP), allocatable :: A(:,:,:) !< interpolation operator(s)
  end type HP_CoarseningOperator_1D

  ! constructor interface
  interface HP_CoarseningOperator_1D
    procedure New_HP_CoarseningOperator
  end interface

  !-----------------------------------------------------------------------------
  !> Options for initializing the hp-coarsening operator

  type, public :: HP_CoarseningOptions_1D
    character(len=2) :: nodes  = 'L' !< node set {'E','G','L','RL','RR'}
    integer          :: po_f   = -1  !< polynomial order of fine mesh
    integer          :: po_c   = -1  !< polynomial order of coarse mesh
    integer          :: mode   = -1  !< coarsening mode {0,1,2}
    character        :: method = 'I' !< L² projection 'P' or interpolation 'I'
    integer          :: smooth =  0  !< discontinuity handling {0,1,2}
    integer          :: filter =  0  !< exponential filter order, if > 0
  end type HP_CoarseningOptions_1D

contains

  !-----------------------------------------------------------------------------
  !> hp-coarsening operator constructor

  function New_HP_CoarseningOperator(opt) result(this)
    class(HP_CoarseningOptions_1D), intent(in) :: opt
    type(HP_CoarseningOperator_1D) :: this

    call Init_HP_CoarseningOperator(this, opt)

  end function New_HP_CoarseningOperator

  !-----------------------------------------------------------------------------
  !> Build hp-coarsening operator

  subroutine Init_HP_CoarseningOperator(this, opt)
    class(HP_CoarseningOperator_1D), intent(inout) :: this
    class(HP_CoarseningOptions_1D), intent(in) :: opt

    type(StandardElementOperators_1D) :: sop_c

    real(RNP), allocatable :: x_f(:), x_c(:), x_q(:), w_q(:)
    integer :: po_c, po_f, po_q, smooth

    ! initialization ...........................................................

    this % po_f   = opt % po_f
    this % po_c   = opt % po_c
    this % mode   = opt % mode
    this % nodes  = opt % nodes
    this % method = opt % method

    select case(this % nodes)
    case('E','L')
      smooth = max(min(opt % smooth, 2), 0)
    case default
      smooth = max(min(opt % smooth, 1), 0)
    end select

    if (allocated(this % A)) deallocate(this % A)

    ! shorthands
    po_f = this % po_f
    po_c = this % po_c
    po_q = max(po_c, po_f)

    select case(this % mode)
    case(0)
      return
    case(1)
      allocate(this % A(0:po_c,0:po_f,1), source = ZERO)
    case(2)
      allocate(this % A(0:po_c,0:po_f,2), source = ZERO)
    end select

    ! coarse element operators
    sop_c = StandardElementOperators_1D(po_c, this % nodes)

    ! collocation points .......................................................

    allocate(x_c(0:po_c), source = sop_c % x)

    select case(this % nodes)
    case('E')
      block
        integer :: i
        allocate(x_f(0:po_f), source = [(i * TWO/po_f - ONE, i = 0, po_f)])
      end block
    case('G')
      allocate(x_f(0:po_f), source = GaussPoints(po_f))
    case('L')
      allocate(x_f(0:po_f), source = LobattoPoints(po_f))
    case('RL')
      allocate(x_f(0:po_f), source = RadauPoints(po_f, right = .false.))
    case('RR')
      allocate(x_f(0:po_f), source = RadauPoints(po_f, right = .true.))
    end select

    ! quadrature points and weights ............................................

    if (this % mode > 0 .and. this % method == 'P') then

      ! Gauss points and weights for exact integration
      allocate(x_q(0:po_q), source = GaussPoints(po_q))
      allocate(w_q(0:po_q), source = GaussWeights(x_q))

    end if

    ! projection operators .....................................................

    select case(this % mode)

    case(1)

      ! 1:1 projection
      select case(this % method)
      case('P')
        call Build_L2_Projection_Operator_1
      case('I')
        call Build_Interpolation_Operator_1
      end select

    case(2)

      ! 2:1 projection
      select case(this % method)
      case('P')
        call Build_L2_Projection_Operator_2
      case('I')
        call Build_Interpolation_Operator_2
      end select
      call Add_Smoothing

    end select

    if (opt % filter > 0) then
      block
        real(RNP) :: Af(0:po_c,0:po_c), pf
        integer :: k

        pf = real(opt%filter, RNP)
        select case(this % method)
        case('P')
          call sop_c % Get_ErfcLogFilter(pf, Af, modes = 'L', left_half=.true.)
        case('I')
          call sop_c % Get_ErfcLogFilter(pf, Af, modes = 'B', left_half=.true.)
        end select

        do k = 1, this % mode
          this % A(:,:,k) = matmul(Af, this % A(:,:,k))
        end do

      end block
    end if

  contains

    !---------------------------------------------------------------------------
    !> Build 1:1 L² projection

    subroutine Build_L2_Projection_Operator_1

      real(RNP), allocatable :: l_c(:,:), l_f(:,:)
      real(RNP), allocatable :: MA(:,:), M_inv(:,:)
      real(RNP), allocatable :: VL(:,:), ML_inv(:)
      integer :: i, j, k, q

      ! coarse basis functions at quadrature points ............................

      allocate(l_c(0:po_q, 0:po_c))
      do i = 0, po_c
      do q = 0, po_q
        select case(this % nodes)
        case('E')
          l_c(q,i) = LagrangePolynomial(i, x_c, x_q(q))
        case('G')
          l_c(q,i) = GaussPolynomial(i, x_c, x_q(q))
        case('L')
          l_c(q,i) = LobattoPolynomial(i, x_c, x_q(q))
        case('RL','RR')
          l_c(q,i) = RadauPolynomial(i, x_c, x_q(q))
        end select
      end do
      end do

      ! fine basis functions at quadrature points ..............................

      allocate(l_f(0:po_q, 0:po_f))
      do j = 0, po_f
      do q = 0, po_q
        select case(this % nodes)
        case('E')
          l_f(q,j) = LagrangePolynomial(j, x_f, x_q(q))
        case('G')
          l_f(q,j) = GaussPolynomial(j, x_f, x_q(q))
        case('L')
          l_f(q,j) = LobattoPolynomial(j, x_f, x_q(q))
        case('RL','RR')
          l_f(q,j) = RadauPolynomial(j, x_f, x_q(q))
        case default
          l_f(q,j) = LagrangePolynomial(j, x_f, x_q(q))
        end select
      end do
      end do

      ! mass weighted projection operator ......................................

      allocate(MA(0:po_c, 0:po_f), source = ZERO)

      do j = 0, po_f
      do i = 0, po_c
        MA(i,j) = sum(w_q * l_c(:,i) * l_f(:,j))
      end do
      end do

      ! L² projection operator .................................................

      ! Vandermode matrix to operator basis functions
      allocate(VL(0:po_c,0:po_c))
      call sop_c % Get_Legendre_VDM(VL)

      ! inverse diagonal Legendre mass matrix
      allocate(ML_inv(0:po_c))
      do k = 0, po_c
        ML_inv(k) = HALF + k
      end do

      ! inverse mass matrix
      allocate(M_inv(0:po_c,0:po_c), source = ZERO)
      do j = 0, po_c
      do i = 0, po_c
        M_inv(i,j) = sum(VL(i,:) * ML_inv * VL(j,:))
      end do
      end do

      ! projection operator
      this % A(:,:,1) = matmul(M_inv, MA)

    end subroutine Build_L2_Projection_Operator_1

    !---------------------------------------------------------------------------
    !> Build 1:1 projection based on embedded interpolation

    subroutine Build_Interpolation_Operator_1

      integer :: i, j

      select case(this % nodes)

      case('E') ! Nodal with equidistant spacing
        do i = 0, po_c
        do j = 0, po_f
          this % A(i,j,1) = LagrangePolynomial(j, x_f, x_c(i))
        end do
        end do

      case('G') ! Gauss
        do i = 0, po_c
        do j = 0, po_f
          this % A(i,j,1) = GaussPolynomial(j, x_f, x_c(i))
        end do
        end do

      case('L') ! Lobatto
        do i = 0, po_c
        do j = 0, po_f
          this % A(i,j,1) = LobattoPolynomial(j, x_f, x_c(i))
        end do
        end do

      case('RL','RR') ! Radau
        do i = 0, po_c
        do j = 0, po_f
          this % A(i,j,1) = RadauPolynomial(j, x_f, x_c(i))
        end do
        end do

      end select

    end subroutine Build_Interpolation_Operator_1

    !---------------------------------------------------------------------------
    !> Build 2:1 L² projection

    subroutine Build_L2_Projection_Operator_2

      real(RNP), allocatable :: l_c(:,:,:), l_f(:,:)
      real(RNP), allocatable :: MA(:,:,:), M_inv(:,:)
      real(RNP), allocatable :: VL(:,:), ML_inv(:)
      real(RNP) :: x_q1, x_q2
      integer :: i, j, k, q

      ! coarse basis functions at quadrature points in fine elements 1 and 2 ...

      allocate(l_c(0:po_q, 0:po_c, 2))
      do i = 0, po_c
      do q = 0, po_q
        x_q1 = HALF*(x_q(q) - 1)
        x_q2 = HALF*(x_q(q) + 1)
        select case(this % nodes)
        case('E')
          l_c(q,i,1) = LagrangePolynomial(i, x_c, x_q1)
          l_c(q,i,2) = LagrangePolynomial(i, x_c, x_q2)
        case('G')
          l_c(q,i,1) = GaussPolynomial(i, x_c, x_q1)
          l_c(q,i,2) = GaussPolynomial(i, x_c, x_q2)
        case('L')
          l_c(q,i,1) = LobattoPolynomial(i, x_c, x_q1)
          l_c(q,i,2) = LobattoPolynomial(i, x_c, x_q2)
        case('RL','RR')
          l_c(q,i,1) = RadauPolynomial(i, x_c, x_q1)
          l_c(q,i,2) = RadauPolynomial(i, x_c, x_q2)
        end select
      end do
      end do

      ! fine basis functions at quadrature points ..............................

      allocate(l_f(0:po_q, 0:po_f))
      do j = 0, po_f
      do q = 0, po_q
        select case(this % nodes)
        case('E')
          l_f(q,j) = LagrangePolynomial(j, x_f, x_q(q))
        case('G')
          l_f(q,j) = GaussPolynomial(j, x_f, x_q(q))
        case('L')
          l_f(q,j) = LobattoPolynomial(j, x_f, x_q(q))
        case('RL','RR')
          l_f(q,j) = RadauPolynomial(j, x_f, x_q(q))
        case default
          l_f(q,j) = LagrangePolynomial(j, x_f, x_q(q))
        end select
      end do
      end do

      ! mass weighted projection operator ......................................

      allocate(MA(0:po_c, 0:po_f, 2), source = ZERO)

      do j = 0, po_f
      do i = 0, po_c
        MA(i,j,1) = HALF * sum(w_q * l_c(:,i,1) * l_f(:,j))
        MA(i,j,2) = HALF * sum(w_q * l_c(:,i,2) * l_f(:,j))
      end do
      end do

      ! L² projection operator .................................................

      ! Vandermode matrix to operator basis functions
      allocate(VL(0:po_c,0:po_c))
      call sop_c % Get_Legendre_VDM(VL)

      ! inverse diagonal Legendre mass matrix
      allocate(ML_inv(0:po_c))
      do k = 0, po_c
        ML_inv(k) = HALF + k
      end do

      ! inverse mass matrix
      allocate(M_inv(0:po_c,0:po_c), source = ZERO)
      do j = 0, po_c
      do i = 0, po_c
        M_inv(i,j) = sum(VL(i,:) * ML_inv * VL(j,:))
      end do
      end do

      ! projection operator
      this % A(:,:,1) = matmul(M_inv, MA(:,:,1))
      this % A(:,:,2) = matmul(M_inv, MA(:,:,2))


    end subroutine Build_L2_Projection_Operator_2

    !---------------------------------------------------------------------------
    !> Build 2:1 projection based on embedded interpolation

    subroutine Build_Interpolation_Operator_2

      real(RNP) :: x_c2f(2)
      integer   :: o2_f, ph_c
      integer   :: i1, i2, j

      ph_c = po_c / 2
      o2_f = po_c - ph_c

      do i1 = 0, ph_c

        i2 = i1 + o2_f

        x_c2f(1) = 2*x_c(i1) + ONE ! coarse point mapped to left  element
        x_c2f(2) = 2*x_c(i2) - ONE ! coarse point mapped to right element

        select case(this % nodes)
        case('E') ! Nodal with equidistant spacing
          do j = 0, po_f
            this % A(i1,j,1) = LagrangePolynomial(j, x_f, x_c2f(1))
            this % A(i2,j,2) = LagrangePolynomial(j, x_f, x_c2f(2))
          end do
        case('G') ! Gauss
          do j = 0, po_f
            this % A(i1,j,1) = GaussPolynomial(j, x_f, x_c2f(1))
            this % A(i2,j,2) = GaussPolynomial(j, x_f, x_c2f(2))
          end do
        case('L') ! Lobatto
          do j = 0, po_f
            this % A(i1,j,1) = LobattoPolynomial(j, x_f, x_c2f(1))
            this % A(i2,j,2) = LobattoPolynomial(j, x_f, x_c2f(2))
          end do
        case('RL','RR') ! Radau: take care of asymmetry
          if (x_c2f(1) <= ONE) then
            do j = 0, po_f
              this % A(i1,j,1) = RadauPolynomial(j, x_f, x_c2f(1))
            end do
          end if
          if (x_c2f(2) >= -ONE) then
            do j = 0, po_f
              this % A(i2,j,2) = RadauPolynomial(j, x_f, x_c2f(2))
            end do
          end if
        end select

      end do

      select case(this % nodes)
      case('E','G','L')
        ! averaging at interface
        if (mod(po_c,2) == 0) then
          do j = 0, po_f
            this % A(ph_c,j,1) = HALF * this % A(ph_c,j,1)
            this % A(ph_c,j,2) = HALF * this % A(ph_c,j,2)
          end do
        end if
      end select

    end subroutine Build_Interpolation_Operator_2

    !---------------------------------------------------------------------------
    !> Build combined smoothing + projection operator

    subroutine Add_Smoothing

      real(RNP), allocatable :: B(:,:), L(:), R(:)
      real(RNP) :: AB
      integer :: i, j

      if (smooth == 0) return

      ! smothing operator ......................................................

      allocate(B(0:po_f,2), source = ZERO)

      select case(smooth)
      case(1)
        ! linear blending
        B(:,1) = -(ONE + x_f) / 4
        B(:,2) =  (ONE - x_f) / 4
      case(2)
        ! interface coefficient averaging
        B(po_f,1) = -HALF
        B( 0  ,2) =  HALF
      end select

      ! trace operators ........................................................

      allocate(L(0:po_f), R(0:po_f))

      do i = 0, po_f
        select case(this % nodes)
        case('E') ! Nodal with equidistant spacing
          L(i) = LagrangePolynomial(i, x_f, -ONE)
          R(i) = LagrangePolynomial(i, x_f,  ONE)
        case('G') ! Gauss
          L(i) = GaussPolynomial(i, x_f, -ONE)
          R(i) = GaussPolynomial(i, x_f,  ONE)
        case('L') ! Lobatto
          L(i) = LobattoPolynomial(i, x_f, -ONE)
          R(i) = LobattoPolynomial(i, x_f,  ONE)
        case('RL','RR') ! Radau: take care of asymmetry
          L(i) = RadauPolynomial(i, x_f, -ONE)
          R(i) = RadauPolynomial(i, x_f,  ONE)
        end select
      end do

      ! combined operator ......................................................

      associate(A => this % A)
        do i = 0, po_c
          AB =  sum(A(i,:,1) * B(:,1) + A(i,:,2) * B(:,2))
          do j = 0, po_f
            A(i,j,1) = A(i,j,1) + AB * R(j)
            A(i,j,2) = A(i,j,2) - AB * L(j)
          end do
        end do
      end associate

    end subroutine Add_Smoothing

    !---------------------------------------------------------------------------

  end subroutine Init_HP_CoarseningOperator

  !=============================================================================

end module HP__Coarsening_Operator__1D
