!> summary:  Coarse-to-fine hp-interpolation operators
!> author:   Joerg Stiller
!> date:     2023/08/19
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module HP__Refinement_Operator__1D
  use Kind_Parameters
  use Constants
  use Gauss_Jacobi
  use Lagrange_Interpolation
  use Execution_Control

  implicit none
  private

  !-----------------------------------------------------------------------------
  !> hp-refinement operator
  !>
  !> Provides 1D operators for interpolating data from a coarse mesh to a finer
  !> a coarser one. Pure p-refinement as well as hp-refinement are supported. In
  !> the latter case, a regular 2:1 refinement of the coarse element is assumed.
  !> When initializing set `mode` to 1 for p-refinement and 2 for hp-refinement.
  !> `mode = 0` indicates the identity, i.e., no interpolation is performed and,
  !> hence, no interpolation operator is provided.

  type, public :: HP_RefinementOperator_1D
    character(len=2) :: nodes  = 'L'   !< nodal basis type
    integer          :: po_c   = -1    !< polynomial order of coarse mesh
    integer          :: po_f   = -1    !< polynomial order of fine mesh
    integer          :: mode   = -1    !< refinement mode
    real(RNP), allocatable :: A(:,:,:) !< interpolation operator(s)
  end type HP_RefinementOperator_1D

  ! constructor interface
  interface HP_RefinementOperator_1D
    procedure New_HP_RefinementOperator
  end interface

  !-----------------------------------------------------------------------------
  !> Options for initializing the hp-refinement operator

  type, public :: HP_RefinementOptions_1D
    character(len=2) :: nodes  = 'L' !< node set {'E','G','L','RL','RR'}
    integer          :: po_c   = -1  !< polynomial order of coarse mesh
    integer          :: po_f   = -1  !< polynomial order of fine mesh
    integer          :: mode   = -1  !< refinement mode {0,1,2}
  end type HP_RefinementOptions_1D

contains

  !-----------------------------------------------------------------------------
  !> hp-refinement operator constructor

  function New_HP_RefinementOperator(opt) result(this)
    class(HP_RefinementOptions_1D), intent(in) :: opt
    type(HP_RefinementOperator_1D) :: this

    call Init_HP_RefinementOperator(this, opt)

  end function New_HP_RefinementOperator

  !-----------------------------------------------------------------------------
  !> Build hp-refinement

  subroutine Init_HP_RefinementOperator(this, opt)
    class(HP_RefinementOperator_1D), intent(inout) :: this
    class(HP_RefinementOptions_1D), intent(in) :: opt

    real(RNP), allocatable :: x_f(:), x_c(:)
    integer :: po_c, po_f
    integer :: i, j

    ! initialization ...........................................................

    this % nodes = opt % nodes
    this % po_c  = opt % po_c
    this % po_f  = opt % po_f
    this % mode  = opt % mode

    ! shorthands
    po_c = this % po_c
    po_f = this % po_f

    if (allocated(this % A)) deallocate(this % A)

    allocate(this % A(0:po_f, 0:po_c, this%mode))

    if (this % mode == 0) return

    ! collocation points .......................................................

    select case(this % nodes)
    case('E')
      allocate(x_c(0:po_c), source = [(i * TWO/po_c - ONE, i = 0, po_c)])
      allocate(x_f(0:po_f), source = [(i * TWO/po_f - ONE, i = 0, po_f)])
    case('G')
      allocate(x_c(0:po_c), source = GaussPoints(po_c))
      allocate(x_f(0:po_f), source = GaussPoints(po_f))
    case('L')
      allocate(x_c(0:po_c), source = LobattoPoints(po_c))
      allocate(x_f(0:po_f), source = LobattoPoints(po_f))
    case('RL')
      allocate(x_c(0:po_c), source = RadauPoints(po_c, right = .false.))
      allocate(x_f(0:po_f), source = RadauPoints(po_f, right = .false.))
    case('RR')
      allocate(x_c(0:po_c), source = RadauPoints(po_c, right = .true.))
      allocate(x_f(0:po_f), source = RadauPoints(po_f, right = .true.))
    end select

    ! interpolation matrices ...................................................

    select case(this % mode)

    case(1)
      ! 1:1 interpolation
      select case(this % nodes)
      case('E') ! Nodal with equidistant spacing
        do i = 0, po_f
        do j = 0, po_c
          this % A(i,j,1) = LagrangePolynomial(j, x_c, x_f(i))
        end do
        end do
      case('G') ! Gauss
        do i = 0, po_f
        do j = 0, po_c
          this % A(i,j,1) = GaussPolynomial(j, x_c, x_f(i))
        end do
        end do
      case('L') ! Lobatto
        do i = 0, po_f
        do j = 0, po_c
          this % A(i,j,1) = LobattoPolynomial(j, x_c, x_f(i))
        end do
        end do
      case('RL','RR') ! Radau
        do i = 0, po_f
        do j = 0, po_c
          this % A(i,j,1) = RadauPolynomial(j, x_c, x_f(i))
        end do
        end do
      end select

    case(2)
      ! 1:2 interpolation
      select case(this % nodes)
      case('E') ! Nodal with equidistant spacing
        do i = 0, po_f
        do j = 0, po_c
          this % A(i,j,1) = LagrangePolynomial(j, x_c, HALF * (x_f(i) - ONE))
          this % A(i,j,2) = LagrangePolynomial(j, x_c, HALF * (x_f(i) + ONE))
        end do
        end do
      case('G') ! Gauss
        do i = 0, po_f
        do j = 0, po_c
          this % A(i,j,1) = GaussPolynomial(j, x_c, HALF * (x_f(i) - ONE))
          this % A(i,j,2) = GaussPolynomial(j, x_c, HALF * (x_f(i) + ONE))
        end do
        end do
      case('L') ! Lobatto
        do i = 0, po_f
        do j = 0, po_c
          this % A(i,j,1) = LobattoPolynomial(j, x_c, HALF * (x_f(i) - ONE))
          this % A(i,j,2) = LobattoPolynomial(j, x_c, HALF * (x_f(i) + ONE))
        end do
        end do
      case('RL','RR') ! Radau
        do i = 0, po_f
        do j = 0, po_c
          this % A(i,j,1) = RadauPolynomial(j, x_c, HALF * (x_f(i) - ONE))
          this % A(i,j,2) = RadauPolynomial(j, x_c, HALF * (x_f(i) + ONE))
        end do
        end do
      end select

    end select

  end subroutine Init_HP_RefinementOperator

  !=============================================================================

end module HP__Refinement_Operator__1D
