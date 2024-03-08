!> summary:  Embedded interpolation operator
!> author:   Joerg Stiller
!> date:     2018/03/09
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> Provides one-dimensional interpolation operators for the following `basis`
!> types
!>
!>   - `'G'`   Lagrangian with Gauss points
!>   - `'L'`   Lagrangian with on Lobatto points
!>   - `'RL'`  Lagrangian with left-sided Radau points
!>   - `'RR'`  Lagrangian with right-sided Radau points
!>   - `'N'`   Lagrangian with arbitrary nodes
!>
!===============================================================================

module Embedded_Interpolation__1D
  use Kind_Parameters, only: RNP
  use Execution_Control, only: Error
  use Gauss_Jacobi
  use Lagrange_Interpolation
  use Standard_Operators__1D
  implicit none
  private

  !-----------------------------------------------------------------------------
  !> Embedded interpolation operator

  type, public :: EmbeddedInterpolation_1D

    integer :: no = -1 !< number of original points per direction
    integer :: ni = -1 !< number of interpolated points per direction
    real(RNP), allocatable :: A(:,:)  !< 1D interpolation operator

  contains

    generic :: Init_EmbeddedInterpolation_1D => Init_Points, Init_StdOps
    procedure, private :: Init_Points
    procedure, private :: Init_StdOps

  end type EmbeddedInterpolation_1D

  ! constructor interface
  interface EmbeddedInterpolation_1D
    module procedure New_Points
    module procedure New_StdOps
  end interface

contains

  !=============================================================================
  ! Constructors

  !-----------------------------------------------------------------------------
  !> New interpolation operator from basis type and points

  type(EmbeddedInterpolation_1D) function New_Points(basis, xo, xi) result(this)
    character(*), intent(in) :: basis  !< basis type
    real(RNP),    intent(in) :: xo(0:) !< basis points
    real(RNP),    intent(in) :: xi(1:) !< interpolation points in [-1,1]

    call Init_Points(this, basis, xo, xi)

  end function New_Points

  !-----------------------------------------------------------------------------
  !> New interpolation operator from standard operators and interpolation points

  type(EmbeddedInterpolation_1D) function New_StdOps(eop, xi) result(this)
    class(StandardOperators_1D), intent(in) :: eop !< standard operators
    real(RNP), intent(in) :: xi(1:) !< interpolation points in [-1,1]

    call Init_StdOps(this, eop, xi)

  end function New_StdOps

  !=============================================================================
  ! Type-bound procedures

  !-----------------------------------------------------------------------------
  !> Initialization for given nodes and interpolation points

  subroutine Init_Points(this, basis, xo, xi)
    class(EmbeddedInterpolation_1D), intent(inout) :: this
    character(*), intent(in) :: basis  !< basis type
    real(RNP),    intent(in) :: xo(0:) !< basis points
    real(RNP),    intent(in) :: xi(1:) !< interpolation points in [-1,1]

    integer :: j, k, po

    if (allocated(this % A)) deallocate(this % A)

    associate(no => this % no, ni => this % ni)

      po = ubound(xo, 1)
      no = size(xo)
      ni = size(xi)

      allocate(this % A(ni,0:po))

      select case(basis)
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
        call Error('Init_Points', 'Invalid basis', 'Embedded_Interpolation__1D')
      end select

    end associate

  end subroutine Init_Points

  !-----------------------------------------------------------------------------
  !> Initialization for given 1D standard operators and interpolation points

  subroutine Init_StdOps(this, eop, xi)
    class(EmbeddedInterpolation_1D), intent(inout) :: this
    class(StandardOperators_1D), intent(in) :: eop  !< standard operators
    real(RNP), intent(in) :: xi(1:) !< interpolation points in [-1,1]

    call Init_Points(this, eop%basis, eop%x, xi)

  end subroutine Init_StdOps

  !=============================================================================

end module Embedded_Interpolation__1D
