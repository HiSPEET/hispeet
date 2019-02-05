!> summary:  Embedded interpolation operator
!> author:   Joerg Stiller
!> date:     2018/03/09
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Embedded interpolation operator
!===============================================================================

module Embedded_Interpolation
  use Kind_Parameters, only: RNP
  use Gauss_Jacobi
  use Standard_Operators_1D
  implicit none
  private

  !-----------------------------------------------------------------------------
  !> Embedded interpolation operator

  type, public :: InterpolationOperator

    integer :: no = -1 !< number of original points per direction
    integer :: ni = -1 !< number of interpolated points per direction
    real(RNP), allocatable :: iop(:,:)  !< 1D interpolation operator

  contains

    generic :: Init_InterpolationOperator => Init_PP, Init_PX, Init_SP, Init_SX
    procedure, private :: Init_PP
    procedure, private :: Init_PX
    procedure, private :: Init_SP
    procedure, private :: Init_SX

  end type InterpolationOperator

  ! constructor interface
  interface InterpolationOperator
    module procedure New_PP
    module procedure New_PX
    module procedure New_SP
    module procedure New_SX
  end interface

contains

!===============================================================================
! Constructors

!-------------------------------------------------------------------------------
!> New InterpolationOperator from orders of source and interpolant

type(InterpolationOperator) function New_PP(po, pi) result(this)
  integer, intent(in) :: po  !< polynomial order of the original
  integer, intent(in) :: pi  !< polynomial order of the interpolant

  call Init_PP(this, po, pi)

end function New_PP

!-------------------------------------------------------------------------------
!> New InterpolationOperator from order of source and interpolation points

type(InterpolationOperator) function New_PX(po, xi) result(this)
  integer,   intent(in) :: po     !< polynomial order of the original
  real(RNP), intent(in) :: xi(0:) !< interpolation points in [-1,1]

  call Init_PX(this, po, xi)

end function New_PX

!-------------------------------------------------------------------------------
!> New InterpolationOperator from 1D standard operators and order of interpolant

type(InterpolationOperator) function New_SP(eop, pi) result(this)
  class(StandardOperators1D), intent(in) :: eop !< standard operators
  integer,                    intent(in) :: pi  !< order of interpolant

  call Init_SP(this, eop, pi)

end function New_SP

!-------------------------------------------------------------------------------
!> New InterpolationOperator from 1D standard operators and interpolation points

type(InterpolationOperator) function New_SX(eop, xi) result(this)
  class(StandardOperators1D), intent(in) :: eop    !< standard operators
  real(RNP),                  intent(in) :: xi(0:) !< points in [-1,1]

  call Init_SX(this, eop, xi)

end function New_SX

!===============================================================================
! Type-bound procedures

!-------------------------------------------------------------------------------
!> Initialization for given orders of source and interpolant

subroutine Init_PP(this, po, pi)
  class(InterpolationOperator), intent(inout) :: this
  integer, intent(in) :: po  !< polynomial order of the original
  integer, intent(in) :: pi  !< polynomial order of the interpolant

  real(RNP) :: xo(0:po), xi(0:pi)
  integer   :: j, k

  if (this % no > 0) call Delete_InterpolationOperator(this)

  xo = GLL_Points(po)
  xi = GLL_Points(pi)

  this % no = size(xo)
  this % ni = size(xi)

  allocate(this % iop(0:pi,0:po))

  do k = 0, po
  do j = 0, pi
    this % iop(j,k) = GLL_Polynomial(k, xo, xi(j))
  end do
  end do

end subroutine Init_PP

!-------------------------------------------------------------------------------
!> Initialization for given order of source and interpolation points

subroutine Init_PX(this, po, xi)
  class(InterpolationOperator), intent(inout) :: this
  integer,   intent(in) :: po     !< polynomial order of the original
  real(RNP), intent(in) :: xi(0:) !< interpolation points in [-1,1]

  real(RNP) :: xo(0:po)
  integer   :: j, k, pi

  if (this % no > 0) call Delete_InterpolationOperator(this)

  xo = GLL_Points(po)
  pi = ubound(xi,1)

  this % no = size(xo)
  this % ni = size(xi)

  allocate(this % iop(0:pi,0:po))

  do k = 0, po
  do j = 0, pi
    this % iop(j,k) = GLL_Polynomial(k, xo, xi(j))
  end do
  end do

end subroutine Init_PX

!-------------------------------------------------------------------------------
!> Initialization for given 1D standard operators and order of interpolant

subroutine Init_SP(this, eop, pi)
  class(InterpolationOperator), intent(inout) :: this
  class(StandardOperators1D),   intent(in)    :: eop  !< standard operators
  integer,                      intent(in)    :: pi   !< order of interpolant

  real(RNP) :: xi(0:pi)
  integer   :: j, k

  if (this % no > 0) call Delete_InterpolationOperator(this)

  associate(po => eop%po, xo => eop%x)

    xi = GLL_Points(pi)

    this % no = size(xo)
    this % ni = size(xi)

    allocate(this % iop(0:pi,0:po))

    do k = 0, po
    do j = 0, pi
      this % iop(j,k) = GLL_Polynomial(k, xo, xi(j))
    end do
    end do

  end associate

end subroutine Init_SP

!-------------------------------------------------------------------------------
!> Initialization for given 1D standard operators and interpolation points

subroutine Init_SX(this, eop, xi)
  class(InterpolationOperator), intent(inout) :: this
  class(StandardOperators1D),   intent(in)    :: eop    !< standard operators
  real(RNP),                    intent(in)    :: xi(0:) !< points in [-1,1]

  integer   :: j, k, pi

  if (this % no > 0) call Delete_InterpolationOperator(this)

  associate(po => eop%po, xo => eop%x)

    pi = ubound(xi,1)

    this % no = size(xo)
    this % ni = size(xi)

    allocate(this % iop(0:pi,0:po))

    do k = 0, po
    do j = 0, pi
      this % iop(j,k) = GLL_Polynomial(k, xo, xi(j))
    end do
    end do

  end associate

end subroutine Init_SX

!-------------------------------------------------------------------------------
!> Delete InterpolationOperator object

subroutine Delete_InterpolationOperator(this)
  class(InterpolationOperator), intent(inout) :: this

  if (allocated(this % iop)) deallocate(this % iop)

end subroutine Delete_InterpolationOperator

!===============================================================================

end module Embedded_Interpolation
