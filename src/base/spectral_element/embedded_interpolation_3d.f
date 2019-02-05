!> summary:  Embedded interpolation of 3D mesh variables
!> author:   Joerg Stiller
!> date:     2018/11/02
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Embedded interpolation of 3D mesh variables
!===============================================================================

module Embedded_Interpolation_3D
  use Kind_Parameters, only: RNP
  use Standard_Operators_1D
  use Embedded_Interpolation
  use TPO_AAA
  implicit none
  private

  !-----------------------------------------------------------------------------
  !> Embedded interpolation of 3D mesh variables

  type, extends(InterpolationOperator), public :: InterpolationOperator3D

  contains

    generic :: Apply  =>  Interpolate_S, Interpolate_A
    procedure, private :: Interpolate_S
    procedure, private :: Interpolate_A

  end type InterpolationOperator3D

  ! constructor interface
  interface InterpolationOperator3D
    module procedure New_PP
    module procedure New_PX
    module procedure New_SP
    module procedure New_SX
  end interface

contains

!===============================================================================
! Constructors

!-------------------------------------------------------------------------------
!> New InterpolationOperator3D from orders of source and interpolant

type(InterpolationOperator3D) function New_PP(po, pi) result(this)
  integer, intent(in) :: po  !< polynomial order of the original
  integer, intent(in) :: pi  !< polynomial order of the interpolant

  call this % Init_InterpolationOperator(po, pi)

end function New_PP

!-------------------------------------------------------------------------------
!> New InterpolationOperator from order of source and interpolation points

type(InterpolationOperator3D) function New_PX(po, xi) result(this)
  integer,   intent(in) :: po     !< polynomial order of the original
  real(RNP), intent(in) :: xi(0:) !< interpolation points in [-1,1]

  call this % Init_InterpolationOperator(po, xi)

end function New_PX

!-------------------------------------------------------------------------------
!> New InterpolationOperator from 1D standard operators and order of interpolant

type(InterpolationOperator3D) function New_SP(eop, pi) result(this)
  class(StandardOperators1D), intent(in) :: eop !< standard operators
  integer,                    intent(in) :: pi  !< order of interpolant

  call this % Init_InterpolationOperator(eop, pi)

end function New_SP

!-------------------------------------------------------------------------------
!> New InterpolationOperator from 1D standard operators and interpolation points

type(InterpolationOperator3D) function New_SX(eop, xi) result(this)
  class(StandardOperators1D), intent(in) :: eop    !< standard operators
  real(RNP),                  intent(in) :: xi(0:) !< points in [-1,1]

  call this % Init_InterpolationOperator(eop, xi)

end function New_SX

!===============================================================================
! Type-bound procedures

!-------------------------------------------------------------------------------
!> Interpolation of scalar variables

subroutine Interpolate_S(this, uo, ui)
  class(InterpolationOperator3D), intent(in) :: this
  real(RNP), intent(in)  :: uo(:,:,:,:) !< original mesh variable
  real(RNP), intent(out) :: ui(:,:,:,:) !< interpolated mesh variable

  call TPO_AAA_Eval(this%ni, this%no, size(uo,4), this%iop, uo, ui)

end subroutine Interpolate_S

!-------------------------------------------------------------------------------
!> Interpolation of array variables

subroutine Interpolate_A(this, uo, ui)
  class(InterpolationOperator3D), intent(in) :: this
  real(RNP), intent(in)  :: uo(:,:,:,:,:) !< original mesh variable
  real(RNP), intent(out) :: ui(:,:,:,:,:) !< interpolated mesh variable

  call TPO_AAA_Eval(this%ni, this%no, size(uo,4)*size(uo,5), this%iop, uo, ui)

end subroutine Interpolate_A

!===============================================================================

end module Embedded_Interpolation_3D
