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
  use TPO__AAA_3d
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
    module procedure New_SX
  end interface

contains

!===============================================================================
! Constructor

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

  call TPO_AAA(this%A, uo, ui)

end subroutine Interpolate_S

!-------------------------------------------------------------------------------
!> Interpolation of array variables

subroutine Interpolate_A(this, uo, ui)
  class(InterpolationOperator3D), intent(in) :: this
  real(RNP), intent(in)  :: uo(:,:,:,:,:) !< original mesh variable
  real(RNP), intent(out) :: ui(:,:,:,:,:) !< interpolated mesh variable

  integer :: c

  do c = 1, size(uo,5)
    call TPO_AAA(this%A, uo(:,:,:,:,c), ui(:,:,:,:,c))
  end do

end subroutine Interpolate_A

!===============================================================================

end module Embedded_Interpolation_3D
