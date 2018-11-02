!> summary:  Embedded interpolation of 3D mesh variables
!> author:   Joerg Stiller
!> date:     2018/11/02
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Embedded interpolation of 3D mesh variables
!===============================================================================

module Embedded_Interpolation_3D
  use Kind_Parameters, only: RNP
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

contains

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
