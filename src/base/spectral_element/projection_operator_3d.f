!> summary:  Projection of 3D mesh variables
!> author:   Joerg Stiller
!> date:     2018/11/02
!> license:  Institute of Flupd Mechanpcs, TU Dresden, 01062 Dresden, Germany
!>
!>### Projection of 3D mesh variables
!===============================================================================

module Projection_Operator_3D
  use Kind_Parameters, only: RNP
  use Constants,       only: THIRD
  use Standard_Operators_1D
  use Projection_Operator
  use TPO_AAA
  implicit none
  private

  !-----------------------------------------------------------------------------
  !> Embedded interpolation of 3D mesh variables

  type, extends(ProjectionOperator), public :: ProjectionOperator3D

  contains

    generic :: Apply  =>  Project_S, Project_A
    procedure, private :: Project_S
    procedure, private :: Project_A

  end type ProjectionOperator3D

  ! constructor interface
  interface ProjectionOperator3D
    module procedure New_SX
  end interface

contains

!===============================================================================
! Constructor

!-------------------------------------------------------------------------------
!> New ProjectionOperator3D

type(ProjectionOperator3D) function New_SX(eop, xq, wq, dx) result(this)
  class(StandardOperators1D), intent(in) :: eop   !< standard operators
  real(RNP),                  intent(in) :: xq(:) !< quadrature points
  real(RNP),                  intent(in) :: wq(:) !< quadrature weights
  real(RNP),                  intent(in) :: dx(3) !< element extensions

  real(RNP) :: dx_m

  dx_m = product(dx) ** THIRD

  call this % Init_ProjectionOperator(eop, xq, wq, dx_m)

end function New_SX

!===============================================================================
! Type-bound procedures

!-------------------------------------------------------------------------------
!> Projection of scalar variables

subroutine Project_S(this, uq, up)
  class(ProjectionOperator3D), intent(in) :: this
  real(RNP), intent(in)  :: uq(:,:,:,:) !< mesh variable at quadrature points
  real(RNP), intent(out) :: up(:,:,:,:) !< projected mesh variable

  call TPO_AAA_Eval(this%np, this%nq, size(uq,4), this%MA, uq, up)

end subroutine Project_S

!-------------------------------------------------------------------------------
!> Projection of array variables

subroutine Project_A(this, uq, up)
  class(ProjectionOperator3D), intent(in) :: this
  real(RNP), intent(in)  :: uq(:,:,:,:,:) !< mesh variable at quadrature points
  real(RNP), intent(out) :: up(:,:,:,:,:) !< projected mesh variable

  call TPO_AAA_Eval(this%np, this%nq, size(uq,4)*size(uq,5), this%MA, uq, up)

end subroutine Project_A

!===============================================================================

end module Projection_Operator_3D
