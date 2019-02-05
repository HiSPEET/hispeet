!> summary:  Legendre projection of 3D mesh variables
!> author:   Joerg Stiller
!> date:     2019/02/05
!> license:  Institute of Flupd Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Legendre projection of 3D mesh variables
!===============================================================================

module Legendre_Projection_3D
  use Kind_Parameters, only: RNP
  use Standard_Operators_1D
  use Legendre_Projection
  use TPO_AAA
  implicit none
  private

  public :: LegendreProjectionOperator3D

  !-----------------------------------------------------------------------------
  !> Embedded interpolation of 3D mesh variables

  type, extends(LegendreProjectionOperator) :: LegendreProjectionOperator3D

  contains

    generic :: Apply  =>  Project_S, Project_A
    procedure, private :: Project_S
    procedure, private :: Project_A

  end type LegendreProjectionOperator3D

  ! constructor interface
  interface LegendreProjectionOperator3D
    module procedure New_LPO
  end interface

contains

!===============================================================================
! Constructors

!-------------------------------------------------------------------------------
!> New LegendreProjectionOperator3D

type(LegendreProjectionOperator3D) function New_LPO(eop, pp) result(this)
  class(StandardOperators1D), intent(in) :: eop  !< standard operators
  integer,                    intent(in) :: pp   !< order to project to

  call this % Init_LegendreProjectionOperator(eop, pp)

end function New_LPO

!===============================================================================
! Type-bound procedures

!-------------------------------------------------------------------------------
!> Legendre projection of scalar variables

subroutine Project_S(this, uo, up)
  class(LegendreProjectionOperator3D), intent(in) :: this
  real(RNP), intent(in)  :: uo(:,:,:,:) !< original mesh variable
  real(RNP), intent(out) :: up(:,:,:,:) !< projected mesh variable

  call TPO_AAA_Eval(this%np, this%no, size(uo,4), this%pop, uo, up)

end subroutine Project_S

!-------------------------------------------------------------------------------
!> Legendre projection of array variables

subroutine Project_A(this, uo, up)
  class(LegendreProjectionOperator3D), intent(in) :: this
  real(RNP), intent(in)  :: uo(:,:,:,:,:) !< original mesh variable
  real(RNP), intent(out) :: up(:,:,:,:,:) !< projected mesh variable

  call TPO_AAA_Eval(this%np, this%no, size(uo,4)*size(uo,5), this%pop, uo, up)

end subroutine Project_A

!===============================================================================

end module Legendre_Projection_3D
