!> summary:  Projection of 3D mesh variables
!> author:   Joerg Stiller
!> date:     2018/11/02
!> license:  Institute of Flupd Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> @remark
!>   - assumes Cartesian mesh
!>   - element-wise projection, disregarding continuity
!>   - to be replaced by more general version
!===============================================================================

module Projection_Operator__3D
  use Kind_Parameters, only: RNP
  use Constants,       only: THIRD
  use Standard_Element_Operators__1D
  use Projection_Operator__1D
  use TPO__AAA__3D
  implicit none
  private

  !-----------------------------------------------------------------------------
  !> Embedded interpolation of 3D mesh variables

  type, extends(ProjectionOperator_1D), public :: ProjectionOperator_3D

  contains

    generic :: Apply  =>  Project_S, Project_A
    procedure, private :: Project_S
    procedure, private :: Project_A

  end type ProjectionOperator_3D

  ! constructor interface
  interface ProjectionOperator_3D
    module procedure New_SX
  end interface

contains

  !=============================================================================
  ! Constructor

  !-----------------------------------------------------------------------------
  !> New ProjectionOperator_3D

  type(ProjectionOperator_3D) function New_SX(eop, xq, wq, dx) result(this)
    class(StandardElementOperators_1D), intent(in) :: eop !< standard operators
    real(RNP), intent(in) :: xq(:) !< quadrature points
    real(RNP), intent(in) :: wq(:) !< quadrature weights
    real(RNP), intent(in) :: dx(3) !< element extensions

    real(RNP) :: dx_m

    dx_m = product(dx) ** THIRD

    call this % Init_ProjectionOperator_1D(eop, xq, wq, dx_m)

  end function New_SX

  !=============================================================================
  ! Type-bound procedures

  !-----------------------------------------------------------------------------
  !> Projection of scalar variables

  subroutine Project_S(this, uq, up)
    class(ProjectionOperator_3D), intent(in) :: this
    real(RNP), intent(in)  :: uq(:,:,:,:) !< mesh variable at quadrature points
    real(RNP), intent(out) :: up(:,:,:,:) !< projected mesh variable

    call TPO_AAA(this%MA, uq, up)

  end subroutine Project_S

  !-----------------------------------------------------------------------------
  !> Projection of array variables

  subroutine Project_A(this, uq, up)
    class(ProjectionOperator_3D), intent(in) :: this
    real(RNP), intent(in)  :: uq(:,:,:,:,:) !< variable at quadrature points
    real(RNP), intent(out) :: up(:,:,:,:,:) !< projected variable

    integer :: c

    do c = 1, size(uq,5)
      call TPO_AAA(this%MA, uq(:,:,:,:,c), up(:,:,:,:,c))
    end do

  end subroutine Project_A

  !=============================================================================

end module Projection_Operator__3D
