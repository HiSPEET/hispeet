!> summary:  Element operators for continuous cuboidal elements
!> author:   Joerg Stiller
!> date:     2018/10/22
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Element operators for continuous cuboidal elements
!===============================================================================

module CART__CG_Element_Operators
  use Kind_Parameters,   only: RNP
  use Execution_Control, only: Error
  use CG_Element_Operators_1D
  implicit none
  private

  public :: CG_ElementOperators3D

  !-----------------------------------------------------------------------------
  !> Element operators for continuous cuboidal elements

  type, extends(CG_ElementOperators1D) :: CG_ElementOperators3D

    real(RNP) :: dx(3)  !< element extensions

  contains

    procedure :: Init_CG_ElementOperators3D
    procedure :: Get_1D_StiffnessMatrix

  end type CG_ElementOperators3D

contains

!-------------------------------------------------------------------------------
!> Initialize a new CG element operators

subroutine Init_CG_ElementOperators3D(this, po, dx)
  class(CG_ElementOperators3D), intent(inout) :: this
  integer,   intent(in) :: po           !< polynomial order
  real(RNP), intent(in) :: dx(3)        !< element extensions

  call this % Init_CG_ElementOperators1D(CG_ElementOptions1D(po))
  this % dx = dx

end subroutine Init_CG_ElementOperators3D

!-------------------------------------------------------------------------------
!> 1D stiffness matrix for given direction and boundary conditions
!>
!> The element stiffness matrix `Le` represents the nontrivial row entries
!> of the global 1D stiffness matrix corresponding to the given element. It
!> must be dimensioned as `Le(0:po,0:po,-1:1)`, where `po = this%po` is the
!> polynomial order. The third index refers to the preceding (-1), current (0)
!> and succeeding (1) element, respectively.

subroutine Get_1D_StiffnessMatrix(this, direction, bc, Le)
  class(CG_ElementOperators3D), intent(in) :: this
  integer,   intent(in)  :: direction     !< coordinate direction {1,2,3}
  character, intent(in)  :: bc(2)         !< boundary conditions {'','D','N'}
  real(RNP), intent(out) :: Le(0:,0:,-1:) !< 1D element stiffness matrix

  real(RNP) :: dx(-1:1)

  select case(direction)
  case(1:3)
    dx = this % dx(direction)
    call this % GetStiffnessMatrix(dx, bc, Le)
  case default
    call Error('Get_1D_StiffnessMatrix',    &
               'direction must 1,2 or 3',   &
               'CART__CG_Element_Operators' )
  end select

end subroutine Get_1D_StiffnessMatrix

!===============================================================================

end module CART__CG_Element_Operators
