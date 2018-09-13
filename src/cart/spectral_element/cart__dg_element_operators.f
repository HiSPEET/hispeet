!> summary:  Element operators for discontinuous cuboidal elements
!> author:   Joerg Stiller
!> date:     2016/03/18
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Element operators for discontinuous cuboidal elements
!===============================================================================

module CART__DG_Element_Operators
  use Kind_Parameters,   only: RNP
  use Execution_Control, only: Error
  use DG_Element_Operators_1D
  implicit none
  private

  public :: DG_ElementOperators

  !-----------------------------------------------------------------------------
  !> Element operators for discontinuous cuboidal elements

  type, extends(DG_ElementOperators1D) :: DG_ElementOperators

    ! discretization parameters
    real(RNP) :: dx(3)       !< element extensions
    real(RNP) :: penalty = 1 !< penalty parameter > 1
    real(RNP) :: mu(3)       !< penalty coefficients

  contains

    generic :: New => New_ElementOperators
    procedure, private :: New_ElementOperators

    procedure :: Adjust => Adjust_ElementOperators
    procedure :: Get_1D_StiffnessMatrix

  end type DG_ElementOperators

contains

!-------------------------------------------------------------------------------
!> Initialize a new DG element operators

subroutine New_ElementOperators(this, po, dx, penalty)
  class(DG_ElementOperators), intent(inout) :: this
  integer,   intent(in) :: po           !< polynomial order
  real(RNP), intent(in) :: dx(3)        !< element extensions
  real(RNP), intent(in) :: penalty      !< penalty parameter > 1

  ! initialize standard operators ..............................................

  call this%New(po)

  ! discretization parameters ..................................................

  this % dx      = dx
  this % penalty = penalty
  this % mu(1)   = this % PenaltyFactor([dx(1), dx(1)], penalty)
  this % mu(2)   = this % PenaltyFactor([dx(2), dx(2)], penalty)
  this % mu(3)   = this % PenaltyFactor([dx(3), dx(3)], penalty)

end subroutine New_ElementOperators

!-------------------------------------------------------------------------------
!> Adjust DG element operators to given element dimensions and/or penalty

subroutine Adjust_ElementOperators(this, dx, penalty)
  class(DG_ElementOperators), intent(inout) :: this
  real(RNP),     optional, intent(in)    :: dx(3)   !< element extensions
  real(RNP),     optional, intent(in)    :: penalty !< penalty parameter > 1

  if (present(dx)) then
    this % dx = dx
  end if

  if (present(penalty)) then
    this % penalty = penalty
  end if

  associate(dx => this%dx, penalty => this%penalty)
    this % mu(1) = this % PenaltyFactor([dx(1), dx(1)], penalty)
    this % mu(2) = this % PenaltyFactor([dx(2), dx(2)], penalty)
    this % mu(3) = this % PenaltyFactor([dx(3), dx(3)], penalty)
  end associate

end subroutine Adjust_ElementOperators

!-------------------------------------------------------------------------------
!> 1D stiffness matrix for given direction and boundary conditions
!>
!> The element stiffness matrix `Le` represents the nontrivial row entries
!> of the global 1D stiffness matrix corresponding to the given element. It
!> must be dimensioned as `Le(0:po,0:po,-1:1)`, where `po = this%po` is the
!> polynomial order. The third index refers to the preceding (-1), current (0)
!> and succeeding (1) element, respectively.

subroutine Get_1D_StiffnessMatrix(this, direction, bc, Le)
  class(DG_ElementOperators), intent(in) :: this
  integer,   intent(in)  :: direction     !< coordinate direction {1,2,3}
  character, intent(in)  :: bc(2)         !< boundary conditions {'','D','N'}
  real(RNP), intent(out) :: Le(0:,0:,-1:) !< 1D element stiffness matrix

  real(RNP) :: dx(-1:1)

  select case(direction)
  case(1:3)
    dx = this % dx(direction)
    call this % GetStiffnessMatrix(dx, this%penalty, bc, Le)
  case default
    call Error('Get_1D_StiffnessMatrix',    &
               'direction must 1,2 or 3',   &
               'CART__DG_Element_Operators' )
  end select

end subroutine Get_1D_StiffnessMatrix

!===============================================================================

end module CART__DG_Element_Operators
