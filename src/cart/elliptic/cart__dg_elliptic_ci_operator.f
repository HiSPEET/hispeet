!> summary:  Cartesian elliptic operator, constant+isotropic
!> author:   Joerg Stiller
!> date:     2017/07/02
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Cartesian elliptic operator, constant+isotropic
!===============================================================================

module CART__DG_Elliptic_CI_Operator

  use Kind_Parameters, only: RNP
  use Constants,       only: ZERO, ONE
  use Array_Assignments

  use CART__DG_Element_Operators
  use CART__Mesh_Partition
  use CART__DG_Elliptic_CIS_Operator, EllipticOperator_S => EllipticOperator
  use CART__DG_Elliptic_CIU_Operator, EllipticOperator_U => EllipticOperator

  implicit none
  private

  public :: EllipticOperator

contains

!-------------------------------------------------------------------------------
!> Elliptic operator, v = lambda M u + nu L u

subroutine EllipticOperator(mesh, eops, lambda, nu, bc, u, v)

  ! arguments ..................................................................

  class(MeshPartition),         intent(in) :: mesh  !< mesh partition
  class(DG_ElementOperators3D), intent(in) :: eops  !< DG element oprators

  real(RNP), intent(in)  :: lambda      !< Helmholtz parameter
  real(RNP), intent(in)  :: nu          !< diffusivity
  character, intent(in)  :: bc(:)       !< boundary conditions {P,D,N}
  real(RNP), intent(in)  :: u(:,:,:,:)  !< approximate solution
  real(RNP), intent(out) :: v(:,:,:,:)  !< result

  ! evaluation .................................................................

  if (mesh%structured) then
   ! call EllipticOperator_S(mesh, eops, lambda, nu, bc, u, v)
    call EllipticOperator_U(mesh, eops, lambda, nu, bc, u, v)
  else
    call EllipticOperator_U(mesh, eops, lambda, nu, bc, u, v)
  end if

end subroutine EllipticOperator

!===============================================================================

end module CART__DG_Elliptic_CI_Operator
