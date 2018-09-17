!> summary:  Cartesian diffusion operator, constant+isotropic
!> author:   Joerg Stiller
!> date:     2017/07/02
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Cartesian diffusion operator, constant+isotropic
!===============================================================================

module CART__DG_Diffusion_CI_Operator

  use Kind_Parameters, only: RNP
  use Constants,       only: ZERO, ONE
  use Array_Assignments

  use CART__DG_Element_Operators
  use CART__Mesh_Partition
  use CART__DG_Diffusion_CIS_Operator, DiffusionOperator_S => DiffusionOperator
  use CART__DG_Diffusion_CIU_Operator, DiffusionOperator_U => DiffusionOperator

  implicit none
  private

  public :: DiffusionOperator

contains

!-------------------------------------------------------------------------------
!> Diffusion operator, v = lambda M u + nu L u

subroutine DiffusionOperator(mesh, eops, lambda, nu, bc, u, v)

  ! arguments ..................................................................

  class(MeshPartition),       intent(in) :: mesh  !< mesh partition
  class(DG_ElementOperators), intent(in) :: eops  !< DG element oprators

  real(RNP), intent(in)  :: lambda      !< Helmholtz parameter
  real(RNP), intent(in)  :: nu          !< diffusivity
  character, intent(in)  :: bc(:)       !< boundary conditions {P,D,N}
  real(RNP), intent(in)  :: u(:,:,:,:)  !< approximate solution
  real(RNP), intent(out) :: v(:,:,:,:)  !< result

  ! evaluation .................................................................

  if (mesh%structured) then
    call DiffusionOperator_S(mesh, eops, lambda, nu, bc, u, v)
  ! call DiffusionOperator_U(mesh, eops, lambda, nu, bc, u, v)
  else
    call DiffusionOperator_U(mesh, eops, lambda, nu, bc, u, v)
  end if

end subroutine DiffusionOperator

!===============================================================================

end module CART__DG_Diffusion_CI_Operator
