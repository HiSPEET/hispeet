!> summary:  Evaluation of Cartesian diffusion residual, constant+isotropic
!> author:   Joerg Stiller
!> date:     2016/12/12
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Evaluation of Cartesian diffusion residual, constant+isotropic
!===============================================================================

module CART__DG_Diffusion_CI_Residual

  use Kind_Parameters, only: RNP
  use Constants,       only: ZERO, ONE
  use Array_Assignments

  use CART__DG_Element_Operators
  use CART__Mesh_Partition
  use CART__DG_Diffusion_CI_Operator

  implicit none
  private

  public :: DiffusionResidual

contains

!-------------------------------------------------------------------------------
!> Diffusion residual, r = f - (lambda M u + nu L u)

subroutine DiffusionResidual(mesh, eop, lambda, nu, bc, u, f, r)

  ! arguments ..................................................................

  class(MeshPartition),       intent(in) :: mesh !< mesh partition
  class(DG_ElementOperators), intent(in) :: eop  !< DG element operators

  real(RNP), intent(in)  :: lambda      !< Helmholtz parameter
  real(RNP), intent(in)  :: nu          !< diffusivity
  character, intent(in)  :: bc(:)       !< boundary conditions {P,D,N}
  real(RNP), intent(in)  :: u(:,:,:,:)  !< approximate solution
  real(RNP), intent(in)  :: f(:,:,:,:)  !< right hand side
  real(RNP), intent(out) :: r(:,:,:,:)  !< result

  ! evaluation .................................................................

  call DiffusionOperator(mesh, eop, lambda, nu, bc, u, r)

  call MergeArrays(-ONE, r, ONE, f)

end subroutine DiffusionResidual

!===============================================================================

end module CART__DG_Diffusion_CI_Residual
