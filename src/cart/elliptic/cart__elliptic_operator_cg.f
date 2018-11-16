!> summary:  3D Cartesian elliptic operator for CG-SEM
!> author:   Joerg Stiller
!> date:     2018/11/05
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### 3D Cartesian elliptic operator for CG-SEM
!===============================================================================

module CART__Elliptic_Operator_CG
  use Kind_Parameters, only: RNP
  use CG_Element_Operators_1D
  use CART__Elliptic_Operator

! remove if not needed with complete version
use CART__Mesh_Partition
use CART__Boundary_Variable

  implicit none
  private

  type, extends(EllipticOperator3D) :: EllipticOperator3D_CG
    type(CG_ElementOperators1D) :: eop !< 1D OPs for CG-SEM
  contains
    private
    procedure :: Apply_CI
    procedure :: BcToRHS_CI
    procedure :: Residual_CI
    procedure :: ConjugateGradients_CI
    procedure :: OverlappingSchwarz_CI
  end type EllipticOperator3D_CG

!===============================================================================

contains

!===============================================================================
! Dummy procedures

!--------------------------------------------------------------------------
!> Applies the operator to given approximation: const isotropic

subroutine Apply_CI(this, lambda, nu, bc, u, v)
  class(EllipticOperator3D_CG), intent(in) :: this
  real(RNP), intent(in)  :: lambda         !< Helmholtz parameter
  real(RNP), intent(in)  :: nu             !< diffusivity
  character, intent(in)  :: bc(:)          !< boundary conditions {P,D,N}
  real(RNP), intent(in)  :: u(0:,0:,0:,:)  !< approximate solution
  real(RNP), intent(out) :: v(0:,0:,0:,:)  !< result

  v = 0

end subroutine Apply_CI

!--------------------------------------------------------------------------
!> Adds the boundary contributions of the right hand side: const isotropic

subroutine BcToRHS_CI(this, nu, bv, f)
  class(EllipticOperator3D_CG), intent(in) :: this
  real(RNP),              intent(in)    :: nu    !< diffusivity
  type(BoundaryVariable), intent(in)    :: bv(:) !< boundary conditions
  real(RNP),              intent(inout) :: f     !< RHS
end subroutine BcToRHS_CI

!--------------------------------------------------------------------------
!> Computes the residual to given approximation: const isotropic

subroutine Residual_CI(this, lambda, nu, bc, u, f, r)
  class(EllipticOperator3D_CG), intent(in) :: this
  real(RNP), intent(in)  :: lambda         !< Helmholtz parameter
  real(RNP), intent(in)  :: nu             !< diffusivity
  character, intent(in)  :: bc(:)          !< BC types {P,D,N}
  real(RNP), intent(in)  :: u(0:,0:,0:,:)  !< approximate solution
  real(RNP), intent(in)  :: f(0:,0:,0:,:)  !< right hand side
  real(RNP), intent(out) :: r(0:,0:,0:,:)  !< result
end subroutine Residual_CI

!--------------------------------------------------------------------------
!> Performs iteration sweeps starting from given approx: const isotropic

subroutine ConjugateGradients_CI(this, lambda, nu, bc, u, f, i_max, r_red, r_max, ni)
  class(EllipticOperator3D_CG), intent(in) :: this
  real(RNP), intent(in)    :: lambda         !< Helmholtz parameter
  real(RNP), intent(in)    :: nu             !< diffusivity
  character, intent(in)    :: bc(:)          !< BC types {P,D,N}
  real(RNP), intent(inout) :: u(0:,0:,0:,:)  !< approximate solution
  real(RNP), intent(in)    :: f(0:,0:,0:,:)  !< right hand side
  integer,   intent(in)    :: i_max          !< max num iterations
  real(RNP), optional, intent(in)  :: r_red  !< min residual reduction
  real(RNP), optional, intent(in)  :: r_max  !< max admissible residual
  integer,   optional, intent(out) :: ni     !< exec num iterations
end subroutine ConjugateGradients_CI

!--------------------------------------------------------------------------
!> Performs iteration sweeps starting from given approx: const isotropic

subroutine OverlappingSchwarz_CI(this, lambda, nu, bc, u, f, i_max, r_red, r_max, ni)
  class(EllipticOperator3D_CG), intent(in) :: this
  real(RNP), intent(in)    :: lambda         !< Helmholtz parameter
  real(RNP), intent(in)    :: nu             !< diffusivity
  character, intent(in)    :: bc(:)          !< BC types {P,D,N}
  real(RNP), intent(inout) :: u(0:,0:,0:,:)  !< approximate solution
  real(RNP), intent(in)    :: f(0:,0:,0:,:)  !< right hand side
  integer,   intent(in)    :: i_max          !< max num iterations
  real(RNP), optional, intent(in)  :: r_red  !< min residual reduction
  real(RNP), optional, intent(in)  :: r_max  !< max admissible residual
  integer,   optional, intent(out) :: ni     !< exec num iterations
end subroutine OverlappingSchwarz_CI

end module CART__Elliptic_Operator_CG
