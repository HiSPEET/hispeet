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
    procedure :: Apply
    procedure :: BcToRHS
    procedure :: Residual
    procedure :: ConjugateGradients
    procedure :: SchwarzMethod
  end type EllipticOperator3D_CG

!===============================================================================

contains

!===============================================================================
! Dummy procedures

!--------------------------------------------------------------------------
!> Applies the operator to given approximation: const isotropic

subroutine Apply(this, u, v)
  class(EllipticOperator3D_CG), intent(in) :: this
  real(RNP), intent(in)  :: u(0:,0:,0:,:)  !< approximate solution
  real(RNP), intent(out) :: v(0:,0:,0:,:)  !< result

  v = 0

end subroutine Apply

!--------------------------------------------------------------------------
!> Adds the boundary contributions of the right hand side

subroutine BcToRHS(this, bv, f)
  class(EllipticOperator3D_CG), intent(in)    :: this
  type(BoundaryVariable),       intent(in)    :: bv(:)      !< BC
  real(RNP),                    intent(inout) :: f(:,:,:,:) !< RHS
end subroutine BcToRHS

!--------------------------------------------------------------------------
!> Computes the residual to given approximation

subroutine Residual(this, u, f, r)
  class(EllipticOperator3D_CG), intent(in) :: this
  real(RNP), intent(in)  :: u(0:,0:,0:,:)  !< approximate solution
  real(RNP), intent(in)  :: f(0:,0:,0:,:)  !< right hand side
  real(RNP), intent(out) :: r(0:,0:,0:,:)  !< result
end subroutine Residual

!--------------------------------------------------------------------------
!> Performs iteration sweeps starting from given approximation

subroutine ConjugateGradients(this, u, f, i_max, r_red, r_max, ni)
  class(EllipticOperator3D_CG), intent(in) :: this
  real(RNP), intent(inout) :: u(0:,0:,0:,:)  !< approximate solution
  real(RNP), intent(in)    :: f(0:,0:,0:,:)  !< right hand side
  integer,   intent(in)    :: i_max          !< max num iterations
  real(RNP), optional, intent(in)  :: r_red  !< min residual reduction
  real(RNP), optional, intent(in)  :: r_max  !< max admissible residual
  integer,   optional, intent(out) :: ni     !< exec num iterations
end subroutine ConjugateGradients

!--------------------------------------------------------------------------
!> Performs iteration sweeps starting from given approximation

subroutine SchwarzMethod(this, u, f, i_max, r_red, r_max, ni)
  class(EllipticOperator3D_CG), intent(in) :: this
  real(RNP), intent(inout) :: u(0:,0:,0:,:)  !< approximate solution
  real(RNP), intent(in)    :: f(0:,0:,0:,:)  !< right hand side
  integer,   intent(in)    :: i_max          !< max num iterations
  real(RNP), optional, intent(in)  :: r_red  !< min residual reduction
  real(RNP), optional, intent(in)  :: r_max  !< max admissible residual
  integer,   optional, intent(out) :: ni     !< exec num iterations
end subroutine SchwarzMethod

!===============================================================================

end module CART__Elliptic_Operator_CG
