!> summary:  3D Cartesian elliptic operator for interior penalty (IP) DG-SEM
!> author:   Joerg Stiller
!> date:     2018/11/05
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### 3D Cartesian elliptic operator for interior penalty (IP) DG-SEM
!===============================================================================

module CART__Elliptic_Operator_IP
  use Kind_Parameters, only: RNP
  use IP_Element_Operators_1D
  use CART__Elliptic_Operator

! remove if not needed with complete version
use CART__Mesh_Partition
use CART__Boundary_Variable

  implicit none
  private

  public :: EllipticOperator3D_IP

  type, extends(EllipticOperator3D) :: EllipticOperator3D_IP

    type(IP_ElementOperators1D) :: eop    !< 1D operators for IP-DG-SEM
    real(RNP), allocatable :: nu_f(:,:,:) !< max diffusivity on faces

  contains
    private

    procedure :: Apply
    procedure :: BcToRHS
    procedure :: Residual
    procedure :: ConjugateGradients

    !procedure :: OverlappingSchwarz_CI

  end type EllipticOperator3D_IP

  !=============================================================================
  ! Separate procedures

  interface

    !---------------------------------------------------------------------------
    !> Application of the IP/DG elliptic operator

    module subroutine Apply(this, bc, u, v)
      class(EllipticOperator3D_IP), intent(in) :: this
      character, intent(in)  :: bc(:)          !< boundary conditions {P,D,N}
      real(RNP), intent(in)  :: u(0:,0:,0:,:)  !< approximate solution
      real(RNP), intent(out) :: v(0:,0:,0:,:)  !< result
    end subroutine Apply

    !---------------------------------------------------------------------------
    !> Adds the boundary contributions of the right hand side

    module subroutine BcToRHS(this, bv, f)
      class(EllipticOperator3D_IP), intent(in)    :: this
      type(BoundaryVariable),       intent(in)    :: bv(:)      !< BC
      real(RNP),                    intent(inout) :: f(:,:,:,:) !< RHS
    end subroutine BcToRHS

    !--------------------------------------------------------------------------
    !> Computes the residual to given approximation

    module subroutine Residual(this, bc, u, f, r)
      class(EllipticOperator3D_IP), intent(in) :: this
      character, intent(in)  :: bc(:)          !< BC types {P,D,N}
      real(RNP), intent(in)  :: u(0:,0:,0:,:)  !< approximate solution
      real(RNP), intent(in)  :: f(0:,0:,0:,:)  !< right hand side
      real(RNP), intent(out) :: r(0:,0:,0:,:)  !< result
    end subroutine Residual

  end interface

  !=============================================================================
  !> Element-centered overlapping Schwarz method

  interface

    !---------------------------------------------------------------------------
    !> Element-centered overlapping Schwarz method with constant coefficients

    module subroutine OverlappingSchwarz(this, bc, u, f, i_max, r_red, r_max, ni)
      class(EllipticOperator3D_IP), intent(in) :: this
      character, intent(in)    :: bc(:)          !< BC types {P,D,N}
      real(RNP), intent(inout) :: u(0:,0:,0:,:)  !< approximate solution
      real(RNP), intent(in)    :: f(0:,0:,0:,:)  !< right hand side
      integer,   intent(in)    :: i_max          !< max num iterations
      real(RNP), optional, intent(in)  :: r_red  !< min residual reduction
      real(RNP), optional, intent(in)  :: r_max  !< max admissible residual
      integer,   optional, intent(out) :: ni     !< exec num iterations
    end subroutine OverlappingSchwarz

  end interface

!===============================================================================

contains

!===============================================================================
! Dummy procedures

!--------------------------------------------------------------------------
!> Performs iteration sweeps starting from given approx: const isotropic

subroutine ConjugateGradients(this, bc, u, f, i_max, r_red, r_max, ni)
  class(EllipticOperator3D_IP), intent(in) :: this
  character, intent(in)    :: bc(:)          !< BC types {P,D,N}
  real(RNP), intent(inout) :: u(0:,0:,0:,:)  !< approximate solution
  real(RNP), intent(in)    :: f(0:,0:,0:,:)  !< right hand side
  integer,   intent(in)    :: i_max          !< max num iterations
  real(RNP), optional, intent(in)  :: r_red  !< min residual reduction
  real(RNP), optional, intent(in)  :: r_max  !< max admissible residual
  integer,   optional, intent(out) :: ni     !< exec num iterations
end subroutine ConjugateGradients

!===============================================================================

end module CART__Elliptic_Operator_IP
