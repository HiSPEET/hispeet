!> summary:  Abstract 3D Cartesian elliptic operator
!> author:   Joerg Stiller
!> date:     2018/11/05
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Abstract 3D Cartesian elliptic operator
!===============================================================================

module CART__Elliptic_Operator
  use Kind_Parameters, only: RNP
  use Standard_Operators_1D
  use CART__Mesh_Partition
  use CART__Boundary_Variable
  use CART__Schwarz_Operator

  implicit none
  private

  public :: EllipticOperator3D

  !-----------------------------------------------------------------------------
  !> Abstract type accommodating 3D Cartesian elliptic operators

  type, abstract :: EllipticOperator3D

    class(MeshPartition), pointer :: mesh => null()

    real(RNP) :: lambda = 0                  !< Helmholtz parameter
    real(RNP), allocatable :: nu_ci          !< constant isotropic diffusivity
    real(RNP), allocatable :: nu_vi(:,:,:,:) !< variable isotropic diffusivity
    character, allocatable :: bc(:)          !< boundary conditions {P,D,N}

    class(StandardOperators1D), allocatable :: eop     !< standard SE operators
    type(SchwarzOperator3D),    allocatable :: schwarz !< Schwarz operator

  contains

    generic :: SetProblem => SetProblem_CI, SetProblem_VI
    procedure(SetProblem_CI), deferred :: SetProblem_CI  ! should be private
    procedure(SetProblem_VI), deferred :: SetProblem_VI  ! but fails with ifort

    procedure(Apply),     deferred :: Apply
    procedure(BcToRHS),   deferred :: BcToRHS
    procedure(Residual),  deferred :: Residual
    procedure(Iteration), deferred :: ConjugateGradients
    procedure(Iteration), deferred :: SchwarzMethod

  end type EllipticOperator3D

  abstract interface

    !--------------------------------------------------------------------------
    !> (Re)Set problem parameters for constant isotropic viscosity

    subroutine SetProblem_CI(this, lambda, nu, bc)
      import
      class(EllipticOperator3D), intent(inout) :: this
      real(RNP), intent(in) :: lambda !< Helmholtz parameter
      real(RNP), intent(in) :: nu     !< diffusivity
      character, intent(in) :: bc(:)  !< BC {'D','N','P'}
    end subroutine SetProblem_CI

    !--------------------------------------------------------------------------
    !> (Re)Set problem parameters for variable isotropic viscosity

    subroutine SetProblem_VI(this, lambda, nu, bc)
      import
      class(EllipticOperator3D), intent(inout) :: this
      real(RNP), intent(in) :: lambda         !< Helmholtz parameter
      real(RNP), intent(in) :: nu(0:,0:,0:,:) !< diffusivity
      character, intent(in) :: bc(:)          !< BC {'D','N','P'}
    end subroutine SetProblem_VI

    !--------------------------------------------------------------------------
    !> Applies the operator to given approximation: const isotropic

    subroutine Apply(this, u, v)
      import
      class(EllipticOperator3D), intent(in) :: this
      real(RNP), intent(in)  :: u(0:,0:,0:,:)  !< approximate solution
      real(RNP), intent(out) :: v(0:,0:,0:,:)  !< result
    end subroutine Apply

    !--------------------------------------------------------------------------
    !> Adds the boundary contributions of the right hand side

    subroutine BcToRHS(this, bv, c, f)
      import
      class(EllipticOperator3D), intent(in)    :: this
      type(BoundaryVariable),    intent(in)    :: bv(:)      !< boundary values
      integer,         optional, intent(in)    :: c          !< component [1]
      real(RNP),                 intent(inout) :: f(:,:,:,:) !< RHS
    end subroutine BcToRHS

    !--------------------------------------------------------------------------
    !> Computes the residual to given approximation: const isotropic

    subroutine Residual(this, u, f, r)
      import
      class(EllipticOperator3D), intent(in) :: this
      real(RNP), intent(in)  :: u(0:,0:,0:,:)  !< approximate solution
      real(RNP), intent(in)  :: f(0:,0:,0:,:)  !< right hand side
      real(RNP), intent(out) :: r(0:,0:,0:,:)  !< result
    end subroutine Residual

    !--------------------------------------------------------------------------
    !> Performs iteration sweeps starting from given approximation

    subroutine Iteration(this, u, f, i_max, r_red, r_max, ni)
      import
      class(EllipticOperator3D), intent(in) :: this
      real(RNP), intent(inout) :: u(0:,0:,0:,:)  !< approximate solution
      real(RNP), intent(in)    :: f(0:,0:,0:,:)  !< right hand side
      integer,   intent(in)    :: i_max          !< max num iterations
      real(RNP), optional, intent(in)  :: r_red  !< min residual reduction
      real(RNP), optional, intent(in)  :: r_max  !< max admissible residual
      integer,   optional, intent(out) :: ni     !< exec num iterations
    end subroutine Iteration

  end interface

!===============================================================================

end module CART__Elliptic_Operator
