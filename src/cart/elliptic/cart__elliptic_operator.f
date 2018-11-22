!> summary:  Abstract 3D Cartesian elliptic operator
!> author:   Joerg Stiller
!> date:     2018/11/05
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Abstract 3D Cartesian elliptic operator
!===============================================================================

module CART__Elliptic_Operator
  use Kind_Parameters, only: RNP
  use CART__Mesh_Partition
  use CART__Boundary_Variable

  implicit none
  private

  public :: EllipticOperator3D

  !-----------------------------------------------------------------------------
  !> Abstract type accommodating 3D Cartesian elliptic operators

  type, abstract :: EllipticOperator3D
    class(MeshPartition), pointer :: mesh => null()
    ! to be extended by
    ! - element operators
    ! - conjugate gradient method
    ! - Schwarz operator/method
  contains
    private

    generic, public :: Apply => Apply_CI
    procedure(Apply_CI), deferred :: Apply_CI

    generic, public :: BcToRHS => BcToRHS_CI
    procedure(BcToRHS_CI), deferred :: BcToRHS_CI

    generic, public :: Residual => Residual_CI
    procedure(Residual_CI), deferred :: Residual_CI

    generic, public :: ConjugateGradients => ConjugateGradients_CI
    procedure(Iteration_CI), deferred :: ConjugateGradients_CI

    generic, public :: OverlappingSchwarz => OverlappingSchwarz_CI
    procedure(Iteration_CI), deferred :: OverlappingSchwarz_CI

  end type EllipticOperator3D

  abstract interface

    !--------------------------------------------------------------------------
    !> Applies the operator to given approximation: const isotropic

    subroutine Apply_CI(this, lambda, nu, bc, u, v)
      import
      class(EllipticOperator3D), intent(in) :: this
      real(RNP), intent(in)  :: lambda         !< Helmholtz parameter
      real(RNP), intent(in)  :: nu             !< diffusivity
      character, intent(in)  :: bc(:)          !< boundary conditions {P,D,N}
      real(RNP), intent(in)  :: u(0:,0:,0:,:)  !< approximate solution
      real(RNP), intent(out) :: v(0:,0:,0:,:)  !< result
    end subroutine Apply_CI

    !--------------------------------------------------------------------------
    !> Adds the boundary contributions of the right hand side: const isotropic

    subroutine BcToRHS_CI(this, nu, bv, f)
      import
      class(EllipticOperator3D), intent(in)    :: this
      real(RNP),                 intent(in)    :: nu         !< diffusivity
      type(BoundaryVariable),    intent(in)    :: bv(:)      !< BC
      real(RNP),                 intent(inout) :: f(:,:,:,:) !< RHS
    end subroutine BcToRHS_CI

    !--------------------------------------------------------------------------
    !> Computes the residual to given approximation: const isotropic

    subroutine Residual_CI(this, lambda, nu, bc, u, f, r)
      import
      class(EllipticOperator3D), intent(in) :: this
      real(RNP), intent(in)  :: lambda         !< Helmholtz parameter
      real(RNP), intent(in)  :: nu             !< diffusivity
      character, intent(in)  :: bc(:)          !< BC types {P,D,N}
      real(RNP), intent(in)  :: u(0:,0:,0:,:)  !< approximate solution
      real(RNP), intent(in)  :: f(0:,0:,0:,:)  !< right hand side
      real(RNP), intent(out) :: r(0:,0:,0:,:)  !< result
    end subroutine Residual_CI

    !--------------------------------------------------------------------------
    !> Performs iteration sweeps starting from given approx: const isotropic

    subroutine Iteration_CI(this, lambda, nu, bc, u, f, i_max, r_red, r_max, ni)
      import
      class(EllipticOperator3D), intent(in) :: this
      real(RNP), intent(in)    :: lambda         !< Helmholtz parameter
      real(RNP), intent(in)    :: nu             !< diffusivity
      character, intent(in)    :: bc(:)          !< BC types {P,D,N}
      real(RNP), intent(inout) :: u(0:,0:,0:,:)  !< approximate solution
      real(RNP), intent(in)    :: f(0:,0:,0:,:)  !< right hand side
      integer,   intent(in)    :: i_max          !< max num iterations
      real(RNP), optional, intent(in)  :: r_red  !< min residual reduction
      real(RNP), optional, intent(in)  :: r_max  !< max admissible residual
      integer,   optional, intent(out) :: ni     !< exec num iterations
    end subroutine Iteration_CI

  end interface

!===============================================================================

end module CART__Elliptic_Operator
