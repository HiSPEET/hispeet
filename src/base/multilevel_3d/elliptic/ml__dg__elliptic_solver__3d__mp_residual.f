submodule(ML__DG__Elliptic_Solver__3D) MP_Residual
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Residual with constant diffusivity
  !>
  !> Homogeneous boundary conditions are assumed if `bv` is absent.
  !> This case is used for implementing the correction scheme

  module subroutine Residual_C(this, l, bc, lambda, nu, f, bv, u, r)
    class(ML_DG_EllipticSolver_3D), intent(in) :: this
    integer,   intent(in) :: l
    character, intent(in) :: bc(:)
    real(RNP), intent(in) :: lambda
    real(RNP), intent(in) :: nu
    real(RNP), contiguous, intent(in) :: f(:,:,:,:)
    class(BoundaryVariable_3D), optional, intent(in) :: bv(:)
    real(RNP), contiguous, intent(in)  :: u(:,:,:,:)
    real(RNP), contiguous, intent(out) :: r(:,:,:,:)

    call  this % elliptic_op(l) % Residual(bc, lambda, nu, f, bv, u, r)

  end subroutine Residual_C

  !-----------------------------------------------------------------------------
  !> Residual with variable diffusivity
  !>
  !> Homogeneous boundary conditions are assumed if `bv` is absent.
  !> This case is used for implementing the correction scheme

  module subroutine Residual_V(this, l, bc, lambda, nu, f, bv, u, r)
    class(ML_DG_EllipticSolver_3D), intent(in) :: this
    integer,   intent(in) :: l
    character, intent(in) :: bc(:)
    real(RNP), intent(in) :: lambda
    real(RNP), contiguous, intent(in) :: nu(:,:,:,:)
    real(RNP), contiguous, intent(in) :: f(:,:,:,:)
    class(BoundaryVariable_3D), optional, intent(in) :: bv(:)
    real(RNP), contiguous, intent(in)  :: u(:,:,:,:)
    real(RNP), contiguous, intent(out) :: r(:,:,:,:)

    call  this % elliptic_op(l) % Residual(bc, lambda, nu, f, bv, u, r)

  end subroutine Residual_V

  !=============================================================================

end submodule MP_Residual
