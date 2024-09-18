submodule(ML__DG__Elliptic_Solver__3D) MP_Residual
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Residual with constant diffusivity

  module subroutine Residual_C(this, l, lambda, nu, f, bv, u, r)
    class(ML_DG_EllipticSolver_3D), intent(in) :: this
    integer,   intent(in) :: l
    real(RNP), intent(in) :: lambda
    real(RNP), intent(in) :: nu
    real(RNP), contiguous, intent(in) :: f(:,:,:,:)
    class(BoundaryVariable_3D), intent(in) :: bv(:)
    real(RNP), contiguous, intent(in)  :: u(:,:,:,:)
    real(RNP), contiguous, intent(out) :: r(:,:,:,:)

    call  this % elliptic_op(l) % Residual(lambda, nu, f, bv, u, r)

  end subroutine Residual_C

  !-----------------------------------------------------------------------------
  !> Residual with variable diffusivity

  module subroutine Residual_V(this, l, lambda, nu, f, bv, u, r)
    class(ML_DG_EllipticSolver_3D), intent(in) :: this
    integer,   intent(in) :: l
    real(RNP), intent(in) :: lambda
    real(RNP), contiguous, intent(in) :: nu(:,:,:,:)
    real(RNP), contiguous, intent(in) :: f(:,:,:,:)
    class(BoundaryVariable_3D), intent(in) :: bv(:)
    real(RNP), contiguous, intent(in)  :: u(:,:,:,:)
    real(RNP), contiguous, intent(out) :: r(:,:,:,:)

    call  this % elliptic_op(l) % Residual(lambda, nu, f, bv, u, r)

  end subroutine Residual_V

  !=============================================================================

end submodule MP_Residual
