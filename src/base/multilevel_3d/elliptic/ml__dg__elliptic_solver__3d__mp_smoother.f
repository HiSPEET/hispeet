submodule(ML__DG__Elliptic_Solver__3D) MP_Smoother
  use Boundary_Variable__3D
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Smoother with constant diffusivity

  module subroutine Smoother_C(this, lambda, nu, u, f, bv, n_s)
    class(ML_DG_EllipticSolver_3D), intent(in) :: this
    real(RNP), intent(in) :: lambda
    real(RNP), intent(in) :: nu
    real(RNP), contiguous, intent(inout) :: u(:,:,:,:)
    real(RNP), contiguous, intent(in) :: f(:,:,:,:)
    class(BoundaryVariable_3D), intent(in) :: bv(:)
    integer, intent(in) :: n_s

    call Smoother_X(this, lambda, nu, null(), u, f, bv, n_s)

  end subroutine Smoother_C

  !-----------------------------------------------------------------------------
  !> Smoother with variable diffusivity

  module subroutine Smoother_V(this, lambda, nu, u, f, bv, n_s)
    class(ML_DG_EllipticSolver_3D), intent(in) :: this
    real(RNP), intent(in) :: lambda
    real(RNP), contiguous, intent(in) :: nu(:,:,:,:)
    real(RNP), contiguous, intent(inout) :: u(:,:,:,:)
    real(RNP), contiguous, intent(in) :: f(:,:,:,:)
    class(BoundaryVariable_3D), intent(in) :: bv(:)
    integer, intent(in) :: n_s

    call Smoother_X(this, lambda, null(), nu, u, f, bv, n_s)

  end subroutine Smoother_V

  !-----------------------------------------------------------------------------
  !> Generic smoother with constant or variable diffusivity

  subroutine Smoother_X(this, lambda, nu_c, nu_v, u, f, bv, n_s)
    class(ML_DG_EllipticSolver_3D), intent(in) :: this
    real(RNP), intent(in) :: lambda
    real(RNP), optional, intent(in) :: nu_c
    real(RNP), contiguous, optional, intent(in) :: nu_v(:,:,:,:)
    real(RNP), contiguous, intent(inout) :: u(:,:,:,:)
    real(RNP), contiguous, intent(in) :: f(:,:,:,:)
    class(BoundaryVariable_3D), intent(in) :: bv(:)
    integer, intent(in) :: n_s

    select case(this % smooth_method)
    case(SOLVER_CG)
      ! Flexible CG
      call this % elliptic_op(1) % CG_Method_X &
             (lambda, nu_c, nu_v, u, f, bv, n_s)
    case(SOLVER_WS)
      ! Weighted Additive Schwarz
      call this % elliptic_op(1) % Schwarz_Method_X &
             (lambda, nu_c, nu_v, u, f, bv, n_s)
    case(SOLVER_SPCG)
      ! Schwarz-preconditioned flexible CG
      call this % elliptic_op(1) % SchwarzPCG_Method_X &
             (lambda, nu_c, nu_v, u, f, bv, n_s)
    end select

  end subroutine Smoother_X

  !=============================================================================

end submodule MP_Smoother
