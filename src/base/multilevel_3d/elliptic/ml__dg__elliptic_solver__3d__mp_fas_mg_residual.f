submodule(ML__DG__Elliptic_Solver__3D) MP_FAS_MG_Residual
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> FAS-MG residual with constant diffusivity

  module subroutine FAS_MG_Residual_C(this, bc, lambda, nu, f, bv, u, r)
    class(ML_DG_EllipticSolver_3D), intent(in)    :: this
    character,                      intent(in)    :: bc(:)
    real(RNP),                      intent(in)    :: lambda
    real(RNP),                      intent(in)    :: nu
    class(ML_MeshVariable_3D),      intent(in)    :: f
    class(ML_BoundaryVariable_3D),  intent(in)    :: bv
    class(ML_MeshVariable_3D),      intent(in)    :: u
    class(ML_MeshVariable_3D),      intent(inout) :: r

    call FAS_MG_Residual_X(this, bc, lambda, nu, null(), f, bv, u, r)

  end subroutine FAS_MG_Residual_C

  !-----------------------------------------------------------------------------
  !> FAS-MG residual with variable diffusivity

  module subroutine FAS_MG_Residual_V(this, bc, lambda, nu, f, bv, u, r)
    class(ML_DG_EllipticSolver_3D), intent(in)    :: this
    character,                      intent(in)    :: bc(:)
    real(RNP),                      intent(in)    :: lambda
    class(ML_MeshVariable_3D),      intent(in)    :: nu
    class(ML_MeshVariable_3D),      intent(in)    :: f
    class(ML_BoundaryVariable_3D),  intent(in)    :: bv
    class(ML_MeshVariable_3D),      intent(in)    :: u
    class(ML_MeshVariable_3D),      intent(inout) :: r

    call FAS_MG_Residual_X(this, bc, lambda, null(), nu, f, bv, u, r)

  end subroutine FAS_MG_Residual_V

  !-----------------------------------------------------------------------------
  !> Generic FAS-MG residual with constant or variable diffusivity

  module subroutine FAS_MG_Residual_X(this, bc, lambda, nu_0, nu_v, f, bv, u, r)
    class(ML_DG_EllipticSolver_3D),      intent(in)    :: this
    character,                           intent(in)    :: bc(:)
    real(RNP),                           intent(in)    :: lambda
    real(RNP),                 optional, intent(in)    :: nu_0
    class(ML_MeshVariable_3D), optional, intent(in)    :: nu_v
    class(ML_MeshVariable_3D),           intent(in)    :: f
    class(ML_BoundaryVariable_3D),       intent(in)    :: bv
    class(ML_MeshVariable_3D),           intent(in)    :: u
    class(ML_MeshVariable_3D),           intent(inout) :: r

    integer :: l

    do l = 1, size(this%elliptic_op)
      associate( bv_l   => bv % level(l  ) % var            &
               , f_l    => f  % level(l  ) % val(:,:,:,:,1) &
               , u_l    => u  % level(l  ) % val(:,:,:,:,1) &
               , r_l    => r  % level(l  ) % val(:,:,:,:,1) )

        if (present(nu_0)) then
          call this % Residual(l, bc, lambda, nu_0, f_l, bv_l, u_l, r_l)
        else
          associate(nu_l => nu_v % level(l) % val(:,:,:,:,1))
            call this % Residual(l, bc, lambda, nu_l, f_l, bv_l, u_l, r_l)
          end associate
        end if

      end associate
    end do

  end subroutine FAS_MG_Residual_X

  !=============================================================================

end submodule MP_FAS_MG_Residual
