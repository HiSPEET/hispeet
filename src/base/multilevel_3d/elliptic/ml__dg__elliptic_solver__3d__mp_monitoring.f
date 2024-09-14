submodule(ML__DG__Elliptic_Solver__3D) MP_Monitoring
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Monitoring with given residual

  module subroutine Monitoring_R(this, l, step, r)
    class(ML_DG_EllipticSolver_3D), intent(in) :: this
    integer,               intent(in) :: l          !< level
    character(len=*),      intent(in) :: step       !< current step
    real(RNP), contiguous, intent(in) :: r(:,:,:,:) !< residual

    real(RNP) :: rr

    associate(mesh => this % elliptic_op(l) % sem % mesh)

      if (log_level_multigrid_cycle < 1 .or. mesh%part < 0) return

      rr = ScalarProduct(r, r, mesh%comm_parts)

      !$omp master
      if (mesh%part == 0) then
        print '(2X,A,I2,3A,ES10.3)', 'level ', l, ': r[', step, '] =', sqrt(rr)
      end if
      !$omp end master

    end associate

  end subroutine Monitoring_R

  !-----------------------------------------------------------------------------
  !> Monitoring with constant diffusivity

  module subroutine Monitoring_C(this, l, step, lambda, nu, f, bv, u)
    class(ML_DG_EllipticSolver_3D), intent(in) :: this
    integer,                    intent(in) :: l          !< level
    character(len=*),           intent(in) :: step       !< current step
    real(RNP),                  intent(in) :: lambda     !< λ
    real(RNP),                  intent(in) :: nu         !< diffusivity
    real(RNP), contiguous,      intent(in) :: f(:,:,:,:) !< RHS
    class(BoundaryVariable_3D), intent(in) :: bv(:)      !< boundary values
    real(RNP), contiguous,      intent(in) :: u(:,:,:,:) !< operand

    real(RNP), dimension(:,:,:,:), allocatable, save :: r

    associate(mesh => this % elliptic_op(l) % sem % mesh)

      if (log_level_multigrid_cycle < 1 .or. mesh%part < 0) return

      !$omp master
      allocate(r, mold = u)
      !$omp end master
      !$omp barrier

      call this % Residual(l, lambda, nu, f, bv, u, r)
      call this % Monitoring(l, step, r)

      !$omp master
      deallocate(r)
      !$omp end master

    end associate

  end subroutine Monitoring_C

  !-----------------------------------------------------------------------------
  !> Monitoring with variable diffusivity

  module subroutine Monitoring_V(this, l, step, lambda, nu, f, bv, u)
    class(ML_DG_EllipticSolver_3D), intent(in) :: this
    integer,                    intent(in) :: l           !< level
    character(len=*),           intent(in) :: step        !< current step
    real(RNP),                  intent(in) :: lambda      !< λ
    real(RNP), contiguous,      intent(in) :: nu(:,:,:,:) !< diffusivity
    real(RNP), contiguous,      intent(in) :: f(:,:,:,:)  !< RHS
    real(RNP), contiguous,      intent(in) :: u(:,:,:,:)  !< operand
    class(BoundaryVariable_3D), intent(in) :: bv(:)       !< boundary values

    real(RNP), dimension(:,:,:,:), allocatable, save :: r

    associate(mesh => this % elliptic_op(l) % sem % mesh)

      if (log_level_multigrid_cycle < 1 .or. mesh%part < 0) return

      !$omp master
      allocate(r, mold = u)
      !$omp end master
      !$omp barrier

      call this % Residual(l, lambda, nu, f, bv, u, r)
      call this % Monitoring(l, step, r)

      !$omp master
      deallocate(r)
      !$omp end master

    end associate

  end subroutine Monitoring_V

  !=============================================================================

end submodule MP_Monitoring
