submodule(ML__DG__Elliptic_Solver__3D) MP_CS_MG_Solver
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> CS-MG solver for problems with global refinement and constant diffusivity

  module subroutine CS_MG_Solver_C(this, lambda, nu, u, f, bv, ni, r_2)
    class(ML_DG_EllipticSolver_3D), intent(in) :: this
    real(RNP), intent(in) :: lambda                 !< Helmholtz parameter
    real(RNP), intent(in) :: nu                     !< diffusivity
    class(ML_MeshVariable_3D), intent(inout) :: u   !< approx/final solution
    class(ML_MeshVariable_3D), intent(inout) :: f   !< RHS
    class(ML_BoundaryVariable_3D), intent(in) :: bv !< boundary values
    integer,   optional, intent(out) :: ni          !< num executed cycles
    real(RNP), optional, intent(out) :: r_2         !< Euclidean residual norm

    call CS_MG_Solver_X(this, lambda, nu, null(), u, f, bv, ni, r_2)

  end subroutine CS_MG_Solver_C

  !-----------------------------------------------------------------------------
  !> CS-MG solver for problems with global refinement and variable diffusivity

  module subroutine CS_MG_Solver_V(this, lambda, nu, u, f, bv, ni, r_2)
    class(ML_DG_EllipticSolver_3D), intent(in) :: this
    real(RNP), intent(in) :: lambda                 !< Helmholtz parameter
    class(ML_MeshVariable_3D), intent(in) :: nu     !< diffusivity
    class(ML_MeshVariable_3D), intent(inout) :: u   !< approx/final solution
    class(ML_MeshVariable_3D), intent(inout) :: f   !< RHS
    class(ML_BoundaryVariable_3D), intent(in) :: bv !< boundary values
    integer, optional, intent(out) :: ni            !< num executed cycles
    real(RNP), optional, intent(out) :: r_2         !< Euclidean residual norm

    call CS_MG_Solver_X(this, lambda, null(), nu, u, f, bv, ni, r_2)

  end subroutine CS_MG_Solver_V

  !-----------------------------------------------------------------------------
  !> Generic CS-MG solver for problems with constant or variable diffusivity
  !>
  !> Either `nu_0` or `nu_v` must be given.

  module subroutine CS_MG_Solver_X(this, lambda, nu_0, nu_v, u, f, bv, ni, r_2)
    class(ML_DG_EllipticSolver_3D), intent(in) :: this
    real(RNP), intent(in) :: lambda                 !< Helmholtz parameter
    real(RNP), intent(in) :: nu_0                   !< constant diffusivity
    class(ML_MeshVariable_3D), intent(in) :: nu_v   !< variable diffusivity
    class(ML_MeshVariable_3D), intent(inout) :: u   !< approx/final solution
    class(ML_MeshVariable_3D), intent(inout) :: f   !< RHS
    class(ML_BoundaryVariable_3D), &
                          target, intent(in) :: bv  !< boundary values
    integer, optional, intent(out) :: ni            !< num executed cycles
    real(RNP), optional, intent(out) :: r_2         !< Euclidean residual norm

    optional :: nu_0, nu_v

    ! internal variables :::::::::::::::::::::::::::::::::::::::::::::::::::::::

    type(ML_MeshVariable_3D), allocatable, save :: r
    type(BoundaryVariable_3D), pointer, save :: bv_l(:), bv_1(:)
    integer, save :: l_top
    logical, save :: converged

    real(RNP) :: rr, r_max, r_new, r_old
    logical :: check_convergence
    integer :: e, l, m, n

    associate( sem    => this % ml_op % sem      &
             , iop_cf => this % ml_op % iop_cf_x &
             , iop_fc => this % ml_op % iop_fc_x &
             , ell_op => this % elliptic_op      )

      ! initialization :::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      check_convergence = max(this%r_red, this%r_max) > 0

      !$omp master
      l_top = size(sem)
      call r % Init(this%ml_op, nc=1)
      bv_l => bv % level(l_top) % var
      if (l_top == 1) then
        bv_1 => bv_l
      else
        bv_1 => null()
      end if
      !$omp end master

      do l = 1, l_top
        associate(po => sem(l) % std_op % po)
          !$omp do
          do e = 1, sem(l) % mesh % n_elem
            r % level(l) % val(0:po,0:po,0:po,e,1) = ZERO
          end do
          !$omp end do nowait
        end associate
      end do
      !$omp barrier

      ! termination conditions
      if (check_convergence) then
        l = l_top
        associate( f_l => f % level(l) % val(:,:,:,:,1) &
                 , u_l => u % level(l) % val(:,:,:,:,1) &
                 , r_l => r % level(l) % val(:,:,:,:,1) )

          ! initial residual
          if (present(nu_0)) then
            call this % Residual(l, lambda, nu_0, f_l, bv_l, u_l, r_l)
          else
            associate(nu_l => nu_v % level(l) % val(:,:,:,:,1))
              call this % Residual(l, lambda, nu_l, f_l, bv_l, u_l, r_l)
            end associate
          end if

          rr = ScalarProduct(r_l, r_l)
          r_old  = sqrt(rr)
          r_max  = max(r_old * this%r_red, this%r_max)
          !$omp master
          converged = r_old < r_max
          call XMPI_Bcast(converged, root = 0, comm = sem(1)%mesh%comm_world)
          !$omp end master
          !$omp barrier
        end associate
      else
        !$omp single
        converged = .false.
        !$omp end single
      end if

      if (converged) then
        !$omp master
        deallocate(r)
        if (present(ni)) ni = 0
        if (present(r_2)) r_2 = r_old
        !$omp end master
        return
      end if

      ! V cycles :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      V_OUTER: do m = 1, this % i_max

        V_DOWN: do l = l_top, 2, -1

          associate( mesh_l => sem(l  ) % mesh                 &
                   , mesh_p => sem(l-1) % mesh                 &
                   , f_l    => f % level(l  ) % val(:,:,:,:,1) &
                   , u_l    => u % level(l  ) % val(:,:,:,:,1) &
                   , r_l    => r % level(l  ) % val(:,:,:,:,1) &
                   , f_p    => f % level(l-1) % val(:,:,:,:,1) )

            ! pre-smoothing and residual computation .........................

            if (l < l_top) then
              call SetArray(u_l, ZERO)
            end if

            if (present(nu_0)) then
              call this % Monitoring(l, '0', lambda, nu_0, f_l, bv_l, u_l)
            else
              associate(nu_l => nu_v % level(l) % val(:,:,:,:,1))
                call this % Monitoring(l, '0', lambda, nu_l, f_l, bv_l, u_l)
              end associate
            end if

            n = this % ns_1
            if (present(nu_0)) then
              call this % Smoother(l, lambda, nu_0, u_l, f_l, bv_l, n)
              call this % Residual(l, lambda, nu_0, f_l, bv_l, u_l, r_l)
            else
              associate(nu_l => nu_v % level(l) % val(:,:,:,:,1))
                call this % Smoother(l, lambda, nu_l, u_l, f_l, bv_l, n)
                call this % Residual(l, lambda, nu_l, f_l, bv_l, u_l, r_l)
              end associate
            end if
            call this % Monitoring(l, '1', r_l)

            ! nullify boundary variables for l < t_top .......................

            !$omp master
            bv_l => null()
            !$omp end master

            ! restriction ....................................................

            call ChildToParentRestriction_3D &
                     (mesh_l, mesh_p, iop_cf(l-1), r_l, f_p)

          end associate
        end do V_DOWN

        V_COARSE: associate( f_1 => f % level(1) % val(:,:,:,:,1) &
                           , u_1 => u % level(1) % val(:,:,:,:,1) )

          ! coarse grid solver ...............................................

          if (present(nu_0)) then
            call this % Monitoring(1, '0', lambda, nu_0, f_1, bv_1, u_1)
            call this % CoarseSolver(lambda, nu_0, u_1, f_1, bv_1)
            call this % Monitoring(1, 's', lambda, nu_0, f_1, bv_1, u_1)
          else
            associate(nu_1 => nu_v % level(1) % val(:,:,:,:,1))
              call this % Monitoring(1, '0', lambda, nu_1, f_1, bv_1, u_1)
              call this % CoarseSolver(lambda, nu_1, u_1, f_1, bv_1)
              call this % Monitoring(1, 's', lambda, nu_1, f_1, bv_1, u_1)
            end associate
          end if

        end associate V_COARSE

        V_UP: do l = 2, l_top

          associate( mesh_l => sem(l  ) % mesh                 &
                   , mesh_p => sem(l-1) % mesh                 &
                   , f_l    => f % level(l  ) % val(:,:,:,:,1) &
                   , u_l    => u % level(l  ) % val(:,:,:,:,1) &
                   , v_l    => r % level(l  ) % val(:,:,:,:,1) &
                   , u_p    => u % level(l-1) % val(:,:,:,:,1) )

            ! set boundary variables for l = t_top .......................

            if (l == l_top) then
              !$omp master
              bv_l => bv % level(l) % var
              !$omp end master
            end if

            ! prolongation .................................................

            call ParentToChildInterpolation_3D &
                     (mesh_p, mesh_l, iop_cf(l-1), u_p, v_l)

            call MergeArrays(ONE, u_l, ONE, v_l)

            ! post-smoothing .................................................

            if (l < l_top .or. m == this % i_max) then
              n = this % ns_2
            else
              n = this % ns_c
            end if

            if (present(nu_0)) then
              call this % Monitoring(l, 'c', lambda, nu_0, f_l, bv_l, u_l)
              call this % Smoother(l, lambda, nu_0, u_l, f_l, bv_l, n)
              call this % Monitoring(l, '2', lambda, nu_0, f_l, bv_l, u_l)
            else
              associate(nu_l => nu_v % level(l) % val(:,:,:,:,1))
                call this % Monitoring(l, 'c', lambda, nu_l, f_l, bv_l, u_l)
                call this % Smoother(l, lambda, nu_l, u_l, f_l, bv_l, n)
                call this % Monitoring(l, '2', lambda, nu_l, f_l, bv_l, u_l)
              end associate
            end if

          end associate
        end do V_UP

        ! termination check ....................................................

        if (check_convergence .and. m < this%i_max) then

          l = l_top
          associate( f_l  => f  % level(l) % val(:,:,:,:,1) &
                   , u_l  => u  % level(l) % val(:,:,:,:,1) &
                   , r_l  => r  % level(l) % val(:,:,:,:,1) )

            if (present(nu_0)) then
              call this % Residual(l, lambda, nu_0, f_l, bv_l, u_l, r_l)
            else
              associate(nu_l => nu_v % level(l) % val(:,:,:,:,1))
                call this % Residual(l, lambda, nu_l, f_l, bv_l, u_l, r_l)
              end associate
            end if

            rr = ScalarProduct(r_l, r_l)
            r_new = sqrt(rr)
          end associate
          !$omp master
          converged = r_new <= r_max
          call XMPI_Bcast(converged, root = 0, comm = sem(1)%mesh%comm_world)
          !$omp end master
          !$omp barrier

          if (converged) exit
          r_old = r_new

        end if

      end do V_OUTER

      ! optional output arguments ::::::::::::::::::::::::::::::::::::::::::::::

      if (present(ni)) then
        !$omp master
        ni = min(m, this%i_max)
        !$omp end master
      end if

      if (present(r_2)) then
        if (check_convergence .and. m < this%i_max) then
          !$omp master
          r_2 = r_new
          !$omp end master
        else
          l = l_top
          associate( f_l  => f  % level(l) % val(:,:,:,:,1) &
                   , u_l  => u  % level(l) % val(:,:,:,:,1) &
                   , r_l  => r  % level(l) % val(:,:,:,:,1) )

            if (present(nu_0)) then
              call this % Residual(l, lambda, nu_0, f_l, bv_l, u_l, r_l)
            else
              associate(nu_l => nu_v % level(l) % val(:,:,:,:,1))
                call this % Residual(l, lambda, nu_l, f_l, bv_l, u_l, r_l)
              end associate
            end if

            rr = ScalarProduct(r_l, r_l)
          end associate
          !$omp master
          r_2 = sqrt(rr)
          !$omp end master
        end if
      end if

      ! finalization :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      !$omp master
      deallocate(r)
      !$omp end master

    end associate

  end subroutine CS_MG_Solver_X

  !=============================================================================

end submodule MP_CS_MG_Solver
