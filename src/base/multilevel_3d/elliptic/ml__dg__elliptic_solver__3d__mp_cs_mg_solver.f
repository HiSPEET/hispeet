submodule(ML__DG__Elliptic_Solver__3D) MP_CS_MG_Solver
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> CS-MG solver for problems with global refinement and constant diffusivity
  !>
  !> Use `l_top` to specify a top level lower than `size(this%ml_op%sem)`

  module subroutine CS_MG_Solver_C( this, bc, lambda, nu, bv, f, u &
                                  , i_max, l_top, ni, r_2          )
    class(ML_DG_EllipticSolver_3D), intent(in) :: this
    character, intent(in) :: bc(:)                  !< boundary conditions
    real(RNP), intent(in) :: lambda                 !< Helmholtz parameter
    real(RNP), intent(in) :: nu                     !< diffusivity
    class(ML_BoundaryVariable_3D), intent(in) :: bv !< boundary values
    class(ML_MeshVariable_3D), intent(inout) :: f   !< RHS
    class(ML_MeshVariable_3D), intent(inout) :: u   !< approx/final solution
    integer,   optional, intent(in)  :: i_max       !< overrides max num cycles
    integer,   optional, intent(in)  :: l_top       !< top level
    integer,   optional, intent(out) :: ni          !< num executed cycles
    real(RNP), optional, intent(out) :: r_2         !< Euclidean residual norm

    call CS_MG_Solver_X &
             (this, bc, lambda, nu, null(), bv, f, u, i_max, l_top, ni, r_2)

  end subroutine CS_MG_Solver_C

  !-----------------------------------------------------------------------------
  !> CS-MG solver for problems with global refinement and variable diffusivity
  !>
  !> Use `l_top` to specify a top level lower than `size(this%ml_op%sem)`

  module subroutine CS_MG_Solver_V( this, bc, lambda, nu, bv, f, u &
                                  , i_max, l_top, ni, r_2          )
    class(ML_DG_EllipticSolver_3D), intent(in) :: this
    character, intent(in) :: bc(:)                  !< boundary conditions
    real(RNP), intent(in) :: lambda                 !< Helmholtz parameter
    class(ML_MeshVariable_3D), intent(in) :: nu     !< diffusivity
    class(ML_BoundaryVariable_3D), intent(in) :: bv !< boundary values
    class(ML_MeshVariable_3D), intent(inout) :: f   !< RHS
    class(ML_MeshVariable_3D), intent(inout) :: u   !< approx/final solution
    integer,   optional, intent(in)  :: i_max       !< overrides max num cycles
    integer,   optional, intent(in)  :: l_top       !< top level
    integer,   optional, intent(out) :: ni          !< num executed cycles
    real(RNP), optional, intent(out) :: r_2         !< Euclidean residual norm

    call CS_MG_Solver_X &
             (this, bc, lambda, null(), nu, bv, f, u, i_max, l_top, ni, r_2)

  end subroutine CS_MG_Solver_V

  !-----------------------------------------------------------------------------
  !> Generic CS-MG solver for problems with constant or variable diffusivity
  !>
  !> Either `nu_0` or `nu_v` must be given.

  module subroutine CS_MG_Solver_X( this, bc, lambda, nu_0, nu_v, bv, f, u &
                                  , i_max, l_top, ni, r_2 )

    class(ML_DG_EllipticSolver_3D), intent(in) :: this

    character, intent(in) :: bc(:)
      !< boundary conditions
    real(RNP), intent(in) :: lambda
      !< Helmholtz parameter
    real(RNP), optional, intent(in) :: nu_0
      !< constant diffusivity
    class(ML_MeshVariable_3D), optional, intent(in) :: nu_v
      !< variable diffusivity
    class(ML_BoundaryVariable_3D), target, optional, intent(in) :: bv
      !< boundary values [homogeneous]
    class(ML_MeshVariable_3D), intent(inout) :: f
      !< RHS
    class(ML_MeshVariable_3D), intent(inout) :: u
      !< approx/final solution
    integer, optional, intent(in) :: i_max
      !< overrides preset maximum number of cycles
    integer, optional, intent(in) :: l_top
      !< top level different from size(this%ml_op%sem)
    integer, optional, intent(out) :: ni
      !< number of executed cycles
    real(RNP), optional, intent(out) :: r_2
      !< Euclidean residual norm

    ! internal variables :::::::::::::::::::::::::::::::::::::::::::::::::::::::

    type(ML_MeshVariable_3D), allocatable, save :: r
    type(BoundaryVariable_3D), pointer, save :: bv_l(:)
    logical, save :: converged

    real(RNP) :: rr, r_max, r_new, r_old
    logical :: check_convergence
    integer :: i_max_, l_top_
    integer :: l, m, n

    ! prerequisites ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    if (present(l_top)) then
      l_top_ = min(l_top, size(this % ml_op % sem))
    else
      l_top_ = size(this % ml_op % sem)
    end if

    if (present(i_max)) then
      i_max_ = i_max
    else
      i_max_ = this % i_max
    end if

    check_convergence = i_max_ > 1 .and. max(this%r_red, this%r_max) > 0

    associate( sem    => this % ml_op % sem      &
             , iop_cf => this % ml_op % iop_cf_x &
             , iop_fc => this % ml_op % iop_fc_x &
             , pop_fc => this % ml_op % pop_fc_x )

      ! initialization :::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      !$omp master
      allocate(r)
      call r % Init(this%ml_op, nc=1)
      !$omp end master

      ! termination conditions
      if (check_convergence) then
        l = l_top_
        associate( mesh_l => sem(l) % mesh              &
                 , f_l => f % level(l) % val(:,:,:,:,1) &
                 , u_l => u % level(l) % val(:,:,:,:,1) &
                 , r_l => r % level(l) % val(:,:,:,:,1) )

          !$omp master
          if (present(bv)) then
            bv_l => bv % level(l) % var
          end if
          !$omp end master

          ! initial residual
          if (present(nu_0)) then
            call this % Residual(l, bc, lambda, nu_0, f_l, bv_l, u_l, r_l)
          else
            associate(nu_l => nu_v % level(l) % val(:,:,:,:,1))
              call this % Residual(l, bc, lambda, nu_l, f_l, bv_l, u_l, r_l)
            end associate
          end if

          rr = ScalarProduct(r_l, r_l, mesh_l%comm_parts)
          r_old  = sqrt(rr)
          r_max  = max(r_old * this%r_red, this%r_max)
          !$omp master
          converged = r_old < r_max
          call XMPI_Bcast(converged, root = 0, comm = mesh_l%comm_parts)
          !$omp end master
        end associate
      else
        !$omp master
        converged = .false.
        !$omp end master
      end if
      !$omp barrier

      if (converged) then
        !$omp master
        deallocate(r)
        if (present(ni)) ni = 0
        if (present(r_2)) r_2 = r_old
        !$omp end master
        return
      end if

      ! V cycles :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      V_OUTER: do m = 1, i_max_

        V_DOWN: do l = l_top_, 2, -1

          associate( mesh_l => sem(l  ) % mesh                 &
                   , mesh_p => sem(l-1) % mesh                 &
                   , f_l    => f % level(l  ) % val(:,:,:,:,1) &
                   , u_l    => u % level(l  ) % val(:,:,:,:,1) &
                   , r_l    => r % level(l  ) % val(:,:,:,:,1) &
                   , f_p    => f % level(l-1) % val(:,:,:,:,1) )

            !$omp master
            if (l == l_top_ .and. present(bv)) then
              bv_l => bv % level(l) % var
            else
              bv_l => null()
            end if
            !$omp end master

            ! pre-smoothing and residual computation .........................

            if (l < l_top_) then
              call SetArray(u_l, ZERO)
            end if

            if (present(nu_0)) then
              call this % Monitoring(l, '0', bc, lambda, nu_0, f_l, bv_l, u_l)
            else
              associate(nu_l => nu_v % level(l) % val(:,:,:,:,1))
                call this % Monitoring(l, '0', bc, lambda, nu_l, f_l, bv_l, u_l)
              end associate
            end if

            n = this % NumSmoothingSteps(l, 1)
            if (present(nu_0)) then
              call this % Smoother(l, bc, lambda, nu_0, u_l, f_l, bv_l, n)
              call this % Residual(l, bc, lambda, nu_0, f_l, bv_l, u_l, r_l)
            else
              associate(nu_l => nu_v % level(l) % val(:,:,:,:,1))
                call this % Smoother(l, bc, lambda, nu_l, u_l, f_l, bv_l, n)
                call this % Residual(l, bc, lambda, nu_l, f_l, bv_l, u_l, r_l)
              end associate
            end if
            call this % Monitoring(l, '1', r_l)

            ! restriction ....................................................

            select case(this % fc_restriction)
            case('C')
              ! canonical restriction
              call ChildToParentRestriction_3D &
                       (mesh_l, mesh_p, iop_cf(l-1), r_l, f_p)
            case('P')
              ! L²-projection
              call ChildToParentProjection_3D &
                       (mesh_l, mesh_p, pop_fc(l), r_l, f_p)
            end select

          end associate
        end do V_DOWN

        associate( f_l => f % level(1) % val(:,:,:,:,1) &
                 , u_l => u % level(1) % val(:,:,:,:,1) )

          ! coarse grid solver ...............................................

          !$omp master
          if (l_top_ == 1 .and. present(bv)) then
            bv_l => bv % level(1) % var
          else
            bv_l => null()
          end if
          !$omp end master

          if (present(nu_0)) then
            call this % Monitoring(1, '0', bc, lambda, nu_0, f_l, bv_l, u_l)
            call this % CoarseSolver(bc, lambda, nu_0, u_l, f_l, bv_l)
            call this % Monitoring(1, 's', bc, lambda, nu_0, f_l, bv_l, u_l)
          else
            associate(nu_l => nu_v % level(1) % val(:,:,:,:,1))
              call this % Monitoring(1, '0', bc, lambda, nu_l, f_l, bv_l, u_l)
              call this % CoarseSolver(bc, lambda, nu_l, u_l, f_l, bv_l)
              call this % Monitoring(1, 's', bc, lambda, nu_l, f_l, bv_l, u_l)
            end associate
          end if

        end associate

        V_UP: do l = 2, l_top_

          associate( mesh_l => sem(l  ) % mesh                 &
                   , mesh_p => sem(l-1) % mesh                 &
                   , f_l    => f % level(l  ) % val(:,:,:,:,1) &
                   , u_l    => u % level(l  ) % val(:,:,:,:,1) &
                   , v_l    => r % level(l  ) % val(:,:,:,:,1) &
                   , u_p    => u % level(l-1) % val(:,:,:,:,1) )

            ! boundary values ..................................................

            !$omp master
            if (l == l_top_ .and. present(bv)) then
              bv_l => bv % level(l) % var
            else
              bv_l => null()
            end if
            !$omp end master

            ! prolongation .....................................................

            call ParentToChildInterpolation_3D &
                     (mesh_p, mesh_l, iop_cf(l-1), u_p, v_l)

            call MergeArrays(ONE, u_l, ONE, v_l)

            ! post-smoothing ...................................................

            if (l < l_top_) then
              n = this % NumSmoothingSteps(l, 2)
            else if (m < i_max_) then
              n = this % ns_c
            else
              n = this % ns_f
            end if

            if (present(nu_0)) then
              call this % Monitoring(l, 'c', bc, lambda, nu_0, f_l, bv_l, u_l)
              call this % Smoother(l, bc, lambda, nu_0, u_l, f_l, bv_l, n)
              call this % Monitoring(l, '2', bc, lambda, nu_0, f_l, bv_l, u_l)
            else
              associate(nu_l => nu_v % level(l) % val(:,:,:,:,1))
                call this % Monitoring(l, 'c', bc, lambda, nu_l, f_l, bv_l, u_l)
                call this % Smoother(l, bc, lambda, nu_l, u_l, f_l, bv_l, n)
                call this % Monitoring(l, '2', bc, lambda, nu_l, f_l, bv_l, u_l)
              end associate
            end if

          end associate
        end do V_UP

        ! termination check ....................................................

        if (check_convergence .and. m < i_max_) then

          l = l_top_
          associate( f_l  => f % level(l) % val(:,:,:,:,1) &
                   , u_l  => u % level(l) % val(:,:,:,:,1) &
                   , r_l  => r % level(l) % val(:,:,:,:,1) )

            !$omp master
            if (present(bv)) then
              bv_l => bv % level(l) % var
            end if
            !$omp end master

            if (present(nu_0)) then
              call this % Residual(l, bc, lambda, nu_0, f_l, bv_l, u_l, r_l)
            else
              associate(nu_l => nu_v % level(l) % val(:,:,:,:,1))
                call this % Residual(l, bc, lambda, nu_l, f_l, bv_l, u_l, r_l)
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
        ni = min(m, i_max_)
        !$omp end master
      end if

      if (present(r_2)) then
        if (check_convergence .and. m < i_max_) then
          !$omp master
          r_2 = r_new
          !$omp end master
        else
          l = l_top_
          associate( f_l  => f  % level(l) % val(:,:,:,:,1) &
                   , u_l  => u  % level(l) % val(:,:,:,:,1) &
                   , r_l  => r  % level(l) % val(:,:,:,:,1) )

            !$omp master
            if (present(bv)) then
              bv_l => bv % level(l) % var
            end if
            !$omp end master

            if (present(nu_0)) then
              call this % Residual(l, bc, lambda, nu_0, f_l, bv_l, u_l, r_l)
            else
              associate(nu_l => nu_v % level(l) % val(:,:,:,:,1))
                call this % Residual(l, bc, lambda, nu_l, f_l, bv_l, u_l, r_l)
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
