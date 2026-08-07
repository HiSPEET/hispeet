!> summary:  FAS-MG V-cycle for multilevel elliptic solver
!> author:   Joerg Stiller
!> date:     2025/09/04
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(ML__DG__Elliptic_Solver__3D) MP_FAS_MG_Cycle
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Generic FAS-MG V-cycle for problems with constant or variable diffusivity
  !>
  !> Either `nu_0` or `nu_v` must be given.

  module subroutine FAS_MG_Cycle_X( this, bc, lambda, nu_0, nu_v, bv, f, u &
                                  , n_cyc, l_top, r0_2, ni, r_2 )

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
      !< boundary values
    class(ML_MeshVariable_3D), intent(in) :: f
      !< RHS
    class(ML_MeshVariable_3D), intent(inout) :: u
      !< approx/final solution
    integer, optional, intent(in) :: n_cyc
      !< number of cycles  [this%i_max]
    integer, optional, intent(in) :: l_top
      !< top level  [auto]
    real(RNP), optional, intent(in) :: r0_2
      !< initial Euclidian residual norm, if < 0
    integer, optional, intent(out) :: ni
      !< num executed cycles
    real(RNP), optional, intent(out) :: r_2
      !< Euclidean residual norm

    ! internal variables :::::::::::::::::::::::::::::::::::::::::::::::::::::::

    type(ML_MeshVariable_3D), allocatable, save :: g, r, v
    type(BoundaryVariable_3D), pointer, save :: bv_l(:) => null()
    type(BoundaryVariable_3D), pointer, save :: bv_p(:) => null()
    logical, save :: converged

    real(RNP) :: rr, r_max, r_new, r_old
    logical   :: check_convergence
    integer   :: l_top_, n_cyc_
    integer   :: e, l, m, n

    character(len=:), allocatable :: prefix
    logical :: logging

    ! start logging ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    if (log_level > 0) then
      associate(proc => this % ml_op % sem(1) % mesh % proc)
        logging = proc == 0 .or. log_level > 1
        prefix  = LoggingPrefix('FAS_MG_Cycle_X', proc)
      end associate
    else
      logging = .false.
    end if

    if (logging) then
      print '(2A)', prefix, 'start'
    end if

    associate( sem    => this % ml_op % sem      &
             , iop_cf => this % ml_op % iop_cf_x &
             , iop_fc => this % ml_op % iop_fc_x &
             , pop_fc => this % ml_op % pop_fc_x &
             , ell    => this % elliptic_op      )

      ! initialization :::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      if (present(n_cyc)) then
        n_cyc_ = n_cyc
        check_convergence = .false.
      else
        n_cyc_ = this % i_max
        check_convergence = max(this%r_red, this%r_max) > 0
      end if
      if (n_cyc_ < 1) return

      if (present(l_top)) then
        l_top_ = min(l_top, size(sem))
      else
        l_top_ = size(sem)
      end if

      !$omp master
      allocate(g, r, v)
      call g % Init(this%ml_op, nc = 1, l_top = l_top_)
      call r % Init(this%ml_op, nc = 1, l_top = l_top_)
      call v % Init(this%ml_op, nc = 1, l_top = l_top_)
      !$omp end master
      !$omp barrier

      call ML_SetArray_3D(g, f, l_top = l_top_)

      ! termination condition
      if (check_convergence) then
        if (present(r0_2)) then
          r_old = r0_2
        else
          call this % FAS_MG_Residual_X( bc, lambda, nu_0, nu_v, bv, f, u, r &
                                       , l_top_ )
          rr = ML_ScalarProduct_3D(r, r, l_top = l_top_)
          r_old = sqrt(rr)
        end if
        r_max = max(r_old * this%r_red, this%r_max)
        !$omp master
        converged = r_old < r_max
        call XMPI_Bcast(converged, root = 0, comm = sem(1)%mesh%comm_world)
        !$omp end master
      else
        !$omp master
        converged = .false.
        !$omp end master
      end if
      !$omp barrier

      if (converged) then
        !$omp master
        if (present(ni)) ni = 0
        if (present(r_2)) r_2 = r_old
        deallocate(g, r, v)
        !$omp end master
        return
      end if

      ! V cycles :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      V_OUTER: do m = 1, n_cyc_

        V_DOWN: do l = l_top_, 2, -1

          associate( mesh_l => sem(l  ) % mesh                  &
                   , mesh_p => sem(l-1) % mesh                  &
                   , g_l    => g  % level(l  ) % val(:,:,:,:,1) &
                   , u_l    => u  % level(l  ) % val(:,:,:,:,1) &
                   , r_l    => r  % level(l  ) % val(:,:,:,:,1) &
                   , g_p    => g  % level(l-1) % val(:,:,:,:,1) &
                   , u_p    => u  % level(l-1) % val(:,:,:,:,1) &
                   , r_p    => r  % level(l-1) % val(:,:,:,:,1) &
                   , v_p    => v  % level(l-1) % val(:,:,:,:,1) )

            ! boundary values ..................................................

            !$omp master
            if (present(bv)) then
              bv_l => bv % level(l ) % var
              bv_p => bv % level(l-1) % var
            end if
            !$omp end master

            ! pre-smoothing and residual computation ...........................

            if (present(nu_0)) then
              call this % Monitoring(l, '0', bc, lambda, nu_0, g_l, bv_l, u_l)
            else
              associate(nu_l => nu_v % level(l) % val(:,:,:,:,1))
                call this % Monitoring(l, '0', bc, lambda, nu_l, g_l, bv_l, u_l)
              end associate
            end if

            if (l < l_top_ .or. m == 1) then
              n = this % NumSmoothingSteps(l, stage = 1)
            else
              n = 0
            end if

            if (present(nu_0)) then
              call this % Smoother(l, bc, lambda, nu_0, u_l, g_l, bv_l, n)
              call this % Residual(l, bc, lambda, nu_0, g_l, bv_l, u_l, r_l)
            else
              associate(nu_l => nu_v % level(l) % val(:,:,:,:,1))
                call this % Smoother(l, bc, lambda, nu_l, u_l, g_l, bv_l, n)
                call this % Residual(l, bc, lambda, nu_l, g_l, bv_l, u_l, r_l)
              end associate
            end if
            call this % Monitoring(l, '1', r_l)

            ! restriction ....................................................

            ! project solution to regularly refined parent elements
            select case(this % fc_projection)
            case('I')
              ! interpolation
              call ChildToParentProjection_3D &
                       (mesh_l, mesh_p, iop_fc(l), u_l, v_p)
            case('P')
              ! L²-projection
              call ChildToParentProjection_3D &
                       (mesh_l, mesh_p, pop_fc(l), u_l, v_p)
            end select

            ! restrict residual
            select case(this % fc_restriction)
            case('C')
              ! canonical restriction
              call ChildToParentRestriction_3D &
                       (mesh_l, mesh_p, iop_cf(l-1), r_l, r_p)
            case('P')
              ! L²-projection
              call ChildToParentProjection_3D &
                       (mesh_l, mesh_p, pop_fc(l), r_l, r_p)
            end select

            ! parent FAS-RHS .................................................

            do e = 1, mesh_p % n_elem
              if (mesh_p % element(e) % adaptation % refinement >= 1000) then
                ! residual contribution to FAS-RHS in parent twigs
                g_p(:,:,:,e) = r_p(:,:,:,e)
              else
              ! set v_p to solution in leaves
                v_p(:,:,:,e) = u_p(:,:,:,e)
              end if
            end do

            ! apply parent operator to projected solution
            if (present(nu_0)) then
              call ell(l-1) % Apply(bc, lambda, nu_0, bv_p, v_p, r_p)
            else
              associate(nu_p => nu_v % level(l-1) % val(:,:,:,:,1))
                call ell(l-1) % Apply(bc, lambda, nu_p, bv_p, v_p, r_p)
              end associate
            end if

            do e = 1, mesh_p % n_elem
              if (mesh_p % element(e) % adaptation % refinement >= 1000) then
                g_p(:,:,:,e) = g_p(:,:,:,e) + r_p(:,:,:,e)
              end if
            end do

          end associate
        end do V_DOWN

        associate( g_l  => g % level(1) % val(:,:,:,:,1) &
                 , u_l  => u % level(1) % val(:,:,:,:,1) )

          ! coarse grid solver .................................................

          !$omp master
          if (present(bv)) then
            bv_l => bv % level(1) % var
          end if
          !$omp end master

          if (present(nu_0)) then
            call this % Monitoring(1, '0', bc, lambda, nu_0, g_l, bv_l, u_l)
            call this % CoarseSolver(bc, lambda, nu_0, u_l, g_l, bv_l)
            call this % Monitoring(1, 's', bc, lambda, nu_0, g_l, bv_l, u_l)
          else
            associate(nu_l => nu_v % level(1) % val(:,:,:,:,1))
              call this % Monitoring(1, '0', bc, lambda, nu_l, g_l, bv_l, u_l)
              call this % CoarseSolver(bc, lambda, nu_l, u_l, g_l, bv_l)
              call this % Monitoring(1, 's', bc, lambda, nu_l, g_l, bv_l, u_l)
            end associate
          end if

        end associate

        V_UP: do l = 2, l_top_

          associate( mesh_l => sem(l  ) % mesh                 &
                   , mesh_p => sem(l-1) % mesh                 &
                   , g_l    => g % level(l  ) % val(:,:,:,:,1) &
                   , u_l    => u % level(l  ) % val(:,:,:,:,1) &
                   , w_l    => r % level(l  ) % val(:,:,:,:,1) &
                   , u_p    => u % level(l-1) % val(:,:,:,:,1) &
                   , v_p    => v % level(l-1) % val(:,:,:,:,1) )

            ! boundary values ..................................................

            !$omp master
            if (present(bv)) then
              bv_l => bv % level(l) % var
            end if
            !$omp end master

            ! prolongation .....................................................

            do e = 1, mesh_p % n_elem
              if (mesh_p % element(e) % adaptation % refinement >= 1000) then
                v_p(:,:,:,e) = u_p(:,:,:,e) - v_p(:,:,:,e)
              else
                v_p(:,:,:,e) = u_p(:,:,:,e)
              end if
            end do

            call ParentToChildInterpolation_3D &
                     (mesh_p, mesh_l, iop_cf(l-1), v_p, w_l)

            ! update solution on current level
            if (mesh_l % n_elem > 0) then
              n = mesh_l % n_elem_active
              ! apply correction to active elements
              call MergeArrays(ONE, u_l(:,:,:,:n), ONE, w_l(:,:,:,:n))
              ! update frozen elements
              if (mesh_l % n_elem_frozen > 0) then
                call SetArray(u_l(:,:,:,n+1:), w_l(:,:,:,n+1:))
              end if
            end if

            ! post-smoothing ...................................................

            if (l < l_top_) then
              n = this % NumSmoothingSteps(l, stage = 2)
            else if (m < n_cyc_) then
              n = this % ns_c
            else
              n = this % ns_f
            end if

            if (present(nu_0)) then
              call this % Monitoring(l, 'c', bc, lambda, nu_0, g_l, bv_l, u_l)
              call this % Smoother(l, bc, lambda, nu_0, u_l, g_l, bv_l, n)
              call this % Monitoring(l, '2', bc, lambda, nu_0, g_l, bv_l, u_l)
            else
              associate(nu_l => nu_v % level(l) % val(:,:,:,:,1))
                call this % Monitoring(l, 'c', bc, lambda, nu_l, g_l, bv_l, u_l)
                call this % Smoother(l, bc, lambda, nu_l, u_l, g_l, bv_l, n)
                call this % Monitoring(l, '2', bc, lambda, nu_l, g_l, bv_l, u_l)
              end associate
            end if

          end associate
        end do V_UP

        ! termination check ....................................................

        if (check_convergence .and. m < this%i_max) then

          call this % FAS_MG_Residual_X( bc, lambda, nu_0, nu_v, bv, f, u, r &
                                       , l_top_ )
          rr = ML_ScalarProduct_3D(r, r, l_top = l_top_)
          r_new = sqrt(rr)

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
          call this % FAS_MG_Residual_X( bc, lambda, nu_0, nu_v, bv, f &
                                       , u, r, l_top_)
          rr = ML_ScalarProduct_3D(r, r, l_top = l_top_)
          !$omp master
          r_2 = sqrt(rr)
          !$omp end master
        end if
      end if

      ! finalization :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      !$omp master
      deallocate(g, r, v)
      bv_l => null()
      bv_p => null()
      !$omp end master

    end associate

    ! exit logging :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    if (logging) then
      print '(2A)', prefix, 'exit'
    end if

  end subroutine FAS_MG_Cycle_X

  !=============================================================================

end submodule MP_FAS_MG_Cycle
