submodule(ML__DG__Elliptic_Solver__3D) MP_FAS_MG_Solver
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> FAS-MG solver for problems with constant diffusivity

  module subroutine FAS_MG_Solver_C(this, lambda, nu, u, f, bv, n_i, r_2)
    class(ML_DG_EllipticSolver_3D), intent(in) :: this
    real(RNP), intent(in) :: lambda                 !< Helmholtz parameter
    real(RNP), intent(in) :: nu                     !< diffusivity
    class(ML_MeshVariable_3D), intent(inout) :: u   !< approx/final solution
    class(ML_MeshVariable_3D), intent(inout) :: f   !< RHS
    class(ML_BoundaryVariable_3D), intent(in) :: bv !< boundary values
    integer,   optional, intent(out) :: n_i         !< num executed cycles
    real(RNP), optional, intent(out) :: r_2         !< Euclidic residual norm

    call FAS_MG_Solver_X(this, lambda, nu, null(), u, f, bv, n_i, r_2)

  end subroutine FAS_MG_Solver_C

  !-----------------------------------------------------------------------------
  !> FAS-MG solver for problems with variable diffusivity

  module subroutine FAS_MG_Solver_V(this, lambda, nu, u, f, bv, n_i, r_2)
    class(ML_DG_EllipticSolver_3D), intent(in) :: this
    real(RNP), intent(in) :: lambda                 !< Helmholtz parameter
    class(ML_MeshVariable_3D), intent(in) :: nu     !< diffusivity
    class(ML_MeshVariable_3D), intent(inout) :: u   !< approx/final solution
    class(ML_MeshVariable_3D), intent(inout) :: f   !< RHS
    class(ML_BoundaryVariable_3D), intent(in) :: bv !< boundary values
    integer, optional, intent(out) :: n_i           !< num executed cycles
    real(RNP), optional, intent(out) :: r_2         !< Euclidic residual norm

    call FAS_MG_Solver_X(this, lambda, null(), nu, u, f, bv, n_i, r_2)

  end subroutine FAS_MG_Solver_V

  !-----------------------------------------------------------------------------
  !> Generic FAS-MG solver for problems with constant or variable diffusivity
  !>
  !> Either `nu_0` or `nu_v` must be given.

  module subroutine FAS_MG_Solver_X(this, lambda, nu_0, nu_v, u, f, bv, n_i, r_2)
    class(ML_DG_EllipticSolver_3D), intent(in) :: this
    real(RNP), intent(in) :: lambda                 !< Helmholtz parameter
    real(RNP), intent(in) :: nu_0                   !< constant diffusivity
    class(ML_MeshVariable_3D), intent(in) :: nu_v   !< variable diffusivity
    class(ML_MeshVariable_3D), intent(inout) :: u   !< approx/final solution
    class(ML_MeshVariable_3D), intent(inout) :: f   !< RHS
    class(ML_BoundaryVariable_3D), intent(in) :: bv !< boundary values
    integer, optional, intent(out) :: n_i           !< num executed cycles
    real(RNP), optional, intent(out) :: r_2         !< Euclidic residual norm

    optional :: nu_0, nu_v

    ! internal variables :::::::::::::::::::::::::::::::::::::::::::::::::::::::

    type(MPI_Comm), save :: comm_world
    type(ML_MeshVariable_3D), allocatable, save :: r, v
    integer, save :: start_method
    integer, save :: l_top
    logical, save :: converged

    real(RNP) :: rr, r_max, r_new, r_old
    logical :: check_convergence
    integer :: e, l, m, n

    associate( sem    => this % ml_op % sem    &
             , iop_cf => this % ml_op % iop_cf &
             , iop_fc => this % ml_op % iop_fc &
             , pop_fc => this % ml_op % pop_fc &
             , ell_op => this % elliptic_op    )

      ! initialization :::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      check_convergence = max(this%r_red, this%r_max) > 0

      !$omp master
      l_top = size(sem)
      r = ML_MeshVariable_3D(this%ml_op, nc=1)
      v = ML_MeshVariable_3D(this%ml_op, nc=1)
      start_method = this % start_method
      comm_world = sem(1) % mesh % comm_world
      !$omp end master

      do l = 1, l_top
        call SetArray(r % level(l) % val(:,:,:,:,1), ZERO)
        call SetArray(v % level(l) % val(:,:,:,:,1), ZERO)
      end do

      ! termination conditions
      if (check_convergence) then
        call this % FAS_MG_Residual_X(lambda, nu_0, nu_v, f, bv, u, r)
        rr = ML_ScalarProduct_3D(r, r)
        r_old  = sqrt(rr)
        r_max  = max(r_old * this%r_red, this%r_max)
        !$omp master
        converged = r_old < r_max
        call XMPI_Bcast(converged, root = 0, comm = comm_world)
        !$omp end master
        !$omp barrier
      else
        !$omp single
        converged = .false.
        !$omp end single
      end if

      if (converged) then
        !$omp master
        deallocate(r, v)
        if (present(n_i)) n_i = 0
        if (present(r_2)) r_2 = r_old
        !$omp end master
        return
      end if

      ! start ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      START: select case(start_method)

      case(START_CASC)

        ! cascade ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        CASC_ROOT: associate( u_1  => u  % level(1) % val(:,:,:,:,1) &
                            , f_1  => f  % level(1) % val(:,:,:,:,1) &
                            , bv_1 => bv % level(1) % var            )

          if (present(nu_0)) then
            call this % CoarseSolver(lambda, nu_0, u_1, f_1, bv_1)
            call this % Monitoring(1, 's', lambda, nu_0, f_1, bv_1, u_1)
          else
            associate(nu_1 => nu_v % level(1) % val(:,:,:,:,1))
              call this % CoarseSolver(lambda, nu_1, u_1, f_1, bv_1)
              call this % Monitoring(1, 's', lambda, nu_1, f_1, bv_1, u_1)
            end associate
          end if

        end associate CASC_ROOT

        CASC_FINE: do l = 2, l_top
          associate( u_p  => u  % level(l-1) % val(:,:,:,:,1) &
                   , u_l  => u  % level(l  ) % val(:,:,:,:,1) &
                   , f_l  => f  % level(l  ) % val(:,:,:,:,1) &
                   , bv_l => bv % level(l  ) % var            )

            ! interpolation
            call ParentToChildInterpolation_3D( parent = sem(l-1) % mesh &
                                              , child  = sem(l  ) % mesh &
                                              , iop    = iop_cf(l-1)     &
                                              , v_p    = u_p             &
                                              , v_c    = u_l             )
            ! smoothing
            n = this % ns_0
            if (l < l_top) then
              if (present(nu_0)) then
                call this % Monitoring(l, 'p', lambda, nu_0, f_l, bv_l, u_l)
                call this % Smoother(l, lambda, nu_0, u_l, f_l, bv_l, n)
                call this % Monitoring(l, '2', lambda, nu_0, f_l, bv_l, u_l)
              else
                associate(nu_l => nu_v % level(l) % val(:,:,:,:,1))
                  call this % Monitoring(l, 'p', lambda, nu_l, f_l, bv_l, u_l)
                  call this % Smoother(l, lambda, nu_l, u_l, f_l, bv_l, n)
                  call this % Monitoring(l, '2', lambda, nu_l, f_l, bv_l, u_l)
                end associate
              end if
            end if

          end associate
        end do CASC_FINE

      case(START_FMG)

        ! FMG ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        FMG_OUTER: do m = 1, l_top-1

          FMG_DOWN: do l = m, 2, -1

            associate( mesh_l => sem(l  ) % mesh                  &
                     , mesh_p => sem(l-1) % mesh                  &
                     , bv_l   => bv % level(l  ) % var            &
                     , f_l    => f  % level(l  ) % val(:,:,:,:,1) &
                     , u_l    => u  % level(l  ) % val(:,:,:,:,1) &
                     , r_l    => r  % level(l  ) % val(:,:,:,:,1) &
                     , f_p    => f  % level(l-1) % val(:,:,:,:,1) &
                     , u_p    => u  % level(l-1) % val(:,:,:,:,1) &
                     , r_p    => r  % level(l-1) % val(:,:,:,:,1) &
                     , v_p    => v  % level(l-1) % val(:,:,:,:,1) )

              ! pre-smoothing and residual computation .........................

              if (l == m) then
                n = this % ns_0
              else
                n = this % ns_1
              end if

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

              ! restriction ....................................................

              ! project solution to regularly refined parent elements
              select case(this % projection_method)
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
              call ChildToParentRestriction_3D &
                       (mesh_l, mesh_p, iop_cf(l-1), r_l, r_p)

              ! parent RHS .....................................................

              do e = 1, mesh_p % n_elem
                if (mesh_p % element(e) % adaptation % refinement >= 1000) then
                  f_p(:,:,:,e) = r_p(:,:,:,e)
                else
                  v_p(:,:,:,e) = u_p(:,:,:,e)
                end if
              end do

              if (present(nu_0)) then
                call ell_op(l-1) % Apply(lambda, nu_0, v_p, r_p)
              else
                associate(nu_p => nu_v % level(l-1) % val(:,:,:,:,1))
                  call ell_op(l-1) % Apply(lambda, nu_p, v_p, r_p)
                end associate
              end if

              do e = 1, mesh_p % n_elem
                if (mesh_p % element(e) % adaptation % refinement >= 1000) then
                  f_p(:,:,:,e) = f_p(:,:,:,e) + r_p(:,:,:,e)
                end if
              end do

            end associate
          end do FMG_DOWN

          FMG_COARSE: associate( bv_1 => bv % level(1) % var            &
                               , f_1  => f  % level(1) % val(:,:,:,:,1) &
                               , u_1  => u  % level(1) % val(:,:,:,:,1) )

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

          end associate FMG_COARSE

          FMG_UP: do l = 2, m

            associate( mesh_l => sem(l  ) % mesh                  &
                     , mesh_p => sem(l-1) % mesh                  &
                     , bv_l   => bv % level(l  ) % var            &
                     , f_l    => f  % level(l  ) % val(:,:,:,:,1) &
                     , u_l    => u  % level(l  ) % val(:,:,:,:,1) &
                     , w_l    => r  % level(l  ) % val(:,:,:,:,1) &
                     , u_p    => u  % level(l-1) % val(:,:,:,:,1) &
                     , v_p    => v  % level(l-1) % val(:,:,:,:,1) )

                ! prolongation .................................................

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
                if (mesh_p % n_elem > 0) then
                  n = mesh_p % n_elem_active
                  ! apply correction to active elements
                  call MergeArrays(ONE, u_l(:,:,:,:n), ONE, w_l(:,:,:,:n))
                  ! update frozen elements
                  if (mesh_p % n_elem_frozen > 0) then
                    call SetArray(u_l(:,:,:,n+1:), w_l(:,:,:,n+1:))
                  end if
                end if

              ! post-smoothing .................................................

              n = this % ns_2
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
          end do FMG_UP

          ! interpolation to next level
          call ParentToChildInterpolation_3D( parent = sem(m  ) % mesh      &
                                            , child  = sem(m+1) % mesh      &
                                            , iop    = iop_cf(m)            &
                                            , v_p    = u % level(m  ) % val &
                                            , v_c    = u % level(m+1) % val )

        end do FMG_OUTER

      end select START

      ! V cycles :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      V_OUTER: do m = 1, this % i_max

        V_DOWN: do l = l_top, 2, -1

          associate( mesh_l => sem(l  ) % mesh                  &
                   , mesh_p => sem(l-1) % mesh                  &
                   , bv_l   => bv % level(l  ) % var            &
                   , f_l    => f  % level(l  ) % val(:,:,:,:,1) &
                   , u_l    => u  % level(l  ) % val(:,:,:,:,1) &
                   , r_l    => r  % level(l  ) % val(:,:,:,:,1) &
                   , f_p    => f  % level(l-1) % val(:,:,:,:,1) &
                   , u_p    => u  % level(l-1) % val(:,:,:,:,1) &
                   , r_p    => r  % level(l-1) % val(:,:,:,:,1) &
                   , v_p    => v  % level(l-1) % val(:,:,:,:,1) )

            ! pre-smoothing and residual computation .........................

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

            ! restriction ....................................................

            ! project solution to regularly refined parent elements
            select case(this % projection_method)
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
            call ChildToParentRestriction_3D &
                     (mesh_l, mesh_p, iop_cf(l-1), r_l, r_p)

            ! parent RHS .....................................................

            do e = 1, mesh_p % n_elem
              if (mesh_p % element(e) % adaptation % refinement >= 1000) then
                f_p(:,:,:,e) = r_p(:,:,:,e)
              else
                v_p(:,:,:,e) = u_p(:,:,:,e)
              end if
            end do

            if (present(nu_0)) then
              call ell_op(l-1) % Apply(lambda, nu_0, v_p, r_p)
            else
              associate(nu_p => nu_v % level(l-1) % val(:,:,:,:,1))
                call ell_op(l-1) % Apply(lambda, nu_p, v_p, r_p)
              end associate
            end if

             do e = 1, mesh_p % n_elem
               if (mesh_p % element(e) % adaptation % refinement >= 1000) then
                 f_p(:,:,:,e) = f_p(:,:,:,e) + r_p(:,:,:,e)
               end if
             end do

          end associate
        end do V_DOWN

        V_COARSE: associate( bv_1 => bv % level(1) % var            &
                           , f_1  => f  % level(1) % val(:,:,:,:,1) &
                           , u_1  => u  % level(1) % val(:,:,:,:,1) )

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

          associate( mesh_l => sem(l  ) % mesh                    &
                   , mesh_p => sem(l-1) % mesh                    &
                   , bv_l   => bv   % level(l  ) % var            &
                   , f_l    => f    % level(l  ) % val(:,:,:,:,1) &
                   , u_l    => u    % level(l  ) % val(:,:,:,:,1) &
                   , w_l    => r    % level(l  ) % val(:,:,:,:,1) &
                   , u_p    => u    % level(l-1) % val(:,:,:,:,1) &
                   , v_p    => v    % level(l-1) % val(:,:,:,:,1) )

              ! prolongation .................................................

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

          call this % FAS_MG_Residual_X(lambda, nu_0, nu_v, f, bv, u, r)
          rr = ML_ScalarProduct_3D(r, r)
          r_new = sqrt(rr)

          !$omp master
          converged = r_new <= r_max
          call XMPI_Bcast(converged, root = 0, comm = comm_world)
          !$omp end master
          !$omp barrier

          if (converged) exit
          r_old = r_new

        end if

      end do V_OUTER

      ! optional output arguments ::::::::::::::::::::::::::::::::::::::::::::::

      if (present(n_i)) then
        !$omp master
        n_i = min(m, this%i_max)
        !$omp end master
      end if

      if (present(r_2)) then
        if (check_convergence .and. m < this%i_max) then
          !$omp master
          r_2 = r_new
          !$omp end master
        else
          call this % FAS_MG_Residual_X(lambda, nu_0, nu_v, f, bv, u, r)
          rr = ML_ScalarProduct_3D(r, r)
          !$omp master
          r_2 = sqrt(rr)
          !$omp end master
        end if
      end if

      ! finalization :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      !$omp master
      deallocate(r, v)
      !$omp end master

    end associate

  end subroutine FAS_MG_Solver_X

  !=============================================================================

end submodule MP_FAS_MG_Solver
