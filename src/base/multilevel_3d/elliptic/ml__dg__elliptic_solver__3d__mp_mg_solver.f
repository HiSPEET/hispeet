submodule(ML__DG__Elliptic_Solver__3D) MP_MG_Solver
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> FAS-MG solver for problems with constant diffusivity

  module subroutine MG_Solver_C(this, lambda, nu, u, f, bv, n_i, r_2)
    class(ML_DG_EllipticSolver_3D), intent(in) :: this
    real(RNP), intent(in) :: lambda
    real(RNP), intent(in) :: nu
    class(ML_MeshVariable_3D), intent(inout) :: u
    class(ML_MeshVariable_3D), intent(inout) :: f
    class(ML_BoundaryVariable_3D), intent(in) :: bv
    integer, optional, intent(out) :: n_i
    real(RNP), optional, intent(out) :: r_2(:)

    call MG_Solver_X(this, lambda, nu, null(), u, f, bv, n_i, r_2)

  end subroutine MG_Solver_C

  !-----------------------------------------------------------------------------
  !> FAS-MG solver for problems with variable diffusivity

  module subroutine MG_Solver_V(this, lambda, nu, u, f, bv, n_i, r_2)
    class(ML_DG_EllipticSolver_3D), intent(in) :: this
    real(RNP), intent(in) :: lambda
    class(ML_MeshVariable_3D), intent(in) :: nu
    class(ML_MeshVariable_3D), intent(inout) :: u
    class(ML_MeshVariable_3D), intent(inout) :: f
    class(ML_BoundaryVariable_3D), intent(in) :: bv
    integer, optional, intent(out) :: n_i
    real(RNP), optional, intent(out) :: r_2(:)

    call MG_Solver_X(this, lambda, null(), nu, u, f, bv, n_i, r_2)

  end subroutine MG_Solver_V

  !-----------------------------------------------------------------------------
  !> Generic FAS-MG solver for problems with constant or variable diffusivity

  subroutine MG_Solver_X(this, lambda, nu_0, nu_v, u, f, bv, n_i, r_2)
    class(ML_DG_EllipticSolver_3D), intent(in) :: this
    real(RNP), intent(in) :: lambda
    real(RNP), optional, intent(in) :: nu_0
    class(ML_MeshVariable_3D), optional, intent(in) :: nu_v
    class(ML_MeshVariable_3D), intent(inout) :: u
    class(ML_MeshVariable_3D), intent(inout) :: f
    class(ML_BoundaryVariable_3D), intent(in) :: bv
    integer, optional, intent(out) :: n_i
    real(RNP), optional, intent(out) :: r_2(:)

    type(ML_MeshVariable_3D), allocatable, save :: r, v
    integer, save :: start_method
    integer, save :: l_top

    integer :: e, l, m, n

    associate( sem    => this % ml_op % sem    &
             , iop_cf => this % ml_op % iop_cf &
             , iop_fc => this % ml_op % iop_fc &
             , pop_fc => this % ml_op % pop_fc &
             , ell_op => this % elliptic_op    )

      ! initialization :::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      !$omp master
      l_top = size(sem)
      r = ML_MeshVariable_3D(this%ml_op, nc=1)
      v = ML_MeshVariable_3D(this%ml_op, nc=1)
      start_method = this % start_method
      !$omp end master

      do l = 1, l_top
        call SetArray(r % level(l) % val(:,:,:,:,1), ZERO)
        call SetArray(v % level(l) % val(:,:,:,:,1), ZERO)
      end do

      ! start ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      START: select case(start_method)

      case(START_CASC)

        ! cascade ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        CASC_ROOT: associate( u_1  => u    % level(1) % val(:,:,:,:,1) &
                            , f_1  => f    % level(1) % val(:,:,:,:,1) &
                            , nu_1 => nu_v % level(1) % val(:,:,:,:,1) &
                            , bv_1 => bv   % level(1) % var            )

          if (present(nu_0)) then
            call this % CoarseSolver(lambda, nu_0, u_1, f_1, bv_1)
                call this % Monitoring(1, 's', lambda, nu_0, f_1, bv_1, u_1)
          else
            call this % CoarseSolver(lambda, nu_1, u_1, f_1, bv_1)
                call this % Monitoring(1, 's', lambda, nu_1, f_1, bv_1, u_1)
          end if

        end associate CASC_ROOT

        CASC_FINE: do l = 2, l_top
          associate( u_p  => u    % level(l-1) % val(:,:,:,:,1) &
                   , u_l  => u    % level(l  ) % val(:,:,:,:,1) &
                   , f_l  => f    % level(l  ) % val(:,:,:,:,1) &
                   , nu_l => nu_v % level(l  ) % val(:,:,:,:,1) &
                   , bv_l => bv   % level(l  ) % var            )

            ! interpolation
            call ParentToChildInterpolation_3D( parent = sem(l-1) % mesh &
                                              , child  = sem(l  ) % mesh &
                                              , iop    = iop_cf(l-1)     &
                                              , v_p    = u_p             &
                                              , v_c    = u_l             )
            ! smoothing
            n = this % ns_2
            if (l < l_top) then
              if (present(nu_0)) then
                call this % Monitoring(l, 'p', lambda, nu_0, f_l, bv_l, u_l)
                call this % Smoother(l, lambda, nu_0, u_l, f_l, bv_l, n)
                call this % Monitoring(l, '2', lambda, nu_0, f_l, bv_l, u_l)
              else
                call this % Monitoring(l, 'p', lambda, nu_l, f_l, bv_l, u_l)
                call this % Smoother(l, lambda, nu_l, u_l, f_l, bv_l, n)
                call this % Monitoring(l, '2', lambda, nu_l, f_l, bv_l, u_l)
              end if
            end if

          end associate
        end do CASC_FINE

      case(START_FMG)

        ! FMG ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        FMG_OUTER: do m = 1, l_top-1

          FMG_DOWN: do l = m, 2, -1

            associate( mesh_l => sem(l  ) % mesh                    &
                     , mesh_p => sem(l-1) % mesh                    &
                     , bv_l   => bv   % level(l  ) % var            &
                     , nu_l   => nu_v % level(l  ) % val(:,:,:,:,1) &
                     , f_l    => f    % level(l  ) % val(:,:,:,:,1) &
                     , u_l    => u    % level(l  ) % val(:,:,:,:,1) &
                     , r_l    => r    % level(l  ) % val(:,:,:,:,1) &
                     , nu_p   => nu_v % level(l-1) % val(:,:,:,:,1) &
                     , f_p    => f    % level(l-1) % val(:,:,:,:,1) &
                     , u_p    => u    % level(l-1) % val(:,:,:,:,1) &
                     , r_p    => r    % level(l-1) % val(:,:,:,:,1) &
                     , v_p    => v    % level(l-1) % val(:,:,:,:,1) )

              ! pre-smoothing and residual computation .........................

              n = this % ns_1
              if (present(nu_0)) then
                call this % Smoother(l, lambda, nu_0, u_l, f_l, bv_l, n)
                call this % Residual(l, lambda, nu_0, f_l, bv_l, u_l, r_l)
              else
                call this % Smoother(l, lambda, nu_l, u_l, f_l, bv_l, n)
                call this % Residual(l, lambda, nu_l, f_l, bv_l, u_l, r_l)
              end if
              call this % Monitoring(l, '1', r_l)

              ! restriction ....................................................

              ! project solution to regularly refined parent elements
              select case(this % projection_method)
              case('I')
                ! interpolation
                call ChildToParentProjection_3D &
                         (mesh_l, mesh_l, iop_fc(l), u_l, v_p)
              case('P')
                ! L²-projection
                call ChildToParentProjection_3D &
                         (mesh_l, mesh_l, pop_fc(l), u_l, v_p)
              end select

              ! restrict residual
              call ChildToParentRestriction_3D &
                       (mesh_l, mesh_p, iop_cf(l-1), r_l, r_p)

              ! parent RHS .....................................................

              do e = 1, mesh_p % n_elem
                if (mesh_p % element(e) % adaptation % refinement >= 1000) then
                  f_p(:,:,:,e) = -r_p(:,:,:,e)
                else
                  v_p(:,:,:,e) = u_p(:,:,:,e)
                end if
              end do

              if (present(nu_0)) then
                call ell_op(l-1) % Apply(lambda, nu_0, v_p, r_p)
              else
                call ell_op(l-1) % Apply(lambda, nu_p, v_p, r_p)
              end if

              do e = 1, mesh_p % n_elem
                if (mesh_p % element(e) % adaptation % refinement >= 1000) then
                  f_p(:,:,:,e) = f_p(:,:,:,e) + r_p(:,:,:,e)
                end if
              end do

            end associate
          end do FMG_DOWN

          FMG_COARSE: associate( bv_1 => bv   % level(1) % var            &
                               , nu_1 => nu_v % level(1) % val(:,:,:,:,1) &
                               , f_1  => f    % level(1) % val(:,:,:,:,1) &
                               , u_1  => u    % level(1) % val(:,:,:,:,1) )

            ! coarse grid solver ...............................................

            if (present(nu_0)) then
              call this % Monitoring(1, '0', lambda, nu_0, f_1, bv_1, u_1)
              call this % CoarseSolver(lambda, nu_0, u_1, f_1, bv_1)
              call this % Monitoring(1, 's', lambda, nu_0, f_1, bv_1, u_1)
            else
              call this % Monitoring(1, '0', lambda, nu_1, f_1, bv_1, u_1)
              call this % CoarseSolver(lambda, nu_1, u_1, f_1, bv_1)
              call this % Monitoring(1, 's', lambda, nu_1, f_1, bv_1, u_1)
            end if

          end associate FMG_COARSE

          FMG_UP: do l = 2, m

            associate( mesh_l => sem(l  ) % mesh                    &
                     , mesh_p => sem(l-1) % mesh                    &
                     , bv_l   => bv   % level(l  ) % var            &
                     , nu_l   => nu_v % level(l  ) % val(:,:,:,:,1) &
                     , f_l    => f    % level(l  ) % val(:,:,:,:,1) &
                     , u_l    => u    % level(l  ) % val(:,:,:,:,1) &
                     , v_l    => v    % level(l  ) % val(:,:,:,:,1) &
                     , v_p    => v    % level(l-1) % val(:,:,:,:,1) )

                ! prolongation .................................................

                call ParentToChildInterpolation_3D &
                         (mesh_p, mesh_l, iop_cf(l-1), v_p, v_l)

                ! update solution on current level
                if (mesh_p % n_elem > 0) then
                  n = mesh_p % n_elem_active
                  ! apply correction to active elements
                  call MergeArrays(ONE, u_l(:,:,:,:n), ONE, v_l(:,:,:,:n))
                  ! update frozen elements
                  if (mesh_p % n_elem_frozen > 0) then
                    call SetArray(u_l(:,:,:,n+1:), v_l(:,:,:,n+1:))
                  end if
                end if

              ! post-smoothing .................................................

              n = this % ns_2
              if (present(nu_0)) then
                call this % Monitoring(l, 'c', lambda, nu_0, f_l, bv_l, u_l)
                call this % Smoother(l, lambda, nu_0, u_l, f_l, bv_l, n)
                call this % Monitoring(l, '2', lambda, nu_0, f_l, bv_l, u_l)
              else
                call this % Monitoring(l, 'c', lambda, nu_l, f_l, bv_l, u_l)
                call this % Smoother(l, lambda, nu_l, u_l, f_l, bv_l, n)
                call this % Monitoring(l, '2', lambda, nu_l, f_l, bv_l, u_l)
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

            associate( mesh_l => sem(l  ) % mesh                    &
                     , mesh_p => sem(l-1) % mesh                    &
                     , bv_l   => bv   % level(l  ) % var            &
                     , nu_l   => nu_v % level(l  ) % val(:,:,:,:,1) &
                     , f_l    => f    % level(l  ) % val(:,:,:,:,1) &
                     , u_l    => u    % level(l  ) % val(:,:,:,:,1) &
                     , r_l    => r    % level(l  ) % val(:,:,:,:,1) &
                     , nu_p   => nu_v % level(l-1) % val(:,:,:,:,1) &
                     , f_p    => f    % level(l-1) % val(:,:,:,:,1) &
                     , u_p    => u    % level(l-1) % val(:,:,:,:,1) &
                     , r_p    => r    % level(l-1) % val(:,:,:,:,1) &
                     , v_p    => v    % level(l-1) % val(:,:,:,:,1) )

              ! pre-smoothing and residual computation .........................

              n = this % ns_1
              if (l < l_top .or. m == 1) then
                if (present(nu_0)) then
                  call this % Smoother(l, lambda, nu_0, u_l, f_l, bv_l, n)
                  call this % Residual(l, lambda, nu_0, f_l, bv_l, u_l, r_l)
                else
                  call this % Smoother(l, lambda, nu_l, u_l, f_l, bv_l, n)
                  call this % Residual(l, lambda, nu_l, f_l, bv_l, u_l, r_l)
                end if
              end if
              call this % Monitoring(l, '1', r_l)

              ! restriction ....................................................

              ! project solution to regularly refined parent elements
              select case(this % projection_method)
              case('I')
                ! interpolation
                call ChildToParentProjection_3D &
                         (mesh_l, mesh_l, iop_fc(l), u_l, v_p)
              case('P')
                ! L²-projection
                call ChildToParentProjection_3D &
                         (mesh_l, mesh_l, pop_fc(l), u_l, v_p)
              end select

              ! restrict residual
              call ChildToParentRestriction_3D &
                       (mesh_l, mesh_p, iop_cf(l-1), r_l, r_p)

              ! parent RHS .....................................................

              do e = 1, mesh_p % n_elem
                if (mesh_p % element(e) % adaptation % refinement >= 1000) then
                  f_p(:,:,:,e) = -r_p(:,:,:,e)
                else
                  v_p(:,:,:,e) = u_p(:,:,:,e)
                end if
              end do

              if (present(nu_0)) then
                call ell_op(l-1) % Apply(lambda, nu_0, v_p, r_p)
              else
                call ell_op(l-1) % Apply(lambda, nu_p, v_p, r_p)
              end if

              do e = 1, mesh_p % n_elem
                if (mesh_p % element(e) % adaptation % refinement >= 1000) then
                  f_p(:,:,:,e) = f_p(:,:,:,e) + r_p(:,:,:,e)
                end if
              end do

            end associate
          end do V_DOWN

          V_COARSE: associate( bv_1 => bv   % level(1) % var            &
                             , nu_1 => nu_v % level(1) % val(:,:,:,:,1) &
                             , f_1  => f    % level(1) % val(:,:,:,:,1) &
                             , u_1  => u    % level(1) % val(:,:,:,:,1) )

            ! coarse grid solver ...............................................

            if (present(nu_0)) then
              call this % Monitoring(1, '0', lambda, nu_0, f_1, bv_1, u_1)
              call this % CoarseSolver(lambda, nu_0, u_1, f_1, bv_1)
              call this % Monitoring(1, 's', lambda, nu_0, f_1, bv_1, u_1)
            else
              call this % Monitoring(1, '0', lambda, nu_1, f_1, bv_1, u_1)
              call this % CoarseSolver(lambda, nu_1, u_1, f_1, bv_1)
              call this % Monitoring(1, 's', lambda, nu_1, f_1, bv_1, u_1)
            end if

          end associate V_COARSE

          V_UP: do l = 2, m

            associate( mesh_l => sem(l  ) % mesh                    &
                     , mesh_p => sem(l-1) % mesh                    &
                     , bv_l   => bv   % level(l  ) % var            &
                     , nu_l   => nu_v % level(l  ) % val(:,:,:,:,1) &
                     , f_l    => f    % level(l  ) % val(:,:,:,:,1) &
                     , u_l    => u    % level(l  ) % val(:,:,:,:,1) &
                     , v_l    => v    % level(l  ) % val(:,:,:,:,1) &
                     , v_p    => v    % level(l-1) % val(:,:,:,:,1) )

                ! prolongation .................................................

                call ParentToChildInterpolation_3D &
                         (mesh_p, mesh_l, iop_cf(l-1), v_p, v_l)

                ! update solution on current level
                if (mesh_p % n_elem > 0) then
                  n = mesh_p % n_elem_active
                  ! apply correction to active elements
                  call MergeArrays(ONE, u_l(:,:,:,:n), ONE, v_l(:,:,:,:n))
                  ! update frozen elements
                  if (mesh_p % n_elem_frozen > 0) then
                    call SetArray(u_l(:,:,:,n+1:), v_l(:,:,:,n+1:))
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
                call this % Monitoring(l, 'c', lambda, nu_l, f_l, bv_l, u_l)
                call this % Smoother(l, lambda, nu_l, u_l, f_l, bv_l, n)
                call this % Monitoring(l, '2', lambda, nu_l, f_l, bv_l, u_l)
              end if

            end associate
          end do V_UP

        end do V_OUTER

      ! finalization :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      !$omp master
      deallocate(r, v)
      !$omp end master

    end associate

  end subroutine MG_Solver_X

  !=============================================================================

end submodule MP_MG_Solver
