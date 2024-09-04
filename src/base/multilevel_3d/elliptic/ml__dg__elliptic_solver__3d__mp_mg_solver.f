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
    integer, save :: l_top

    integer :: l

    associate( sem         => this % ml_op % sem     &
             , iop_cf      => this % ml_op % iop_cf  &
             , iop_fc      => this % ml_op % iop_fc  &
             , pop_fc      => this % ml_op % pop_fc  &
             , elliptic_op => this % elliptic_op     )

      ! initialization :::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      !$omp master
      l_top = size(sem)
      r = ML_MeshVariable_3D(this%ml_op, nc=1, name=['r'])
      v = ML_MeshVariable_3D(this%ml_op, nc=1, name=['v'])
      !$omp end master

      do l = 1, l_top
        call SetArray(v % level(l) % val(:,:,:,:,1), ZERO)
      end do

      ! start ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      START: select case(this % start_method)

      case(START_CASC)

        ! cascade ..............................................................

        CASC_ROOT: associate( u_1  => u    % level(1) % val(:,:,:,:,1) &
                            , f_1  => f    % level(1) % val(:,:,:,:,1) &
                            , nu_1 => nu_v % level(1) % val(:,:,:,:,1) &
                            , bv_1 => bv   % level(1) % var            )

          if (present(nu_0)) then
            call this % CoarseSolver(lambda, nu_0, u_1, f_1, bv_1)
          else
            call this % CoarseSolver(lambda, nu_1, u_1, f_1, bv_1)
          end if

        end associate CASC_ROOT

        CASC_FINE: do l = 2, l_top
          associate( u_p  => u    % level(l-1) % val(:,:,:,:,1) &
                   , u_l  => u    % level(l  ) % val(:,:,:,:,1) &
                   , f_l  => f    % level(l  ) % val(:,:,:,:,1) &
                   , nu_l => nu_v % level(l  ) % val(:,:,:,:,1) &
                   , bv_l => bv   % level(l  ) % var            )

            ! interpolate
            call ParentToChildInterpolation_3D( parent = sem(l-1) % mesh &
                                              , child  = sem(l  ) % mesh &
                                              , iop    = iop_cf(l-1)     &
                                              , v_p    = u_p             &
                                              , v_c    = u_l             )
            ! smooth
            if (l < l_top) then
              if (present(nu_0)) then
                call this % Smoother(lambda, nu_0, u_l, f_l, bv_l, this%ns_2)
              else
                call this % Smoother(lambda, nu_l, u_l, f_l, bv_l, this%ns_2)
              end if
            end if

          end associate
        end do CASC_FINE

      case(START_FMG)


      case default

      end select START

      ! V cycles

      ! finalization

      !$omp master
      deallocate(r, v)
      !$omp end master

    end associate

  end subroutine MG_Solver_X

  !=============================================================================

end submodule MP_MG_Solver
