!> summary:  Cascade and FMG start procedures for multilevel elliptic solver
!> author:   Joerg Stiller
!> date:     2025/09/04
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(ML__DG__Elliptic_Solver__3D) MP_FAS_MG_Start
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Cascade and FMG start for problems with constant or variable diffusivity
  !>
  !> Either `nu_0` or `nu_v` must be given.

  module subroutine FAS_MG_Start_X( this, bc, lambda, nu_0, nu_v, bv, f, u &
                                  , n_cyc )

    class(ML_DG_EllipticSolver_3D), intent(in) :: this

    character, intent(in) :: bc(:)
      !< boundary conditions
    real(RNP), intent(in) :: lambda
      !< Helmholtz parameter
    real(RNP), optional, intent(in) :: nu_0
      !< constant diffusivity
    class(ML_MeshVariable_3D), optional, intent(in) :: nu_v
      !< variable diffusivity
    class(ML_BoundaryVariable_3D), intent(in) :: bv
      !< boundary values
    class(ML_MeshVariable_3D), intent(in) :: f
      !< RHS
    class(ML_MeshVariable_3D), intent(inout) :: u
      !< approx/final solution
    integer, intent(in) :: n_cyc
      !< number of V-cycles before advancing to next level,
      !! `n_cyc = 0` yields cascade and
      !! `n_cyc < 0` return with no change

    integer :: l, n

    if (n_cyc < 0) return

    associate( sem    => this % ml_op % sem      &
             , iop_cf => this % ml_op % iop_cf_x )

      do l = 1, size(sem) - 1
        associate( bv_l => bv % level(l)   % var            &
                 , f_l  => f  % level(l)   % val(:,:,:,:,1) &
                 , u_l  => u  % level(l)   % val(:,:,:,:,1) &
                 , u_c  => u  % level(l+1) % val(:,:,:,:,1) )

          if (l == 1) then

            ! coarse solver ....................................................

            if (present(nu_0)) then
              call this % CoarseSolver(bc, lambda, nu_0, u_l, f_l, bv_l)
              call this % Monitoring(l, 's', bc, lambda, nu_0, f_l, bv_l, u_l)
            else
              associate(nu_l => nu_v % level(l) % val(:,:,:,:,1))
                call this % CoarseSolver(bc, lambda, nu_l, u_l, f_l, bv_l)
                call this % Monitoring(l, 's', bc, lambda, nu_l, f_l, bv_l, u_l)
              end associate
            end if

          else if (n_cyc == 0) then

            ! cascade ..........................................................

            n = this % NumSmoothingSteps(l, stage = 0)
            if (present(nu_0)) then
              call this % Smoother(l, bc, lambda, nu_0, u_l, f_l, bv_l, n)
              call this % Monitoring(l, 's', bc, lambda, nu_0, f_l, bv_l, u_l)
            else
              associate(nu_l => nu_v % level(l) % val(:,:,:,:,1))
                call this % Smoother(l, bc, lambda, nu_l, u_l, f_l, bv_l, n)
                call this % Monitoring(l, 's', bc, lambda, nu_l, f_l, bv_l, u_l)
              end associate
            end if

          else

            ! full multigrid ...................................................

            call this % FAS_MG_Cycle_X( bc, lambda, nu_0, nu_v, bv, f, u &
                                      , n_cyc, l )

          end if

          ! interpolation to next level ........................................

          call ParentToChildInterpolation_3D( parent = sem(l  ) % mesh &
                                            , child  = sem(l+1) % mesh &
                                            , iop    = iop_cf(l)       &
                                            , v_p    = u_l             &
                                            , v_c    = u_c             )

        end associate
      end do

    end associate

  end subroutine FAS_MG_Start_X

  !=============================================================================

end submodule MP_FAS_MG_Start
