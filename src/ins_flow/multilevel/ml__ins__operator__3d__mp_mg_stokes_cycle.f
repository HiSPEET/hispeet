!> summary:  V-cycle for Stokes part in semi-implicit INS solvers
!> author:   Joerg Stiller
!> date:     2025/05/12
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> @remark
!> Auxiliary version skipping bulk viscosity part when restricting
!===============================================================================

submodule (ML__INS__Operator__3D) MP_MG_Stokes_Cycle
  use Child_To_Parent_Projection__3D
  use Child_To_Parent_Restriction__3D
  use Parent_To_Child_Interpolation__3D
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Performs one or more FAS-MG V-cycles for the Stokes part

  module subroutine MG_Stokes_Cycle(this, tau, mu, nu, bv, f, u, n_cyc, l_top)
    class(ML_INS_Operator_3D), intent(in) :: this
    real(RNP), intent(in) :: tau
      !< effective time step
    class(ML_MeshVariable_3D), intent(in) :: mu
      !< bulk viscosity
    class(ML_MeshVariable_3D), intent(in) :: nu
      !< shear viscosity
    class(ML_BoundaryVariable_3D), intent(in) :: bv
      !< boundary values
    class(ML_MeshVariable_3D), intent(inout) :: f
      !< unweighted RHS: f = v₀/τ + f_c + f_s + ...
    class(ML_MeshVariable_3D), intent(inout) :: u
      !< solution
    integer,  optional, intent(in) :: n_cyc
      !< number of cycles    [1]
    integer,  optional, intent(in) :: l_top
      !< top level        [auto]

    ! internal variables :::::::::::::::::::::::::::::::::::::::::::::::::::::::

    type(ML_MeshVariable_3D), allocatable, save :: mm_inv, r, w, z

    character(len=:), allocatable :: log_prefix
    integer :: l_top_, n_cyc_
    integer :: c, e, l, m

    associate( problem => this % problem            &
             , sem     => this % ml_op_u % sem      &
             , iop_cf  => this % ml_op_u % iop_cf_x &
             , iop_fc  => this % ml_op_u % iop_fc_x &
             , pop_fc  => this % ml_op_u % pop_fc_x )

      ! initialization :::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      if (present(n_cyc)) then
        n_cyc_ = n_cyc
      else
        n_cyc_ = 1
      end if

      if (present(l_top)) then
        l_top_ = min(l_top, size(sem))
      else
        l_top_ = size(sem)
      end if

      !$omp master !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      associate(mesh => sem(1)%mesh)
        if (log_level == 1 .and. mesh%part == 0 .or. log_level > 1) then
          log_prefix = LoggingPrefix('MG_Stokes_Cycle', mesh%part, mesh%n_parts)
        end if
      end associate

      if (allocated(log_prefix)) then
        print '(2A)', log_prefix, 'start'
      end if

      allocate(mm_inv, r, w, z)
      call mm_inv % Init(this%ml_op_u, nc = 1         , l_top = l_top_)
      call r      % Init(this%ml_op_u, nc = problem%nc, l_top = l_top_)
      call w      % Init(this%ml_op_u, nc = problem%nc, l_top = l_top_)
      call z      % Init(this%ml_op_u, nc = 1         , l_top = l_top_)

      !$omp end master !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! inverse mass matrix
      do l = 1, l_top_
        associate(mm_inv_l => mm_inv % level(l) % val(:,:,:,:,1))
          call sem(l) % Get_DG_DiagonalMassMatrix(mm_inv_l)
          !$omp do
          do e = 1, sem(l) % mesh % n_elem
            mm_inv_l(:,:,:,e) = 1 / mm_inv_l(:,:,:,e)
          end do
          !$omp end do nowait
        end associate
      end do

      ! V cycles :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      V_OUTER: do m = 1, n_cyc_

        V_DOWN: do l = l_top_, 2, -1

          associate( ins_l    => this % ins_op(l  )                   &
                   , ins_p    => this % ins_op(l-1)                   &
                   , mesh_l   => sem(l  ) % mesh                      &
                   , mesh_p   => sem(l-1) % mesh                      &
                   , bv_l     => bv     % level(l  ) % var            &
                   , mu_l     => mu     % level(l  ) % val(:,:,:,:,1) &
                   , nu_l     => nu     % level(l  ) % val(:,:,:,:,1) &
                   , f_l      => f      % level(l  ) % val            &
                   , u_l      => u      % level(l  ) % val            &
                   , r_l      => r      % level(l  ) % val            &
                   , z_l      => z      % level(l  ) % val(:,:,:,:,1) &
                   , bv_p     => bv     % level(l-1) % var            &
                   , mu_p     => mu     % level(l-1) % val(:,:,:,:,1) &
                   , nu_p     => nu     % level(l-1) % val(:,:,:,:,1) &
                   , f_p      => f      % level(l-1) % val            &
                   , u_p      => u      % level(l-1) % val            &
                   , r_p      => r      % level(l-1) % val            &
                   , w_p      => w      % level(l-1) % val            &
                   , z_p      => z      % level(l-1) % val(:,:,:,:,1) &
                   , mm_inv_p => mm_inv % level(l-1) % val(:,:,:,:,1) )

            if (log_level_multigrid_cycle > 0) then
              !$omp master
              if (mesh_l % part == 0) then
                write(*,'(99(G0,X))') 'MG Stokes cycle',m, 'down l =',l
              end if
              !$omp master
            end if

            ! pre-smoothing and residual computation ...........................

            if (l < l_top_) then
              call ins_l % StokesSolver(tau, f_l, bv_l, mu_l, nu_l, u_l)
            end if

            ! residual -- with bulk diffusion
            call ins_l % GetStokesResidual(tau, f_l, bv_l, mu_l, nu_l, u_l, r_l)
!!!         ! residual -- without bulk diffusion, μ → z ≡ 0
!!!         call ins_l % GetStokesResidual(tau, f_l, bv_l, z_l, nu_l, u_l, r_l)

            ! restriction ......................................................

            ! project solution to regularly refined parent elements
            select case(this % ml_diffusion_opt % fc_projection)
            case('I')
              ! interpolation
              call ChildToParentProjection_3D &
                       (mesh_l, mesh_p, iop_fc(l), u_l, w_p)
            case('P')
              ! L²-projection
              call ChildToParentProjection_3D &
                       (mesh_l, mesh_p, pop_fc(l), u_l, w_p)
            end select

            ! restrict residual
            call ChildToParentRestriction_3D &
                     (mesh_l, mesh_p, iop_cf(l-1), r_l, r_p)

            ! unweighted parent RHS ............................................

            !$omp do
            do e = 1, mesh_p % n_elem
              if (mesh_p % element(e) % adaptation % refinement >= 1000) then
                f_p(:,:,:,e,1:4) = r_p(:,:,:,e,1:4)
              else
                w_p(:,:,:,e,1:4) = u_p(:,:,:,e,1:4)
              end if
            end do
            !$omp end do nowait

            ! Stokes operator -- with bulk diffusion
            call ins_p % ApplyStokesOperator(tau, bv_p, mu_p, nu_p, w_p, r_p)
!!!         ! Stokes operator -- without bulk diffusion, μ → z ≡ 0
!!!         call ins_p % ApplyStokesOperator(tau, bv_p, z_p, nu_p, w_p, r_p)

            !$omp do
            do e = 1, mesh_p % n_elem
              if (mesh_p % element(e) % adaptation % refinement >= 1000) then
                do c = 1, 4
                  f_p(:,:,:,e,c) = mm_inv_p(:,:,:,e) &
                                 * (f_p(:,:,:,e,c) + r_p(:,:,:,e,c))
                end do
              end if
            end do
            !$omp end do nowait

          end associate

        end do V_DOWN

        ! coarse grid solver ...................................................

        if (log_level_multigrid_cycle > 0) then
          !$omp master
          if (this % ins_op(1) % mesh % part == 0) then
            write(*,'(99(G0,X))') 'MG Stokes cycle',m, 'coarse'
          end if
          !$omp master
        end if

        call this % ins_op(1) % StokesSolver( tau                            &
                                            , f  % level(1) % val            &
                                            , bv % level(1) % var            &
                                            , mu % level(1) % val(:,:,:,:,1) &
                                            , nu % level(1) % val(:,:,:,:,1) &
                                            , u  % level(1) % val            )

        V_UP: do l = 2, l_top_

          associate( ins_l  => this % ins_op(l)                 &
                   , mesh_l => sem(l  ) % mesh                  &
                   , mesh_p => sem(l-1) % mesh                  &
                   , bv_l   => bv % level(l  ) % var            &
                   , mu_l   => mu % level(l  ) % val(:,:,:,:,1) &
                   , nu_l   => nu % level(l  ) % val(:,:,:,:,1) &
                   , f_l    => f  % level(l  ) % val            &
                   , u_l    => u  % level(l  ) % val            &
                   , w_l    => r  % level(l  ) % val            &
                   , u_p    => u  % level(l-1) % val            &
                   , w_p    => w  % level(l-1) % val            )

            if (log_level_multigrid_cycle > 0) then
              !$omp master
              if (mesh_l % part == 0) then
                write(*,'(99(G0,X))') 'MG Stokes cycle',m, 'up l =',l
              end if
              !$omp master
            end if

            ! prolongation .....................................................

            !$omp do
            do e = 1, mesh_p % n_elem
              if (mesh_p % element(e) % adaptation % refinement >= 1000) then
                w_p(:,:,:,e,1:4) = u_p(:,:,:,e,1:4) - w_p(:,:,:,e,1:4)
              else
                w_p(:,:,:,e,1:4) = u_p(:,:,:,e,1:4)
              end if
            end do
            !$omp end do nowait

            call ParentToChildInterpolation_3D &
                     (mesh_p, mesh_l, iop_cf(l-1), w_p, w_l)

            !$omp do collapse(2)
            do c = 1, 4
            do e = 1, mesh_l%n_elem_active
              ! correct active elements
              u_l(:,:,:,e,c) = u_l(:,:,:,e,c) + w_l(:,:,:,e,c)
            end do
            end do
            !$omp end do nowait

            !$omp do collapse(2)
            do c = 1, 4
            do e = mesh_l%n_elem_active + 1, mesh_l%n_elem
              ! update frozen elements
              u_l(:,:,:,e,c) = w_l(:,:,:,e,c)
            end do
            end do
            !$omp end do nowait

            ! post-smoothing ...................................................

            call ins_l % StokesSolver(tau, f_l, bv_l, mu_l, nu_l, u_l)

          end associate
        end do V_UP

      end do V_OUTER

      !$omp master !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      deallocate(mm_inv, r, w, z)
      !$omp end master !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    end associate

    if (allocated(log_prefix)) then
      print '(2A)', log_prefix, 'exit'
    end if

  end subroutine MG_Stokes_Cycle

  !=============================================================================

end submodule MP_MG_Stokes_Cycle
