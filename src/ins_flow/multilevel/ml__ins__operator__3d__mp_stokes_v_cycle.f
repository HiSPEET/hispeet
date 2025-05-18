!> summary:  Stokes V-cycle for Stokes part in semi-implicit INS solvers
!> author:   Joerg Stiller
!> date:     2025/05/12
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule (ML__INS__Operator__3D) MP_Stokes_V_Cycle
  use Child_To_Parent_Projection__3D
  use Child_To_Parent_Restriction__3D
  use Parent_To_Child_Interpolation__3D
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Performs one or more FAS-MG V-cycles for the Stokes part

  module subroutine Stokes_V_Cycle(this, tau, mu, nu, bv, f, u, n_cyc, l_top)
    class(ML_INS_Operator_3D), intent(in) :: this
    real(RNP), intent(in) :: tau
      !< effective time step
    class(ML_MeshVariable_3D), optional, intent(in) :: mu
      !< bulk viscosity
    class(ML_MeshVariable_3D), optional, intent(in) :: nu
      !< shear viscosity
    class(ML_BoundaryVariable_3D), intent(in) :: bv
      !< boundary values
    class(ML_MeshVariable_3D), intent(inout) :: f
      !< RHS
    class(ML_MeshVariable_3D), intent(inout) :: u
      !< solution
    integer,  optional, intent(in) :: n_cyc
      !< number of cycles    [1]
    integer,  optional, intent(in) :: l_top
      !< top level        [auto]

    ! internal variables :::::::::::::::::::::::::::::::::::::::::::::::::::::::

    type(ML_MeshVariable_3D), allocatable, save :: mi, r, w
    real(RNP), contiguous, pointer, save :: mu_l(:,:,:,:), nu_l(:,:,:,:)
    real(RNP), contiguous, pointer, save :: mu_p(:,:,:,:), nu_p(:,:,:,:)
    integer, save :: l_top_, n_cyc_

    integer :: c, e, l, m, n

    associate( problem => this % problem            &
             , sem     => this % ml_op_u % sem      &
             , iop_cf  => this % ml_op_u % iop_cf_x &
             , iop_fc  => this % ml_op_u % iop_fc_x &
             , pop_fc  => this % ml_op_u % pop_fc_x )

      ! initialization :::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      !$omp master !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

      mu_l => null()
      nu_l => null()
      mu_p => null()
      nu_p => null()

      allocate(mi, r, w)
      call mi % Init(this%ml_op_u, nc = problem%nc, l_top = l_top_)
      call r  % Init(this%ml_op_u, nc = problem%nc, l_top = l_top_)
      call w  % Init(this%ml_op_u, nc = problem%nc, l_top = l_top_)

      !$omp end master !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! inverse mass matrix
      do l = 1, l_top_
        associate(mi_l => mi % level(l) % val(:,:,:,:,1))
          call sem(l) % Get_DG_DiagonalMassMatrix(mi_l)
          !$omp do
          do e = 1, sem(l) % mesh % n_elem
            mi_l(:,:,:,e) = 1 / mi_l(:,:,:,e)
          end do
          !$omp end do nowait
        end associate
      end do

      ! V cycles :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      V_OUTER: do m = 1, n_cyc_

        V_DOWN: do l = l_top_, 2, -1

          associate( ins_l  => this % ins_op(l)                  &
                   , mesh_l => sem(l  ) % mesh                   &
                   , mesh_p => sem(l-1) % mesh                   &
                   , bv_l   => bv % level(l  ) % var            &
                   , f_l    => f  % level(l  ) % val            &
                   , u_l    => u  % level(l  ) % val            &
                   , r_l    => r  % level(l  ) % val            &
                   , bv_p   => bv % level(l-1) % var            &
                   , f_p    => f  % level(l-1) % val            &
                   , u_p    => u  % level(l-1) % val            &
                   , r_p    => r  % level(l-1) % val            &
                   , w_p    => w  % level(l-1) % val            &
                   , mi_p   => mi % level(l-1) % val(:,:,:,:,1) )

            ! pre-smoothing and residual computation ...........................

            !$omp master
            if (present(mu) .and. present(nu)) then
              mu_l => mu % level(l  ) % val(:,:,:,:,1)
              nu_l => nu % level(l  ) % val(:,:,:,:,1)
              mu_p => mu % level(l-1) % val(:,:,:,:,1)
              nu_p => nu % level(l-1) % val(:,:,:,:,1)
            end if
            !$omp end master
            !!omp barrier after initialization in StokesProjection/FGMRES

            call ins_l % StokesSolver(tau, f_l, bv_l, mu_l, nu_l, u_l)
            call ins_l % GetStokesResidual(tau, f_l, bv_l, mu_l, nu_l, u_l, r_l)

            ! restriction ......................................................

            ! project solution to regularly refined parent elements
            select case(this % fc_project)
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

            ! parent RHS .......................................................

            !$omp do
            do e = 1, mesh_p % n_elem
              if (mesh_p % element(e) % adaptation % refinement >= 1000) then
                f_p(:,:,:,e,1:4) = r_p(:,:,:,e,1:4)
              else
                w_p(:,:,:,e,1:4) = u_p(:,:,:,e,1:4)
              end if
            end do
            !$omp end do nowait

            call ins_l % ApplyStokesOperator(tau, bv_p, mu_p, nu_p, u_p, r_p)

            !$omp do
            do e = 1, mesh_p % n_elem
              if (mesh_p % element(e) % adaptation % refinement < 1000) cycle
              do c = 1, 4
                f_p(:,:,:,e,c) = mi_p(:,:,:,e)*(f_p(:,:,:,e,c) + r_p(:,:,:,e,c))
              end do
            end do
            !$omp end do nowait

          end associate

        end do V_DOWN

        ! coarse grid solver ...................................................

        l = 1

        V_COARSE: associate( ins_l  => this % ins_op(l)      &
                           , bv_l   => bv % level(l  ) % var &
                           , f_l    => f  % level(l  ) % val &
                           , u_l    => u  % level(l  ) % val )

          !$omp master
          if (present(mu) .and. present(nu)) then
            mu_l => mu % level(l) % val(:,:,:,:,1)
            nu_l => nu % level(l) % val(:,:,:,:,1)
          end if
          !$omp end master

          call ins_l % StokesSolver(tau, f_l, bv_l, mu_l, nu_l, u_l)

        end associate V_COARSE

        V_UP: do l = 2, l_top_

          associate( ins_l  => this % ins_op(l)      &
                   , mesh_l => sem(l  ) % mesh       &
                   , mesh_p => sem(l-1) % mesh       &
                   , bv_l   => bv % level(l  ) % var &
                   , f_l    => f  % level(l  ) % val &
                   , u_l    => u  % level(l  ) % val &
                   , w_l    => r  % level(l  ) % val &
                   , u_p    => u  % level(l-1) % val &
                   , w_p    => w  % level(l-1) % val )

          ! prolongation .......................................................

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

          if (mesh_l % n_elem > 0) then
            n = mesh_l % n_elem_active

            !$omp do collapse(2)
            do c = 1, 4
            do e = 1, n
              ! apply correction to active elements
              u_l(:,:,:,e,c) = u_l(:,:,:,e,c) + w_l(:,:,:,e,c)
            end do
            end do
            !$omp end do nowait

            if (mesh_l % n_elem_frozen > 0) then
              !$omp do collapse(2)
              do c = 1, 4
              do e = n+1, mesh_l%n_elem
                ! update frozen elements
                u_l(:,:,:,e,c) =  w_l(:,:,:,e,c)
              end do
              end do
              !$omp end do nowait
            end if

          end if

          ! post-smoothing .....................................................

          !$omp master
          if (present(mu) .and. present(nu)) then
            mu_l => mu % level(l) % val(:,:,:,:,1)
            nu_l => nu % level(l) % val(:,:,:,:,1)
          end if
          !$omp end master

          call ins_l % StokesSolver(tau, f_l, bv_l, mu_l, nu_l, u_l)

          end associate
        end do V_UP

      end do V_OUTER

      !$omp master
      deallocate(mi, r, w)
      !$omp end master

    end associate

  end subroutine Stokes_V_Cycle

  !=============================================================================

end submodule MP_Stokes_V_Cycle
