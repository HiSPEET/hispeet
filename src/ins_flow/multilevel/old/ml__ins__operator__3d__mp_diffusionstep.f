!------------------------------------------------------------------------------!
! This file is part of HiSPEET: High-order Spectral Element Techniques         !
!                                                                              !
! Copyright (C) 2026 by the HiSPEET authors and the                            !
! Chair of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany                 !
!                                                                              !
! HiSPEET is free software: you can redistribute it and/or modify              !
! it under the terms of the GNU General Public License as published by         !
! the Free Software Foundation, either version 3 of the License, or            !
! (at your option) any later version.                                          !
!                                                                              !
! HiSPEET is distributed in the hope that it will be useful,                   !
! but WITHOUT ANY WARRANTY; without even the implied warranty of               !
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                         !
! See the GNU General Public License for more details.                         !
!                                                                              !
! You should have received a copy of the GNU General Public License            !
! along with HiSPEET. If not, see <http://www.gnu.org/licenses/>.              !
!------------------------------------------------------------------------------!

!> summary:  Multilevel Navier-Stokes diffusion step
!> author:   Joerg Stiller
!> date:     2026/07/02
!===============================================================================

submodule (ML__INS__Operator__3D) MP_DiffusionStep
  use Child_To_Parent_Projection__3D
  use Child_To_Parent_Restriction__3D
  use Parent_To_Child_Interpolation__3D
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Solves the implicit viscous subproblem

  module subroutine DiffusionStep(this, tau, mu, nu, bv, f, u, l_top)
    class(ML_INS_Operator_3D),     intent(in)    :: this
    real(RNP),                     intent(in)    :: tau   !< effective time step
    class(ML_MeshVariable_3D),     intent(in)    :: mu    !< bulk viscosity
    class(ML_MeshVariable_3D),     intent(in)    :: nu    !< shear viscosity
    class(ML_BoundaryVariable_3D), intent(in)    :: bv    !< boundary values
    class(ML_MeshVariable_3D),     intent(inout) :: f     !< RHS
    class(ML_MeshVariable_3D),     intent(inout) :: u     !< solution
    integer,             optional, intent(in)    :: l_top !< top level    [auto]

    ! internal variables :::::::::::::::::::::::::::::::::::::::::::::::::::::::

    type(ML_MeshVariable_3D), allocatable, save :: mm_inv ! inverse mass matrix
    type(ML_MeshVariable_3D), allocatable, save :: r      ! residual
    type(ML_MeshVariable_3D), allocatable, save :: w      ! workspace
    type(ML_MeshVariable_3D), allocatable, save :: z      ! zero

    character(len=:), allocatable :: log_prefix
    integer :: l_top_
    integer :: c, e, l, m, n

    associate( problem => this % problem                  &
             , sem     => this % ml_op_u % sem            &
             , iop_cf  => this % ml_op_u % iop_cf_x       &
             , iop_fc  => this % ml_op_u % iop_fc_x       &
             , pop_fc  => this % ml_op_u % pop_fc_x       &
             , n_cyc   => this % ml_diffusion_opt % i_max )

      ! initialization :::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      if (present(l_top)) then
        l_top_ = min(l_top, size(sem))
      else
        l_top_ = size(sem)
      end if

      !$omp master !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      associate(mesh => sem(1)%mesh)
        if (log_level == 1 .and. mesh%part == 0 .or. log_level > 1) then
          log_prefix = LoggingPrefix('DiffusionStep', mesh%part, mesh%n_parts)
        end if
      end associate

      if (allocated(log_prefix)) then
        print '(2A)', log_prefix, 'start'
      end if

      allocate(mm_inv, r, w, z)
      call mm_inv % Init(this%ml_op_u, nc = 1, l_top = l_top_)
      call r      % Init(this%ml_op_u, nc = 3, l_top = l_top_)
      call w      % Init(this%ml_op_u, nc = 3, l_top = l_top_)
      call z      % Init(this%ml_op_u, nc = 1, l_top = l_top_)

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

      V_OUTER: do m = 1, n_cyc

        V_DOWN: do l = l_top_, 2, -1
          associate( ins_l    => this % ins_op(l  )                     &
                   , ins_p    => this % ins_op(l-1)                     &
                   , mesh_l   => sem(l  ) % mesh                        &
                   , mesh_p   => sem(l-1) % mesh                        &
                   , bv_l     => bv     % level(l  ) % var              &
                   , mu_l     => mu     % level(l  ) % val(:,:,:,:, 1 ) &
                   , nu_l     => nu     % level(l  ) % val(:,:,:,:, 1 ) &
                   , f_l      => f      % level(l  ) % val(:,:,:,:,1:3) &
                   , r_l      => r      % level(l  ) % val(:,:,:,:,1:3) &
                   , v_l      => u      % level(l  ) % val(:,:,:,:,1:3) &
                   , bv_p     => bv     % level(l-1) % var              &
                   , mu_p     => mu     % level(l-1) % val(:,:,:,:, 1 ) &
                   , nu_p     => nu     % level(l-1) % val(:,:,:,:, 1 ) &
                   , f_p      => f      % level(l-1) % val(:,:,:,:,1:3) &
                   , r_p      => r      % level(l-1) % val(:,:,:,:,1:3) &
                   , v_p      => u      % level(l-1) % val(:,:,:,:,1:3) &
                   , w_p      => w      % level(l-1) % val(:,:,:,:,1:3) &
                   , z_p      => z      % level(l-1) % val(:,:,:,:, 1 ) &
                   , mm_inv_p => mm_inv % level(l-1) % val(:,:,:,:, 1 ) )

            ! pre-smoothing and residual computation ...........................

            if (l < l_top_ .or. m == 1) then
              n = this % ml_diffusion_opt % ns_1
              call ins_l % DiffusionSolver(tau, mu_l, nu_l, f_l, bv_l, v_l, n)
            end if

            call ins_l % GetDiffusionResidual( tau, mu_l, nu_l, f_l, bv_l, v_l &
                                             , r_l )

            ! restriction ......................................................

            ! project solution to regularly refined parent elements
            select case(this % ml_diffusion_opt % fc_projection)
            case('I')
              ! interpolation
              call ChildToParentProjection_3D &
                       (mesh_l, mesh_p, iop_fc(l), v_l, w_p)
            case('P')
              ! L²-projection
              call ChildToParentProjection_3D &
                       (mesh_l, mesh_p, pop_fc(l), v_l, w_p)
            end select

            ! restrict residual
            call ChildToParentRestriction_3D &
                     (mesh_l, mesh_p, iop_cf(l-1), r_l, r_p)

            ! unweighted parent RHS ............................................

            !$omp do
            do e = 1, mesh_p % n_elem
              if (mesh_p % element(e) % adaptation % refinement >= 1000) then
                f_p(:,:,:,e,1:3) = r_p(:,:,:,e,1:3)
                if (this % ml_diffusion_opt % fc_reset_twigs) then
                  ! reset twigs to projected child solution
                  v_p(:,:,:,e,1:3) = w_p(:,:,:,e,1:3)
                end if
              else
                w_p(:,:,:,e,1:3) = v_p(:,:,:,e,1:3)
              end if
            end do
            !$omp end do nowait

            ! diffusion operator -- with bulk diffusion
            call ins_p % ApplyDiffusionOperator(tau, mu_p, nu_p, bv_p, w_p, r_p)
!!!         ! diffusion operator -- without bulk diffusion, μ → z ≡ 0
!!!         call ins_p % ApplyDiffusionOperator(tau, mu_p, z_p, bv_p, v_p, r_p)

            !$omp do
            do e = 1, mesh_p % n_elem
              if (mesh_p % element(e) % adaptation % refinement >= 1000) then
                ! update RHS in regularly refined parent elements
                do c = 1, 3
                  f_p(:,:,:,e,c) = mm_inv_p(:,:,:,e) * ( f_p(:,:,:,e,c) &
                                                       + r_p(:,:,:,e,c) )
                end do
              end if
            end do
            !$omp end do nowait

          end associate
        end do V_DOWN

        ! coarse grid solver ...................................................

        call this % ins_op(1) % DiffusionSolver                            &
                                  ( tau = tau                              &
                                  , mu  = mu % level(1) % val(:,:,:,:, 1 ) &
                                  , nu  = nu % level(1) % val(:,:,:,:, 1 ) &
                                  , f   = f  % level(1) % val(:,:,:,:,1:3) &
                                  , bv  = bv % level(1) % var              &
                                  , v   = u  % level(1) % val(:,:,:,:,1:3) )

        V_UP: do l = 2, l_top_

          associate( ins_l  => this % ins_op(l)                   &
                   , mesh_l => sem(l  ) % mesh                    &
                   , mesh_p => sem(l-1) % mesh                    &
                   , bv_l   => bv % level(l  ) % var              &
                   , mu_l   => mu % level(l  ) % val(:,:,:,:, 1 ) &
                   , nu_l   => nu % level(l  ) % val(:,:,:,:, 1 ) &
                   , f_l    => f  % level(l  ) % val(:,:,:,:,1:3) &
                   , v_l    => u  % level(l  ) % val(:,:,:,:,1:3) &
                   , w_l    => r  % level(l  ) % val(:,:,:,:,1:3) &
                   , v_p    => u  % level(l-1) % val(:,:,:,:,1:3) &
                   , w_p    => w  % level(l-1) % val(:,:,:,:,1:3) )

            ! prolongation .....................................................

            !$omp do
            do e = 1, mesh_p % n_elem
              if (mesh_p % element(e) % adaptation % refinement >= 1000) then
                w_p(:,:,:,e,1:3) = v_p(:,:,:,e,1:3) - w_p(:,:,:,e,1:3)
              else
                w_p(:,:,:,e,1:3) = v_p(:,:,:,e,1:3)
              end if
            end do
            !$omp end do nowait

            call ParentToChildInterpolation_3D &
                     (mesh_p, mesh_l, iop_cf(l-1), w_p, w_l)

            !$omp do collapse(2)
            do c = 1, 3
            do e = 1, mesh_l%n_elem_active
              ! correct active elements
              v_l(:,:,:,e,c) = v_l(:,:,:,e,c) + w_l(:,:,:,e,c)
            end do
            end do
            !$omp end do nowait

            !$omp do collapse(2)
            do c = 1, 3
            do e = mesh_l%n_elem_active + 1, mesh_l%n_elem
              ! update frozen elements
              v_l(:,:,:,e,c) = w_l(:,:,:,e,c)
            end do
            end do
            !$omp end do nowait

            ! post-smoothing ...................................................

            if (l < l_top_ .or. m == n_cyc) then
              n = this % ml_diffusion_opt % ns_2
            else
              n = this % ml_diffusion_opt % ns_c
            end if

            call ins_l % DiffusionSolver(tau, mu_l, nu_l, f_l, bv_l, v_l, n)

          end associate
        end do V_UP

      end do V_OUTER

      ! finalization :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      !$omp master !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      deallocate(mm_inv, r, w, z)
      !$omp end master !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    end associate

    if (allocated(log_prefix)) then
      print '(2A)', log_prefix, 'exit'
    end if

  end subroutine DiffusionStep

  !=============================================================================

end submodule MP_DiffusionStep
