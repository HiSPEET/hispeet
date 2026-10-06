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
!> date:     2026/09/30
!===============================================================================

submodule (ML__INS__Diffusion__3D) MP_DiffusionStep
  implicit none

  interface Monitoring
    module procedure Monitoring_R
    module procedure Monitoring_V
  end interface

contains

  !-----------------------------------------------------------------------------
  !> Solution of the implicit viscous subproblem

  module subroutine DiffusionStep(this, tau, mu, nu, bv, f, u, l_top)
    class(ML_INS_Diffusion_3D),    intent(in)    :: this
    real(RNP),                     intent(in)    :: tau    !< step size
    class(ML_MeshVariable_3D),     intent(in)    :: mu     !< bulk viscosity
    class(ML_MeshVariable_3D),     intent(in)    :: nu     !< shear viscosity
    class(ML_BoundaryVariable_3D), intent(in)    :: bv     !< boundary values
    class(ML_MeshVariable_3D),     intent(inout) :: f      !< RHS
    class(ML_MeshVariable_3D),     intent(inout) :: u      !< solution
    integer,             optional, intent(in)    :: l_top  !< top level   [auto]

    ! internal variables :::::::::::::::::::::::::::::::::::::::::::::::::::::::

    type(ML_MeshVariable_3D), allocatable, save :: mm_inv ! inverse mass matrix
    type(ML_MeshVariable_3D), allocatable, save :: r      ! residual
    type(ML_MeshVariable_3D), allocatable, save :: w      ! workspace
    logical, save :: converged

    real(RNP) :: rr, r_max, r_new, r_old
    logical   :: check_convergence, logging
    integer   :: l_top_
    integer   :: c, e, l, m, n

    character(len=:), allocatable :: prefix

    associate( ml_ins  => this % ml_ins                       &
             , problem => this % ml_ins % problem             &
             , sem     => this % ml_ins % ml_op_u % sem       &
             , iop_cf  => this % ml_ins % ml_op_u % iop_cf_x  &
             , iop_fc  => this % ml_ins % ml_op_u % iop_fc_x  &
             , pop_fc  => this % ml_ins % ml_op_u % pop_fc_x  )

      ! initialization :::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      !$omp master
      associate(mesh => sem(1)%mesh)
        if (log_level == 1 .and. mesh%part == 0 .or. log_level > 1) then
          prefix = LoggingPrefix('DiffusionStep', mesh%proc)
        end if
      end associate
      !$omp end master

      logging = allocated(prefix)
      if (logging) then
        print '(2A)', prefix, 'start'
      end if

      if (present(l_top)) then
        l_top_ = min(l_top, size(sem))
      else
        l_top_ = size(sem)
      end if

      check_convergence = max(this%r_red, this%r_max) > 0

      !$omp master
      allocate(mm_inv, r, w)
      call mm_inv % Init(ml_ins % ml_op_u, nc = 1, l_top = l_top_)
      call r      % Init(ml_ins % ml_op_u, nc = 3, l_top = l_top_)
      call w      % Init(ml_ins % ml_op_u, nc = 3, l_top = l_top_)
      !$omp end master

      ! termination condition ..................................................

      if (check_convergence) then
        call this % DiffusionResidual(tau, mu, nu, bv, f, u, r, l_top_)
        rr = ML_ScalarProduct_3D(r, r, l_top = l_top_)
        r_old = sqrt(rr)
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
        deallocate(mm_inv, r, w)
        !$omp end master
        return
      end if

      ! inverse mass matrix ....................................................

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

      V_OUTER: do m = 1, this % i_max

        V_DOWN: do l = l_top_, 2, -1
          associate( ins_l    => ml_ins % ins_op(l  )                   &
                   , ins_p    => ml_ins % ins_op(l-1)                   &
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
                   , mm_inv_p => mm_inv % level(l-1) % val(:,:,:,:, 1 ) )

            ! pre-smoothing and residual computation ...........................

            call Monitoring(l, '0', ins_l, tau, mu_l, nu_l, f_l, bv_l, v_l)

            if (l < l_top_ .or. m == 1) then
              call ins_l % DiffusionSolver( tau, mu_l, nu_l, bv_l, f_l, v_l &
                                          , i_max = this%ns_1               &
                                          , r_red = -ONE                    &
                                          , r_max = -ONE                    )
            end if

            call ins_l % GetDiffusionResidual( tau, mu_l, nu_l, bv_l, f_l, v_l &
                                             , r_l )

            call Monitoring(l, '1', ins_l, r_l)

            ! restriction ......................................................

            ! project solution to regularly refined parent elements
            select case(this % fc_projection)
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
                if (this % fc_reset_twigs) then
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

        associate( ins_l    => ml_ins % ins_op(1)               &
                 , bv_l     => bv % level(1) % var              &
                 , mu_l     => mu % level(1) % val(:,:,:,:, 1 ) &
                 , nu_l     => nu % level(1) % val(:,:,:,:, 1 ) &
                 , f_l      => f  % level(1) % val(:,:,:,:,1:3) &
                 , r_l      => r  % level(1) % val(:,:,:,:,1:3) &
                 , v_l      => u  % level(1) % val(:,:,:,:,1:3) )

          call Monitoring(1, '0', ins_l, tau, mu_l, nu_l, f_l, bv_l, v_l)

          call ins_l % DiffusionSolver( tau, mu_l, nu_l, bv_l, f_l, v_l &
                                      , i_max = this % i_crs            &
                                      , r_red = this % r_crs            &
                                      , r_max = this % r_max            )

          call Monitoring(1, 's', ins_l, tau, mu_l, nu_l, f_l, bv_l, v_l)

        end associate

        V_UP: do l = 2, l_top_

          associate( ins_l  => ml_ins % ins_op(l)                 &
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

            if (l < l_top_ .or. m == this % i_max) then
              n = this % ns_2
            else
              n = this % ns_c
            end if

            call Monitoring(l, 'c', ins_l, tau, mu_l, nu_l, f_l, bv_l, v_l)

            call ins_l % DiffusionSolver( tau, mu_l, nu_l, bv_l, f_l, v_l &
                                        , i_max = n                       &
                                        , r_red = -ONE                    &
                                        , r_max = -ONE                    )

            call Monitoring(l, '2', ins_l, tau, mu_l, nu_l, f_l, bv_l, v_l)

          end associate
        end do V_UP

        ! termination check ....................................................

        if (check_convergence .and. m < this%i_max) then

          call this % DiffusionResidual(tau, mu, nu, bv , f, u, r, l_top_)

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

      ! finalization :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      !$omp master
      deallocate(mm_inv, r, w)
      !$omp end master

    end associate

    if (logging) then
      print '(2A)', prefix, 'exit'
    end if

  end subroutine DiffusionStep

  !-----------------------------------------------------------------------------
  !> Monitoring with known residual

  subroutine Monitoring_R(l, step, ins_op, r)
    integer,                intent(in) :: l            !< level
    character(len=*),       intent(in) :: step         !< current step
    class(INS_Operator_3D), intent(in) :: ins_op       !< INS operator
    real(RNP), contiguous,  intent(in) :: r(:,:,:,:,:) !< residual

    real(RNP) :: rr

    if (log_level_multigrid_cycle < 1 .or. ins_op%mesh%part < 0) return

    rr = ScalarProduct(r, r, ins_op%mesh%comm_parts)

    !$omp master
    if (ins_op%mesh%part == 0) then
      print '(2X,A,I2,3A,ES10.3)', 'level ', l, ': r[', step, '] =', sqrt(rr)
    end if
    !$omp end master

  end subroutine Monitoring_R

  !-----------------------------------------------------------------------------
  !> Monitoring with unknown residual

  subroutine Monitoring_V(l, step, ins_op, tau, mu, nu, f, bv, v)
    integer,                    intent(in) :: l            !< level
    character(len=*),           intent(in) :: step         !< current step
    class(INS_Operator_3D),     intent(in) :: ins_op       !< INS operator
    real(RNP),                  intent(in) :: tau          !< step size
    real(RNP), contiguous,      intent(in) :: mu(:,:,:,:)  !< residual
    real(RNP), contiguous,      intent(in) :: nu(:,:,:,:)  !< residual
    class(BoundaryVariable_3D), intent(in) :: bv(:)        !< boundary values
    real(RNP), contiguous,      intent(in) :: f(:,:,:,:,:) !< RHS
    real(RNP), contiguous,      intent(in) :: v(:,:,:,:,:) !< velocity

    real(RNP), allocatable, save :: r(:,:,:,:,:)

    if (log_level_multigrid_cycle < 1 .or. ins_op%mesh%part < 0) return

    !$omp master
    allocate(r, mold = v)
    !$omp end master
    !$omp barrier  !? needed

    call ins_op % GetDiffusionResidual(tau, mu, nu, bv, f, v, r)
    call Monitoring_R(l, step, ins_op, r)

    !$omp master
    deallocate(r)
    !$omp end master

  end subroutine Monitoring_V

  !=============================================================================

end submodule MP_DiffusionStep
