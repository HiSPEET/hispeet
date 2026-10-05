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

!> summary:  Multilevel Stokes cycle
!> author:   Joerg Stiller
!> date:     2026/10/03
!===============================================================================

submodule (ML__INS__Stokes__3D) MP_StokesCycle
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> FAS MG Stokes cycle

  module subroutine StokesCycle( this, tau, mu, nu, bv, f, u &
                               , n_cyc, l_top, r0_2 )

    class(ML_INS_Stokes_3D),       intent(in)    :: this
    real(RNP),                     intent(in)    :: tau  !< step size
    class(ML_MeshVariable_3D),     intent(in)    :: mu   !< bulk viscosity
    class(ML_MeshVariable_3D),     intent(in)    :: nu   !< shear viscosity
    class(ML_BoundaryVariable_3D), intent(in)    :: bv   !< boundary values
    class(ML_MeshVariable_3D),     intent(inout) :: f    !< RHS
    class(ML_MeshVariable_3D),     intent(inout) :: u    !< solution

    integer,   optional, intent(in)  :: n_cyc !< num cycles       [this%i_max]
    integer,   optional, intent(in)  :: l_top !< top level              [auto]
    real(RNP), optional, intent(in)  :: r0_2  !< initial residual norm, if < 0

    ! internal variables :::::::::::::::::::::::::::::::::::::::::::::::::::::::

    type(ML_MeshVariable_3D), allocatable, save :: mm_inv, r, w
    logical, save :: converged

    real(RNP) :: rr, r_max, r_new, r_old
    logical   :: check_convergence
    integer   :: l_top_, n_cyc_
    integer   :: c, e, l, m

    character(len=:), allocatable :: prefix
    logical :: logging

    associate( ml_ins  => this % ml_ins                      &
             , problem => this % ml_ins % problem            &
             , sem     => this % ml_ins % ml_op_u % sem      &
             , iop_cf  => this % ml_ins % ml_op_u % iop_cf_x &
             , iop_fc  => this % ml_ins % ml_op_u % iop_fc_x &
             , pop_fc  => this % ml_ins % ml_op_u % pop_fc_x )

      ! initialization :::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      if (present(n_cyc)) then
        n_cyc_ = n_cyc
        check_convergence = .false.
      else
        n_cyc_ = this % i_max
        check_convergence = max(this%r_red, this%r_max) > 0
      end if

      if (present(l_top)) then
        l_top_ = min(l_top, size(sem))
      else
        l_top_ = size(sem)
      end if
      if (l_top_ < 1) return

      if (log_level > 0) then
        associate(proc => sem(1) % mesh % proc)
          logging = proc == 0 .or. log_level_multigrid_cycle > 0
          prefix  = LoggingPrefix('StokesCycle', proc)
        end associate
      else
        logging = .false.
      end if

      !$omp master !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      allocate(mm_inv, r, w)
      call mm_inv % Init(ml_ins%ml_op_u, nc = 1         , l_top = l_top_)
      call r      % Init(ml_ins%ml_op_u, nc = problem%nc, l_top = l_top_)
      call w      % Init(ml_ins%ml_op_u, nc = problem%nc, l_top = l_top_)
      !$omp end master !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! termination condition
      if (check_convergence) then
        if (present(r0_2)) then
          r_old = r0_2
        else
          call this % StokesResidual(tau, mu, nu, bv, f, u, r, l_top_)
          rr = ML_ScalarProduct_3D(r, r, l_top = l_top_)
          r_old = sqrt(rr)
        end if
        r_max = max(r_old * this%r_red, this%r_max)
        !$omp master
        converged = r_old < r_max
        call XMPI_Bcast(converged, root = 0, comm = sem(1)%mesh%comm_world)
        if (logging) then
          write(*,'(2A,I0,A,ES12.5)') prefix, ': r_2(',0,') = ',r_old
        end if
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

          associate( ins_l    => ml_ins % ins_op(l  )                 &
                   , ins_p    => ml_ins % ins_op(l-1)                 &
                   , mesh_l   => sem(l  ) % mesh                      &
                   , mesh_p   => sem(l-1) % mesh                      &
                   , bv_l     => bv     % level(l  ) % var            &
                   , mu_l     => mu     % level(l  ) % val(:,:,:,:,1) &
                   , nu_l     => nu     % level(l  ) % val(:,:,:,:,1) &
                   , f_l      => f      % level(l  ) % val            &
                   , u_l      => u      % level(l  ) % val            &
                   , r_l      => r      % level(l  ) % val            &
                   , bv_p     => bv     % level(l-1) % var            &
                   , mu_p     => mu     % level(l-1) % val(:,:,:,:,1) &
                   , nu_p     => nu     % level(l-1) % val(:,:,:,:,1) &
                   , f_p      => f      % level(l-1) % val            &
                   , u_p      => u      % level(l-1) % val            &
                   , r_p      => r      % level(l-1) % val            &
                   , w_p      => w      % level(l-1) % val            &
                   , mm_inv_p => mm_inv % level(l-1) % val(:,:,:,:,1) )

            if (logging) then
              write(*,'(A,2(A,I0))') prefix, ': m = ',m,' down l =',l
            end if

            ! pre-smoothing and residual computation ...........................

            if (l < l_top_ .or. m == 1) then
              call ins_l % StokesSolver(tau, mu_l, nu_l, bv_l, f=f_l, u=u_l)
            end if

            ! residual -- with bulk diffusion
            call ins_l % GetStokesResidual(tau, mu_l, nu_l, bv_l, f_l, u_l, r_l)

            ! restriction ......................................................

            ! project solution to regularly refined parent elements
            select case(this % fc_projection)
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
                if (this % fc_reset_twigs) then
                  ! reset twigs to projected child solution
                  u_p(:,:,:,e,1:4) = w_p(:,:,:,e,1:4)
                end if
              else
                w_p(:,:,:,e,1:4) = u_p(:,:,:,e,1:4)
              end if
            end do
            !$omp end do nowait

            ! Stokes operator
            call ins_p % ApplyStokesOperator(tau, bv_p, mu_p, nu_p, w_p, r_p)

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

        if (logging) then
          write(*,'(A,2(A,I0))') prefix, ': m = ',m,' coarse'
        end if

        call ml_ins % ins_op(1) % StokesSolver( tau                            &
                                              , mu % level(1) % val(:,:,:,:,1) &
                                              , nu % level(1) % val(:,:,:,:,1) &
                                              , bv % level(1) % var            &
                                              , f = f % level(1) % val         &
                                              , u = u % level(1) % val         )

        V_UP: do l = 2, l_top_

          associate( ins_l  => ml_ins % ins_op(l)               &
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

            if (logging) then
              write(*,'(A,2(A,I0))') prefix, ': m = ',m,' up l =',l
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

            call ins_l % StokesSolver(tau, mu_l, nu_l, bv_l, f=f_l, u=u_l)

          end associate
        end do V_UP

        ! termination check ....................................................

        if (check_convergence .and. &
               (m < n_cyc_ .or. log_level_multigrid_cycle > 0)) then

          call this % StokesResidual(tau, mu, nu, bv, f, u, r, l_top_)
          rr = ML_ScalarProduct_3D(r, r, l_top = l_top_)
          r_new = sqrt(rr)
          if (logging) then
            write(*,'(2A,I0,A,ES12.5)') prefix, ': r_2(',m,') = ',r_new
          end if

          !$omp master
          converged = r_new <= r_max
          call XMPI_Bcast(converged, root = 0, comm = sem(1)%mesh%comm_world)
          !$omp end master
          !$omp barrier

          if (converged) exit
          r_old = r_new

        end if

      end do V_OUTER

      !$omp master !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      deallocate(mm_inv, r, w)
      !$omp end master !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    end associate

  end subroutine StokesCycle

  !=============================================================================

end submodule MP_StokesCycle
