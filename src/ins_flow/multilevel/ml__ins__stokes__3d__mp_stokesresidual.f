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

!> summary:  Multilevel Stokes residual
!> author:   Joerg Stiller
!> date:     2026/10/04
!===============================================================================

submodule (ML__INS__Stokes__3D) MP_StokesResidual
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Coupled solution of the projection and viscous diffusion subproblems

  subroutine StokesResidual(this, tau, mu, nu, bv, f, u, r, l_top)

    class(ML_INS_Stokes_3D),       intent(in)    :: this
    real(RNP),                     intent(in)    :: tau   !< step size
    class(ML_MeshVariable_3D),     intent(in)    :: mu    !< bulk viscosity
    class(ML_MeshVariable_3D),     intent(in)    :: nu    !< shear viscosity
    class(ML_BoundaryVariable_3D), intent(in)    :: bv    !< boundary values
    class(ML_MeshVariable_3D),     intent(in)    :: f     !< RHS
    class(ML_MeshVariable_3D),     intent(in)    :: u     !< solution
    class(ML_MeshVariable_3D),     intent(inout) :: r     !< residual
    integer,             optional, intent(in)    :: l_top !< top level   [auto]

    ! internal variables :::::::::::::::::::::::::::::::::::::::::::::::::::::::

    type(ML_MeshVariable_3D), allocatable, save :: mm_inv ! inverse mass matrix
    type(ML_MeshVariable_3D), allocatable, save :: g, w

    integer :: c, e, l, l_top_

    associate( ml_ins  => this % ml_ins                       &
             , problem => this % ml_ins % problem             &
             , sem     => this % ml_ins % ml_op_u % sem       &
             , iop_cf  => this % ml_ins % ml_op_u % iop_cf_x  &
             , iop_fc  => this % ml_ins % ml_op_u % iop_fc_x  &
             , pop_fc  => this % ml_ins % ml_op_u % pop_fc_x  )

      ! initialization :::::::::::::::::::::::::::::::::::::::::::::::::::::::

      if (present(l_top)) then
        l_top_ = min(l_top, size(sem))
      else
        l_top_ = size(sem)
      end if

      !$omp master
      allocate(mm_inv, g, w)
      call mm_inv % Init(ml_ins%ml_op_u, nc = 1         , l_top = l_top_)
      call g      % Init(ml_ins%ml_op_u, nc = problem%nc, l_top = l_top_)
      call w      % Init(ml_ins%ml_op_u, nc = problem%nc, l_top = l_top_)
      !$omp end master

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
      !$omp barrier !? needed ???

      ! FAS residual :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      do l = l_top_, 1, -1
        associate( ins_l => ml_ins % ins_op(l)               &
                 , bv_l  => bv % level(l) % var              &
                 , mu_l  => mu % level(l) % val(:,:,:,:, 1 ) &
                 , nu_l  => nu % level(l) % val(:,:,:,:, 1 ) &
                 , f_l   => f  % level(l) % val              &
                 , g_l   => g  % level(l) % val              &
                 , u_l   => u  % level(l) % val              &
                 , r_l   => r  % level(l) % val              )

          ! residual on current level ..........................................

          if (l == l_top_) then
            ! using original RHS on top level
            call ins_l % GetStokesResidual(tau, mu_l, nu_l, bv_l, f_l, u_l, r_l)
          else
            ! using FAS-augmented RHS on lower levels
            call ins_l % GetStokesResidual(tau, mu_l, nu_l, bv_l, g_l, u_l, r_l)
          end if

          if (l == 1) exit

          ! parent FAS-RHS .....................................................

          associate( mesh_l   => sem(l  ) % mesh                      &
                   , mesh_p   => sem(l-1) % mesh                      &
                   , ins_p    => ml_ins % ins_op(l-1)                 &
                   , bv_p     => bv     % level(l-1) % var            &
                   , mu_p     => mu     % level(l-1) % val(:,:,:,:,1) &
                   , nu_p     => nu     % level(l-1) % val(:,:,:,:,1) &
                   , f_p      => f      % level(l-1) % val            &
                   , g_p      => g      % level(l-1) % val            &
                   , r_p      => r      % level(l-1) % val            &
                   , u_p      => u      % level(l-1) % val            &
                   , w_p      => w      % level(l-1) % val            &
                   , mm_inv_p => mm_inv % level(l-1) % val(:,:,:,:,1) )

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
                g_p(:,:,:,e,1:4) = r_p(:,:,:,e,1:4)
              else
                g_p(:,:,:,e,1:4) = f_p(:,:,:,e,1:4)
                w_p(:,:,:,e,1:4) = u_p(:,:,:,e,1:4)
              end if
            end do
            !$omp end do nowait

            call ins_p % ApplyStokesOperator(tau, bv_p, mu_p, nu_p, w_p, r_p)

            !$omp do
            do e = 1, mesh_p % n_elem
              if (mesh_p % element(e) % adaptation % refinement >= 1000) then
                ! update RHS in regularly refined parent elements
                do c = 1, 4
                  g_p(:,:,:,e,c) = mm_inv_p(:,:,:,e) * ( g_p(:,:,:,e,c) &
                                                       + r_p(:,:,:,e,c) )
                end do
              end if
            end do
            !$omp end do nowait

          end associate
        end associate
      end do

      ! finalization :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      !$omp master
      deallocate(mm_inv, g, w)
      !$omp end master

    end associate

  end subroutine StokesResidual

  !=============================================================================

end submodule MP_StokesResidual
