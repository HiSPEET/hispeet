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

!> summary:  Cascade and FMG start procedures for Stokes multigrid solver
!> author:   Joerg Stiller
!> date:     2025/05/19
!===============================================================================

submodule (ML__INS__Operator__3D) MP_Stokes_MG_Start
  use Parent_To_Child_Interpolation__3D
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Cascade and FMG start for the Stokes multigrid solver

  module subroutine Stokes_MG_Start(this, tau, mu, nu, bv, f_d0, f, u, n_cyc)
    class(ML_INS_Operator_3D), intent(in) :: this
    real(RNP), intent(in) :: tau
      !< effective time step
    class(ML_MeshVariable_3D), intent(in) :: mu
      !< bulk viscosity
    class(ML_MeshVariable_3D), intent(in) :: nu
      !< shear viscosity
    class(ML_BoundaryVariable_3D), intent(in) :: bv
      !< boundary values
    class(ML_MeshVariable_3D), intent(in) :: f_d0
      !< unweighted approximate diffusion term
    class(ML_MeshVariable_3D), intent(inout) :: f
      !< unweighted RHS: f = v₀/τ + f_c + f_s + ...
    class(ML_MeshVariable_3D), intent(inout) :: u
      !< solution
    integer,  optional, intent(in) :: n_cyc
      !< number of V-cycles before advancing to next level  [0]

    character(len=:), allocatable :: log_prefix
    integer :: l, l_top

    !$omp master
    associate(mesh => this % ins_op(1) % mesh)
      if (log_level == 1 .and. mesh%part == 0 .or. log_level > 1) then
        log_prefix = LoggingPrefix('Stokes_MG_Start', mesh%part, mesh%n_parts)
      end if
    end associate
    !$omp end master

    if (allocated(log_prefix)) then
      print '(2A)', log_prefix, 'start'
    end if

    l_top = size(u%level)

    do l = 1, l_top

      call this % ins_op(l) % StokesSolver( tau                            &
                                          , f    % level(l)%val            &
                                          , bv   % level(l)%var            &
                                          , mu   % level(l)%val(:,:,:,:,1) &
                                          , nu   % level(l)%val(:,:,:,:,1) &
                                          , u    % level(l)%val            &
                                          , f_d0 % level(l)%val            )

      if (l == l_top) exit

      if (present(n_cyc)) then
        if (n_cyc > 0) then
          call this % Stokes_MG_Cycle(tau, mu, nu, bv, f, u, n_cyc, l_top = l)
        end if
      end if

      call ParentToChildInterpolation_3D                   &
               ( parent = this % ml_op_u % sem(l  ) % mesh &
               , child  = this % ml_op_u % sem(l+1) % mesh &
               , iop    = this % ml_op_u % iop_cf_x(l)     &
               , v_p    = u % level(l  ) % val             &
               , v_c    = u % level(l+1) % val             )

    end do

    if (allocated(log_prefix)) then
      print '(2A)', log_prefix, 'exit'
    end if

  end subroutine Stokes_MG_Start

  !=============================================================================

end submodule MP_Stokes_MG_Start
