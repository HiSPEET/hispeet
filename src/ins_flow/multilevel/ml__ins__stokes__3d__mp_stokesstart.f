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

!> summary:  Multilevel Stokes start method
!> author:   Joerg Stiller
!> date:     2026/10/04
!===============================================================================

submodule (ML__INS__Stokes__3D) MP_StokesStart
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Cascade and FMG start for the Stokes multigrid solver

  module subroutine StokesStart(this, tau, mu, nu, bv, f_d0, f, u)
    class(ML_INS_Stokes_3D),       intent(in)    :: this
    real(RNP),                     intent(in)    :: tau  !< step size
    class(ML_MeshVariable_3D),     intent(in)    :: mu   !< bulk viscosity
    class(ML_MeshVariable_3D),     intent(in)    :: nu   !< shear viscosity
    class(ML_BoundaryVariable_3D), intent(in)    :: bv   !< boundary values
    class(ML_MeshVariable_3D),     intent(in)    :: f_d0 !< approx diffusion RHS
    class(ML_MeshVariable_3D),     intent(inout) :: f    !< RHS
    class(ML_MeshVariable_3D),     intent(inout) :: u    !< solution

    integer :: l, l_top

    character(len=:), allocatable :: prefix
    logical :: logging

    if (this % i_fmg < 0) return

    l_top = size(u%level)

    logging = .false.
    if (log_level_multigrid_cycle > 0) then
      !$omp master
      associate(proc => this % ml_ins % ins_op(1) % mesh % proc)
        logging = proc == 0
        prefix  = LoggingPrefix('StokesStart', proc)
      end associate
      !$omp end master
    end if

    do l = 1, l_top

      if (logging) then
        write(*,'(2A,I0,A,ES12.5)') prefix, 'l = ',l
      end if

      call this % ml_ins % ins_op(l) %                                     &
                              StokesSolver( tau                            &
                                          , mu   % level(l)%val(:,:,:,:,1) &
                                          , nu   % level(l)%val(:,:,:,:,1) &
                                          , bv   % level(l)%var            &
                                          , f_d0 % level(l)%val            &
                                          , f    % level(l)%val            &
                                          , u    % level(l)%val            )

      if (l == l_top) exit

      if (this%i_fmg > 0) then
        call this % StokesCycle(tau, mu, nu, bv, f, u, this%i_fmg, l_top = l)
      end if

      call ParentToChildInterpolation_3D                            &
               ( parent = this % ml_ins % ml_op_u % sem(l  ) % mesh &
               , child  = this % ml_ins % ml_op_u % sem(l+1) % mesh &
               , iop    = this % ml_ins % ml_op_u % iop_cf_x(l)     &
               , v_p    = u % level(l  ) % val                      &
               , v_c    = u % level(l+1) % val                      )

    end do

  end subroutine StokesStart

  !=============================================================================

end submodule MP_StokesStart
