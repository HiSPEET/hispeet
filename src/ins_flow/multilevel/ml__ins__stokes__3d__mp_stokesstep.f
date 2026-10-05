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

!> summary:  Multilevel Stokes step
!> author:   Joerg Stiller
!> date:     2026/10/03
!===============================================================================

submodule (ML__INS__Stokes__3D) MP_StokesStep
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Coupled solution of the projection and viscous diffusion subproblems

  module subroutine StokesStep(this, tau, mu, nu, bv, f_d0, f, u)
    class(ML_INS_Stokes_3D),       intent(in)    :: this
    real(RNP),                     intent(in)    :: tau  !< step size
    class(ML_MeshVariable_3D),     intent(in)    :: mu   !< bulk viscosity
    class(ML_MeshVariable_3D),     intent(in)    :: nu   !< shear viscosity
    class(ML_BoundaryVariable_3D), intent(in)    :: bv   !< boundary values
    class(ML_MeshVariable_3D),     intent(in)    :: f_d0 !< approx diffusion RHS
    class(ML_MeshVariable_3D),     intent(inout) :: f    !< RHS
    class(ML_MeshVariable_3D),     intent(inout) :: u    !< solution

    type(ML_MeshVariable_3D), allocatable, save :: r
    real(RNP) :: r0_2

    ! initial residual norm ....................................................

    if (max(this%r_red, this%r_max) > 0) then

      !$omp master
      allocate(r)
      call r % Init(this%ml_ins%ml_op_u, nc = this%ml_ins%problem%nc)
      !$omp end master
      !$omp barrier !? needed ???

      call this % StokesResidual(tau,mu, nu, bv, f, u, r)
      r0_2 = sqrt(ML_ScalarProduct_3D(r, r))

      !$omp master
      deallocate(r)
      !$omp end master

    else
      r0_2 = -1
    end if

    ! FAS-FMG solver ...........................................................

    call this % StokesStart(tau, mu, nu, bv, f_d0, f, u)
    call this % StokesCycle(tau, mu, nu, bv, f, u, r0_2 = r0_2)

  end subroutine StokesStep

  !=============================================================================

end submodule MP_StokesStep
