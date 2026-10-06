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

!> summary:  Base type of one-step ML time integrators for incompressible flows
!> author:   Joerg Stiller
!> date:     2025/05/12
!===============================================================================

module ML__INS__Integrator__3D
  use Kind_Parameters, only: RNP
  use XMPI
  use INS__Problem__3D
  use ML__Mesh_Variable__3D
  use ML__INS__Operator__3D
  implicit none
  private

  public :: ML_INS_Integrator_3D
  public :: ML_INS_IntegratorOptions_3D

  !-----------------------------------------------------------------------------
  !> Base type for one-step multilevel IMEX INS integrators

  type, abstract :: ML_INS_Integrator_3D
    class(INS_Problem_3D),    pointer :: problem => null() !< flow problem
    type(ML_INS_Operator_3D), pointer :: ml_ins  => null() !< ML INS operator
  contains
    procedure, non_overridable    :: Init_ML_INS_Integrator_3D
    procedure(TimeStep), deferred :: TimeStep
  end type ML_INS_Integrator_3D

  !=============================================================================
  ! deferred module procedures

  abstract interface

    !---------------------------------------------------------------------------
    !> Execution of a multilevel time step

    subroutine TimeStep(this, t, dt, mu, nu, u, first, last)
      import :: RNP, ML_MeshVariable_3D, ML_INS_Integrator_3D
      class(ML_INS_Integrator_3D), intent(inout) :: this
      real(RNP),                 intent(inout) :: t     !< time t₀ → t
      real(RNP),                 intent(in)    :: dt    !< step size ∆t = t-t₀
      class(ML_MeshVariable_3D), intent(inout) :: mu    !< bulk viscosity μ
      class(ML_MeshVariable_3D), intent(inout) :: nu    !< shear viscosity ν
      class(ML_MeshVariable_3D), intent(inout) :: u     !< u(x,t₀) → u(x,t)
      logical,         optional, intent(in)    :: first !< T for first step [F]
      logical,         optional, intent(in)    :: last  !< T for last  step [F]
    end subroutine TimeStep

  end interface

  !-----------------------------------------------------------------------------
  !> Base type for providing multilevel integrator options

  type ML_INS_IntegratorOptions_3D
  contains
    procedure :: Bcast => Bcast_ML_INS_IntegratorOptions_3D
  end type ML_INS_IntegratorOptions_3D

contains

  !=============================================================================
  ! ML_INS_Integrator_3D: type-bound procedures

  !-----------------------------------------------------------------------------
  !> Basic initialization of one-step multilevel IMEX INS integrators

  subroutine Init_ML_INS_Integrator_3D(this, ml_ins, opt)
    class(ML_INS_Integrator_3D),        intent(inout) :: this
    class(ML_INS_Operator_3D),  target, intent(in)    :: ml_ins
    class(ML_INS_IntegratorOptions_3D), intent(in)    :: opt

    this % problem => ml_ins % problem
    this % ml_ins  => ml_ins

  end subroutine Init_ML_INS_Integrator_3D

  !=============================================================================
  ! ML_INS_IntegratorOptions_3D: type-bound procedures

  !-----------------------------------------------------------------------------
  !> MPI broadcasting of multilevel integrator options

  subroutine Bcast_ML_INS_IntegratorOptions_3D(this, root, comm)
    class(ML_INS_IntegratorOptions_3D), intent(inout) :: this
    integer,        intent(in) :: root !< rank of broadcast root
    type(MPI_Comm), intent(in) :: comm !< MPI communicator
  end subroutine Bcast_ML_INS_IntegratorOptions_3D

  !=============================================================================

end module ML__INS__Integrator__3D
