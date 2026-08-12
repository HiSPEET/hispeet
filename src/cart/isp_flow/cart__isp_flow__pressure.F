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

!> summary:  ISP flow: pressure step
!> author:   Joerg Stiller
!> date:     2018/04/17
!===============================================================================

module CART__ISP_Flow__Pressure

  use Kind_Parameters, only: RNP
  use ISP_Flow_Problem
  use CART__ISP_Flow__Operators
  use CART__ISP_Flow__Pressure__EO
  use CART__ISP_Flow__Pressure__MO

  implicit none
  private

  public :: PressureSolver
  public :: ComputePressure

  interface PressureSolver
    module procedure PressureSolver_IBC
    module procedure PressureSolver_CBC
  end interface

contains

  !-----------------------------------------------------------------------------
  !>  Pressure solver with implied boundary conditions

  subroutine PressureSolver_IBC(problem, flow_op, dt, v_i, p, w, i_max)
    class(FlowProblem),   intent(in)    :: problem        !< flow problem
    class(FlowOperators), intent(inout) :: flow_op        !< flow operators
    real(RNP),            intent(in)    :: dt             !< time-step size
    real(RNP),            intent(in)    :: v_i(:,:,:,:,:) !< v*
    real(RNP),            intent(inout) :: p(:,:,:,:)     !< pressure
    real(RNP),            intent(out)   :: w(:,:,:,:,:)   !< workspace
    integer,    optional, intent(in)    :: i_max          !< num MG/CG cycles

    if (flow_op % po_p == flow_op % po_u) then
      call PressureSolver_EO(problem, flow_op, dt, v_i, p, w, i_max)
    else
      call PressureSolver_MO(problem, flow_op, dt, v_i, p, w, i_max)
    end if

  end subroutine PressureSolver_IBC

  !-----------------------------------------------------------------------------
  !>  Pressure solver with consistent boundary conditions

  subroutine PressureSolver_CBC(problem, flow_op, t, dt, v_i, F_v, p, w, i_max)
    class(FlowProblem),   intent(in)    :: problem         !< flow problem
    class(FlowOperators), intent(inout) :: flow_op         !< flow operators
    real(RNP),            intent(in)    :: t               !< time
    real(RNP),            intent(in)    :: dt              !< time-step size
    real(RNP),            intent(in)    :: v_i(:,:,:,:,:)  !< ṽ
    real(RNP),            intent(in)    :: F_v(:,:,:,:,:)  !< ∂ṽ/∂t
    real(RNP),            intent(inout) :: p(:,:,:,:)      !< pressure
    real(RNP),            intent(out)   :: w(:,:,:,:,:)    !< workspace
    integer,    optional, intent(in)    :: i_max           !< num MG/CG cycles

    if (flow_op % po_p == flow_op % po_u) then
      call PressureSolver_EO(problem, flow_op, t, dt, v_i, F_v, p, w, i_max)
    else
      call PressureSolver_MO(problem, flow_op, t, dt, v_i, F_v, p, w, i_max)
    end if

  end subroutine PressureSolver_CBC

  !-----------------------------------------------------------------------------
  !> Pressure computation

  subroutine ComputePressure(problem, flow_op, t, F_v, p, w, i_max)
    class(FlowProblem),   intent(in)    :: problem         !< flow problem
    class(FlowOperators), intent(inout) :: flow_op         !< flow operators
    real(RNP),            intent(in)    :: t               !< time
    real(RNP),            intent(in)    :: F_v(:,:,:,:,:)  !< ∂v/∂t + ∇p
    real(RNP),            intent(inout) :: p(:,:,:,:)      !< pressure
    real(RNP),            intent(out)   :: w(:,:,:,:,:)    !< workspace
    integer,    optional, intent(in)    :: i_max           !< num MG/CG cycles

    if (flow_op % po_p == flow_op % po_u) then
      call ComputePressure_EO(problem, flow_op, t, F_v, p, w, i_max)
    else
      call ComputePressure_MO(problem, flow_op, t, F_v, p, w, i_max)
    end if

  end subroutine ComputePressure

  !=============================================================================

end module CART__ISP_Flow__Pressure
