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

!> summary:  Collection of test cases for incompressible single-phase flow
!> author:   Joerg Stiller
!> date:     2018/05/22
!===============================================================================

module ISP_Flow_Problem__Test_Suite

  use ISP_Flow_Problem
  use ISP_Flow_Problem__Channel           , only: FlowProblem_Channel
  use ISP_Flow_Problem__No_Flow           , only: FlowProblem_No_Flow
  use ISP_Flow_Problem__Poiseuille        , only: FlowProblem_Poiseuille
  use ISP_Flow_Problem__Stokes_DKM        , only: FlowProblem_Stokes_DKM
  use ISP_Flow_Problem__Stokes_GMS        , only: FlowProblem_Stokes_GMS
  use ISP_Flow_Problem__Transition_TG     , only: FlowProblem_Transition_TG
  use ISP_Flow_Problem__Vortex_HW         , only: FlowProblem_Vortex_HW
  use ISP_Flow_Problem__Vortex_TG         , only: FlowProblem_Vortex_TG
  use ISP_Flow_Problem__Vortex_Sheet      , only: FlowProblem_VortexSheet
  use ISP_Flow_Problem__Variable_Viscosity, only: FlowProblem_VariableViscosity

end module ISP_Flow_Problem__Test_Suite
