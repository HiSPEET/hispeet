!> summary:  Collection of test cases for incompressible single-phase flow
!> author:   Joerg Stiller
!> date:     2018/05/22
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module ISP_Flow_Problem__Test_Suite

  use ISP_Flow_Problem
  use ISP_Flow_Problem__Stokes_DKM        , only: FlowProblem_Stokes_DKM
  use ISP_Flow_Problem__Stokes_GMS        , only: FlowProblem_Stokes_GMS
  use ISP_Flow_Problem__Vortex_HW         , only: FlowProblem_Vortex_HW
  use ISP_Flow_Problem__Vortex_TG         , only: FlowProblem_Vortex_TG
  use ISP_Flow_Problem__Vortex_Sheet      , only: FlowProblem_VortexSheet
  use ISP_Flow_Problem__Variable_Viscosity, only: FlowProblem_VariableViscosity

end module ISP_Flow_Problem__Test_Suite
