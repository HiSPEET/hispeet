!> summary:  Base type of one-step ML time integrators for incompressible flows
!> author:   Joerg Stiller
!> date:     2025/05/12
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module ML__INS__Integrator__3D
  use ML__INS__Operator__3D
  implicit none
  private

  public :: ML_INS_Integrator_3D

  !-----------------------------------------------------------------------------
  !> Base type for one-step multilevel IMEX INS integrators

  type, abstract :: ML_INS_Integrator_3D
    type(ML_INS_Operator_3D) :: ml_ins
  end type ML_INS_Integrator_3D

  !=============================================================================

end module ML__INS__Integrator__3D
