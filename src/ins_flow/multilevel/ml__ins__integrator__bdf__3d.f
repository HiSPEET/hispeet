!> summary:  Multilevel IMEX BDF method for incompressible flows
!> author:   Joerg Stiller
!> date:     2025/05/12
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module ML__INS__Integrator__BDF__3D
  use ML__INS__Integrator__3D
  implicit none
  private

  public :: ML_INS_BDF_3D

  !=============================================================================
  ! Types

  !-----------------------------------------------------------------------------
  !> Type providing IMEX BDF solvers for 3D incompressible flows

  type, extends(ML_INS_Integrator_3D) :: ML_INS_BDF_3D
    integer :: order
  contains
    !procedure ::
  end type ML_INS_BDF_3D

  !=============================================================================

end module ML__INS__Integrator__BDF__3D
