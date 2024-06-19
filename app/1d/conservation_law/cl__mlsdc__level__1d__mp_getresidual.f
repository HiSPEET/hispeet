!> summary:  Computation of the collocation residual for the time slice
!> author:   Erik Pfister
!> date:     2023/09/04
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(CL__MLSDC__Level__1D) MP_GetResidual
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Multi-step collocation residual r = g - L(u)

  module subroutine GetResidual(this, dt, t_0, G, u, r)
    class(CL_MLSDC_Level_1D), intent(in) :: this
    real(RNP), intent(in)  :: dt             !< size of the time slice
    real(RNP), intent(in)  :: t_0            !< start time of the slice
    real(RNP), intent(in)  :: g(0:,:,:,0:,:) !< FAS RHS
    real(RNP), intent(in)  :: u(0:,:,:,0:,:) !< variable
    real(RNP), intent(out) :: r(0:,:,:,0:,:) !< residual

    call this % ApplyOperator(dt, t_0, u, r)

    ! refinement condition
    r = g - r

  end subroutine GetResidual

  !=============================================================================

end submodule MP_GetResidual
