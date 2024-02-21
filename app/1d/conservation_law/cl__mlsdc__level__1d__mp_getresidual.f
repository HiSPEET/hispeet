!> summary:  Computation of the collocation residual for the time slice
!> author:   Erik Pfister
!> date:     2023/09/04
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(CL__MLSDC__Level__1D) MP_GetResidual
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Multi-step collocation residual r = f^ - L(u)

  module subroutine GetResidual(this, G, v, r)
    class(CL_MLSDC_Level_1D), intent(in) :: this
    real(RNP), optional, intent(in)  :: G(0:,:,:,0:,:) !< FAS defect correction
    real(RNP), intent(in)            :: v(0:,:,:,0:,:) !< u after operator applied
    real(RNP), intent(out)           :: r(0:,:,:,0:,:) !< residual

    if(present(G)) then
      r = G - v
    else
      r = - v
    end if

  end subroutine GetResidual

  !=============================================================================

end submodule MP_GetResidual
