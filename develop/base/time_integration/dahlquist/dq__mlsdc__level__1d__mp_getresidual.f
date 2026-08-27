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

!> summary:  Computation of the collocation residual for the time slice
!> author:   Erik Pfister
!> date:     2024/04/18
!===============================================================================

submodule(DQ__MLSDC__Level__1D) MP_GetResidual
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Multi-step collocation residual r = G - L(u)

  module subroutine GetResidual(this, incremental, lambda, dt, G, u, r)
    class(DQ_MLSDC_Level_1D), intent(in) :: this
    logical,      intent(in)  :: incremental
    complex(RNP), intent(in)  :: lambda  !< lambda
    real(RNP),    intent(in)  :: dt      !< size of the time slice
    complex(RNP), intent(in)  :: G(0:,:) !< FAS RHS
    complex(RNP), intent(in)  :: u(0:,:) !< variable
    complex(RNP), intent(out) :: r(0:,:) !< residual

    call this % ApplyOperator(incremental, lambda, dt, u, r)

    ! refinement condition
    r = G - r

  end subroutine GetResidual

  !=============================================================================

end submodule MP_GetResidual
