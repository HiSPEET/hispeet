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

!> summary:  Provision of variable viscosity coefficients
!> author:   Joerg Stiller
!> date:     2025/05/21
!===============================================================================

submodule(INS__Operator__3D) MP_GetVariableViscosity
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Variable viscosity coefficients

  module subroutine GetVariableViscosity(this, t, u, mu, nu)
    class(INS_Operator_3D), intent(in)  :: this
    real(RNP),              intent(in)  :: t
    real(RNP),  contiguous, intent(in)  :: u(:,:,:,:,:)
    real(RNP),  contiguous, intent(out) :: mu(:,:,:,:)
    real(RNP),  contiguous, intent(out) :: nu(:,:,:,:)

    associate(x => this % sem_u % metrics % x)

      ! physical shear viscosity ...............................................

      if (this % problem % HasVariableProperties()) then
        call this % problem % GetViscosity(x, t, u, nu)
      else
        call SetArray(nu, this % nu_0)
      end if

      ! subgrid scale viscosity - using mu as workspace ........................

      if (this % sgs_model % model > 0) then
        call this % sgs_model % Get_SGS_Viscosity(this % sem_u, u, mu)
        call MergeArrays(ONE, nu, ONE, mu)
      end if

      ! artificial bulk viscosity ..............................................

      if (this % c_mu > 0) then
        call GetVariableBulkViscosity(this, u, nu, mu)
      else
        call SetArray(mu, this % mu_0)
      end if

    end associate

  end subroutine GetVariableViscosity

  !-----------------------------------------------------------------------------
  !> Variable bulk viscosity
  !>
  !> Element by element computation of the stabilizing bulk viscosity defined by
  !>
  !>      μᵉ = c_μ sqrt(ν² + (vh)²)
  !>
  !> where
  !>
  !>      ν = max νᵉ
  !>      v = avg|vᵉ|
  !>      h = ∆xᵉ / P

  subroutine GetVariableBulkViscosity(this, u, nu, mu)
    class(INS_Operator_3D), intent(in)  :: this
    real(RNP),  contiguous, intent(in)  :: u(:,:,:,:,:)
    real(RNP),  contiguous, intent(in)  :: nu(:,:,:,:)
    real(RNP),  contiguous, intent(out) :: mu(:,:,:,:)

    real(RNP) :: c, dx(3), h_e, nu_e, v_e
    integer :: e

    ! normalization factor
    c = ONE / size(u,1)**3

    !$omp do
    do e = 1, size(mu,4)
      call this % mesh % element(e) % GetCuboidDimensions(dx)
      h_e  = product(dx)**THIRD / this%eop_u%po
      v_e  = c * sum(sqrt(u(:,:,:,e,1)**2 + u(:,:,:,e,2)**2 + u(:,:,:,e,3)**2))
      nu_e = c * sum(nu(:,:,:,e))
      mu(:,:,:,e) = this%c_mu * max(nu_e, v_e * h_e)
    end do

  end subroutine GetVariableBulkViscosity

  !=============================================================================

end submodule MP_GetVariableViscosity
