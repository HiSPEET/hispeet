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
  use Smooth_Mesh_Data__3D
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

      if (this % div_stab > 1) then
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
  !>      μᵉ = c_μ max(ν, vh)
  !>
  !> where ν is the kinematic shear viscosity, v the characteristic velocity and
  !> h = ∆xᵉ/P the mean element node spacing.

  subroutine GetVariableBulkViscosity(this, u, nu, mu)
    class(INS_Operator_3D), intent(in)  :: this
    real(RNP),  contiguous, intent(in)  :: u(:,:,:,:,:)
    real(RNP),  contiguous, intent(in)  :: nu(:,:,:,:)
    real(RNP),  contiguous, intent(out) :: mu(:,:,:,:)

    real(RNP) :: c, dx(3), h_e, nu_e, v_e
    logical   :: local
    integer   :: e

    local = any(this%div_stab == [3,5])

    if (local) then
      ! normalization factor
      c = ONE / size(u,1)**3
    else
      ! use reference values
      nu_e = this % problem % nu_ref
      v_e  = this % problem % v_ref
    end if

    !$omp do
    do e = 1, size(mu,4)
      call this % mesh % element(e) % GetCuboidDimensions(dx)
      h_e = product(dx)**THIRD / this%eop_u%po
      if (local) then
        nu_e = c * sum(nu(:,:,:,e))
        v_e  = c * sum( sqrt( u(:,:,:,e,1)**2   &
                            + u(:,:,:,e,2)**2   &
                            + u(:,:,:,e,3)**2 ) )
      end if
      mu(:,:,:,e) = this%mu_0 * max(nu_e, v_e * h_e)
    end do

    ! optional filtering
    if (any(this%div_stab == [4,5])) then
      call SmoothMeshData_3D( mesh   = this % mesh  &
                            , eop    = this % eop_u &
                            , u      = mu           &
                            , filter = 1            & ! bubble cut-off
                            , order  = 1            ) ! retain multilinear part
    end if

  end subroutine GetVariableBulkViscosity

  !=============================================================================

end submodule MP_GetVariableViscosity
