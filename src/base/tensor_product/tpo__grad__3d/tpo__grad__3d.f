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

!> summary:  3D gradient operator
!> author:   Joerg Stiller
!> date:     2022/10/09
!===============================================================================

module TPO__Grad__3D
  use Kind_Parameters, only: RDP
  use Standard_Element_Operators__1D
  use Mesh__3D
  use Spectral_Element_Mesh__3D

  use TPO__Grad__3D_R ! TPO for regular meshes
  use TPO__Grad__3D_D ! TPO for deformed meshes

  private

  public :: TPO_Grad

  interface TPO_Grad
    module procedure TPO_Grad_RDP
  end interface

contains

  !-----------------------------------------------------------------------------
  !> Automatic gradient TPO with double precision

  subroutine TPO_Grad_RDP(eop, sem, u, up, v)
    class(StandardElementOperators_1D), intent(in) :: eop
    class(SpectralElementMesh_3D), intent(in) :: sem
    real(RDP),           intent(in)  :: u(:,:,:,:)   !< 3D scalar field u
    real(RDP), optional, intent(in)  :: up(:,:,:,:)  !< exterior traces u⁺
    real(RDP),           intent(out) :: v(:,:,:,:,:) !< gradient of u

    if (sem % mesh % regular) then
      call TPO_Grad_R( eop % w, eop % D, sem % mesh % dx, u, up, v )
    else
      call TPO_Grad_D( eop % w            &
                     , eop % D            &
                     , sem % metrics % Jd &
                     , sem % metrics % Ji &
                     , sem % metrics % a  &
                     , sem % metrics % n  &
                     , u, up, v           )
    end if

  end subroutine TPO_Grad_RDP

  !=============================================================================

end module TPO__Grad__3D
