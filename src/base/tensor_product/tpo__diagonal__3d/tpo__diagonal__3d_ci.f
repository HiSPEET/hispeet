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

!> summary:  3D scaled constant isotropic diagonal operator
!> author:   Joerg Stiller
!> date:     2020/06/04
!===============================================================================

module TPO__Diagonal__3D_CI
  use Kind_Parameters, only: RDP
  use TPO__Diagonal__3D_CI__Gen
  implicit none
  private

  public :: TPO_Diagonal_CI_RDP

contains

  !-----------------------------------------------------------------------------
  !> 3D scaled constant isotropic diagonal operator, v = s DxDxD u

  subroutine TPO_Diagonal_CI_RDP(s, D, u, v)
    real(RDP), intent(in)  :: s          !< scaling factor
    real(RDP), intent(in)  :: D(:)       !< diagonal 1D operator
    real(RDP), intent(in)  :: u(:,:,:,:) !< operand
    real(RDP), intent(out) :: v(:,:,:,:) !< result

    integer :: np, ne

    np = size(D)
    ne = size(u,4)

    call TPO_Diagonal_CI_Gen(np, ne, s, D, u, v)

  end subroutine TPO_Diagonal_CI_RDP

  !=============================================================================

end module TPO__Diagonal__3D_CI
