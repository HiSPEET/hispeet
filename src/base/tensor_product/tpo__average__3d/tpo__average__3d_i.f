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

!> summary:  3D generic isotropic averaging operator
!> author:   Joerg Stiller
!> date:     2022/10/02
!===============================================================================

module TPO__Average__3D_I
  use Kind_Parameters, only: RDP
  use TPO__Average__3D_I__Gen
  implicit none
  private

  public :: TPO_Average_I_RDP

contains

  !-----------------------------------------------------------------------------
  !> 3D scaled isotropic averaging operator

  subroutine TPO_Average_I_RDP(w, u, v)
    real(RDP), intent(in)  :: w(:)       !< 1D averaging operator
    real(RDP), intent(in)  :: u(:,:,:,:) !< operand
    real(RDP), intent(out) :: v(:)       !< result

    integer :: np, ne

    np = size(w)
    ne = size(u,4)

    call TPO_Average_I_Gen(np, ne, w, u, v)

  end subroutine TPO_Average_I_RDP

  !=============================================================================

end module TPO__Average__3D_I
