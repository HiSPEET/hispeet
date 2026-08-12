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

!> summary:  3D Cartesian gradient operator
!> author:   Erik Pfister, Joerg Stiller
!> date:     2020/03/24
!===============================================================================

module TPO__Grad__3D_R
  use Kind_Parameters, only: RDP
  use TPO__Grad__3D_R__Hand
  use TPO__Grad__3D_R__XSMM
  implicit none
  private

  public :: TPO_Grad_R

  interface TPO_Grad_R
    module procedure TPO_Grad_R_RDP
  end interface

contains

  !-----------------------------------------------------------------------------
  ! Gradient of a scalar field: 3D Cartesian equidistant (Lobatto)

  subroutine TPO_Grad_R_RDP(Ms, Ds, dx, u, up, v)
    real(RDP), intent(in)  :: Ms(:)         !< 1D standard mass matrix
    real(RDP), intent(in)  :: Ds(:,:)       !< 1D standard stiffness matrix
    real(RDP), intent(in)  :: dx(3)         !< element extensions
    real(RDP), intent(in)  :: u(:,:,:,:)    !< 3D scalar Field
    real(RDP), intent(in)  :: up(:,:,:,:)   !< exterior traces u⁺ at faces
    real(RDP), intent(out) :: v(:,:,:,:,:)  !< element-wise gradient of u

    optional :: Ms, up

    ! local variables ..........................................................

    integer :: np, ne

    ! initialization ...........................................................

    np = size(u,1)
    ne = size(u,4)

    ! evaluation ...............................................................

    select case(np)

#ifdef __LIBXSMM__
    case(4:)
      call TPO_Grad_R_XSMM(np, ne, Ms, Ds, dx, u, up, v)
#endif

    case default
      call TPO_Grad_R_Hand(np, ne, Ms, Ds, dx, u, up, v)

    end select

  end subroutine TPO_Grad_R_RDP

  !=============================================================================

end module TPO__Grad__3D_R
