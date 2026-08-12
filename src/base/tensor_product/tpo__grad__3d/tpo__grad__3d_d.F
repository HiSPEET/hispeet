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

!> summary:   3D Cartesian gradient operator (D)
!> author:    Jerome Michel, Erik Pfister
!> date:      2021/10/12
!===============================================================================

module TPO__Grad__3D_D
  use Kind_Parameters, only: RDP
  use TPO__Grad__3D_D__Hand
  use TPO__Grad__3D_D__XSMM
  implicit none
  private

  public :: TPO_Grad_D

  interface TPO_Grad_D
    module procedure TPO_Grad_D_RDP
  end interface

contains

  !-----------------------------------------------------------------------------
  ! gradient of a vector field: 3D unstructured and deformed

  subroutine TPO_Grad_D_RDP(Ms, Ds, Jd, Ji, a, n, u, up, v)
    real(RDP), intent(in)  :: Ms(:)           !< 1D standard mass matrix
    real(RDP), intent(in)  :: Ds(:,:)         !< 1D standard diff matrix
    real(RDP), intent(in)  :: Jd(:,:,:,:)     !< element Jacobian determinant
    real(RDP), intent(in)  :: Ji(:,:,:,:,:,:) !< inverse Jacobi matrix
    real(RDP), intent(in)  :: a(:,:,:,:)      !< face area coefficients
    real(RDP), intent(in)  :: n(:,:,:,:,:)    !< face unit normal vectors
    real(RDP), intent(in)  :: u(:,:,:,:)      !< 3D scalar field u
    real(RDP), intent(in)  :: up(:,:,:,:)     !< exterior traces u⁺ at faces
    real(RDP), intent(out) :: v(:,:,:,:,:)    !< element-wise gradient of u

    optional :: Ms, Jd, a, n, up

    ! local variables ..........................................................

    integer :: np, ne

    ! initialization ...........................................................

    np = size(u,1)
    ne = size(u,4)

    ! evaluation ...............................................................

    select case(np)

#ifdef __LIBXSMM__
    case(5:)
      call TPO_Grad_D_XSMM(np, ne, Ms, Ds, Jd, Ji, a, n, u, up, v)
#endif

    case default
      call TPO_Grad_D_Hand(np, ne, Ms, Ds, Jd, Ji, a, n, u, up, v)
    end select

  end subroutine TPO_Grad_D_RDP

  !=============================================================================

end module TPO__Grad__3D_D
