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

!> summary:  3D Cartesian elliptic operator with constant isotropic diffusivity
!> author:   Joerg Stiller
!> date:     2020/06/03
!===============================================================================

module TPO__Elliptic__3D_RLCI
  use Kind_Parameters, only: RDP
  use TPO__Elliptic__3D_RLCI__Hand
  use TPO__Elliptic__3D_RLCI__XSMM
  implicit none
  private

  public :: TPO_Elliptic_RLCI

  interface TPO_Elliptic_RLCI
    module procedure TPO_Elliptic_RLCI_RDP
  end interface

contains

  !-----------------------------------------------------------------------------
  !> 3D element elliptic operator for scalars
  !>
  !> The procedure evaluates the elliptic operator with constant coefficients
  !> lambda and nu for equidistant 3D Cartesian elements. Given the scalar
  !> variable u(1:np,1:np,1:np,1:ne), the result is computed as
  !>
  !>     v(i,j,k,e)
  !>
  !>       = lambda * dx(1)*dx(2)*dx(3)/8 * Ms(i) * Ms(j) * Ms(k) * u(i,j,k,e)
  !>
  !>       + nu * dx(2)*dx(3)/(2*dx(1)) * Ms(j) * Ms(k)
  !>                                    * sum( Ms(:) * Ds(:,i) * u(:,j,k,e) )
  !>
  !>       + nu * dx(1)*dx(3)/(2*dx(2)) * Ms(i) * Ms(k)
  !>                                    * sum( Ms(:) * Ds(:,j) * u(i,:,k,e) )
  !>
  !>       + nu * dx(2)*dx(2)/(2*dx(3)) * Ms(i) * Ms(j)
  !>                                    * sum( Ms(:) * Ds(:,k) * u(i,j,:,e) )
  !>
  !> for each point (i,j,k) in every element e.
  !> Ms and Ls are the 1D standard mass and stiffness matrices for polynomial
  !> order P = np-1, respectively, and dx(1:3) the element dimensions.

  subroutine TPO_Elliptic_RLCI_RDP(Ms, Ls, dx, lambda, nu, u, v, Ds, ub, qb)
    real(RDP), intent(in)  :: Ms(:)      !< 1D standard mass matrix
    real(RDP), intent(in)  :: Ls(:,:)    !< 1D standard stiffness matrix
    real(RDP), intent(in)  :: dx(3)      !< element extensions
    real(RDP), intent(in)  :: lambda     !< Helmholtz parameter
    real(RDP), intent(in)  :: nu         !< diffusivity
    real(RDP), intent(in)  :: u(:,:,:,:) !< operand
    real(RDP), intent(out) :: v(:,:,:,:) !< result

    real(RDP), optional, intent(in)    :: Ds(:,:)     !< 1D standard diff matrix
    real(RDP), optional, intent(inout) :: ub(:,:,:,:) !< element boundary values
    real(RDP), optional, intent(inout) :: qb(:,:,:,:) !< element boundary fluxes

    contiguous :: Ms, Ls, u, v, Ds, ub, qb

    select case(size(Ms))

#ifdef __LIBXSMM__
    case(5:)
      call TPO_Elliptic_RLCI_XSMM(Ms, Ls, dx, lambda, nu, u, v, Ds, ub, qb)
#endif

    case default
      call TPO_Elliptic_RLCI_Hand(Ms, Ls, dx, lambda, nu, u, v, Ds, ub, qb)

    end select

  end subroutine TPO_Elliptic_RLCI_RDP

  !=============================================================================

end module TPO__Elliptic__3D_RLCI
