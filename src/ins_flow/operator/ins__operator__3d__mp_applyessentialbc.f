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

!> summary:  Application of essential boundary conditions
!> author:   Joerg Stiller
!> date:     2022/12/06
!===============================================================================

submodule(INS__Operator__3D) MP_ApplyEssentialBC
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Application of essential boundary conditions to velocity trace variables
  !>
  !> Construction of the outer boundary traces of velocity using inner traces
  !> and boundary values. In the absence of the latter, homogeneous conditions
  !> are assumed.

  module subroutine ApplyEssentialBC(this, bv_u, vm, vp)
    class(INS_Operator_3D), intent(in) :: this
    !< Navier-Stokes operator
    class(BoundaryVariable_3D), optional, intent(in) :: bv_u(:)
    !< boundary values depending on BC type
    real(RNP), contiguous, intent(in) :: vm(:,:,:,:,:)
    !< inner velocity traces, `vm(np,np,6,ne,3) = v⁻`
    real(RNP), contiguous, intent(inout) :: vp(:,:,:,:,:)
    !< outer velocity traces, `vp(np,np,6,ne,3) = v⁺` on boundary faces,
    !! values on interior faces remain unchanged

    integer :: b, e, f, m

    do b = 1, this % mesh % n_bound
      associate(boundary => this % mesh % boundary(b))

        select case(this % problem % bc_v(b))

        case('D') ! Dirichlet: v⁺ = 2vᵇ - v⁻

          !$omp do
          do f = 1, boundary % n_face
            e = boundary % face(f) % element_id
            m = boundary % face(f) % element_face
            if (present(bv_u)) then
              associate(vb => bv_u(b) % val)
                vp(:,:,m,e,1:3) = 2 * vb(:,:,f,1:3) - vm(:,:,m,e,1:3)
              end associate
            else
              vp(:,:,m,e,1:3) = -vm(:,:,m,e,1:3)
            end if
          end do

        case('O') ! Outflow: v⁺ = v⁻

          !$omp do
          do f = 1, boundary % n_face
            e = boundary % face(f) % element_id
            m = boundary % face(f) % element_face
            vp(:,:,m,e,1:3) = vm(:,:,m,e,1:3)
          end do

        end select

      end associate
    end do

  end subroutine ApplyEssentialBC

  !=============================================================================

end submodule MP_ApplyEssentialBC
