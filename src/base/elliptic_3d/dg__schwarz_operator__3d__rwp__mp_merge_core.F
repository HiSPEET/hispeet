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

  !-----------------------------------------------------------------------------
  !>  Merge of subdomain core contributions into mesh variable

  module subroutine Merge_Core(this, mesh, us, u)
    class(DG_SchwarzOperator_3D), intent(in) :: this
    class(Mesh_3D), intent(in)    :: mesh          !< mesh partition
    real(RWP),      intent(inout) :: us(:,:,:,:)   !< subdomain solutions
    real(RNP),      intent(inout) :: u(0:,0:,0:,:) !< mesh variable

    integer :: e, i, j, k
    integer :: c0

    associate(po => this%po, no => this%no)

      ! index of first core point within subdomain
      c0 = no + 1

      !$omp do
      do e = 1, mesh % n_elem_active
        do k = 0, po
        do j = 0, po
        do i = 0, po
          u(i,j,k,e) = u(i,j,k,e) + us(c0 + i, c0 + j, c0 + k, e)
        end do
        end do
        end do
      end do

    end associate

  end subroutine Merge_Core
