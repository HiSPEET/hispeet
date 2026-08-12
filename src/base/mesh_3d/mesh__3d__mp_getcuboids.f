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

!> summary:  Generation of mesh points based on the approximate cuboids
!> author:   Joerg Stiller
!> date:     2021/03/17
!===============================================================================

submodule(Mesh__3D) MP_GetCuboids
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Generates linear mesh points based on the approximate cuboids

  module subroutine GetCuboids(mesh, x)
    class(Mesh_3D),         intent(in)  :: mesh         !< mesh partition
    real(RNP), allocatable, intent(out) :: x(:,:,:,:,:) !< mesh points

    integer :: d, e, i, j, k

    !$omp single
    allocate(x(0:1, 0:1, 0:1, mesh%n_elem, 3))
    !$omp end single

    if (mesh % n_elem < 1) return

    !$omp do
    do e = 1, mesh % n_elem
      associate(x_c => mesh % element(e) % geometry % x_c)

        do d = 1, 3
          do k = 0, 1
          do j = 0, 1
          do i = 0, 1
            x(i,j,k,e,d) = x_c(0,d) + x_c(1,d) * (2*i - 1) &
                                    + x_c(2,d) * (2*j - 1) &
                                    + x_c(3,d) * (2*k - 1)
          end do
          end do
          end do
        end do

      end associate
    end do

  end subroutine GetCuboids

  !=============================================================================

end submodule MP_GetCuboids
