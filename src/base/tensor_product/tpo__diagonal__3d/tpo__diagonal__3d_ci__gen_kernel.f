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

!> summary:  3D generic scaled constant isotropic diagonal operator
!> author:   Joerg Stiller
!> date:     2020/06/04
!===============================================================================

!-----------------------------------------------------------------------------
!> 3D generic scaled constant isotropic diagonal operator, v = s DxDxD u

subroutine TPO_Diagonal_CI_Gen_RWP(np, ne, s, D, u, v)
  integer,   intent(in)  :: np             !< number of points per direction
  integer,   intent(in)  :: ne             !< number of elements
  real(RWP), intent(in)  :: s              !< scaling factor
  real(RWP), intent(in)  :: D(np)          !< diagonal 1D operator
  real(RWP), intent(in)  :: u(np,np,np,ne) !< collapsed 3D operand
  real(RWP), intent(out) :: v(np,np,np,ne) !< collapsed 3D result

  integer :: e, i, j, k

  !$acc data present(u,v) copyin(D) async
  !$acc parallel loop collapse(4) async
  !$omp do
  do e = 1, ne
    do k = 1, np
    do j = 1, np
    do i = 1, np
      v(i,j,k,e) = s * D(k) * D(j) * D(i) * u(i,j,k,e)
    end do
    end do
    end do
  end do
  !$omp end do
  !$acc end parallel
  !$acc end data

end subroutine TPO_Diagonal_CI_Gen_RWP
