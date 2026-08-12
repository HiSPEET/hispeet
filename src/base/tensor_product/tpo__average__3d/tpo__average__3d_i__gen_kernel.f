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

!-------------------------------------------------------------------------------
!> 3D generic isotropic averaging operator

subroutine TPO_Average_I_Gen_RWP(np, ne, w, u, v)
  integer,   intent(in)  :: np             !< number of points per direction
  integer,   intent(in)  :: ne             !< number of elements
  real(RWP), intent(in)  :: w(np)          !< 1D averaging operator
  real(RWP), intent(in)  :: u(np,np,np,ne) !< operand
  real(RWP), intent(out) :: v(ne)          !< elementwise average

  real(RWP) :: s
  integer   :: e, i, j, k

  s = 1 / sum(w)**3

  !$omp do
  do e = 1, ne
    v(e) = 0
    do k = 1, np
    do j = 1, np
    do i = 1, np
      v(e) = v(e) + w(k) * w(j) * w(i) * u(i,j,k,e)
    end do
    end do
    end do
    v(e) = s * v(e)
  end do

end subroutine TPO_Average_I_Gen_RWP
