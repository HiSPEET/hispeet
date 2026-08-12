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

program Schwarz_Overlap
  use Kind_Parameters
  use Gauss_Jacobi
  implicit none

  real(RNP) :: delta
  integer   :: p

  write(*,'(/,A,2/,2X,A)',advance='NO') 'Schwarz overlap', 'delta = '
  read *, delta

  write(*,'(/,2A5,/)') 'po', 'no'

  do p = 1, 32
    write(*,'(2I5)') p, count(LobattoPoints(p) <= 2 * delta - 1)
  end do

end program Schwarz_Overlap
