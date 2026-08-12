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

!> summary:   3d generic interpolation operator for regular 1:8 h-refinement
!> author:    Jörg Stiller
!> date:      2024/07/24
!>
!> This variant runs through the directions in the order 1,2,3
!===============================================================================

!-----------------------------------------------------------------------------
!> Generic 3d interpolation operator for regular 1:8 h-refinement

subroutine TPO_1To8_Gen_RWP(A, u, v)
  real(RWP), contiguous, intent(in)  :: A(:,:,:)   !< 1D operator (na1,na2,2)
  real(RWP), contiguous, intent(in)  :: u(:,:,:,:) !< operand
  real(RWP), contiguous, intent(out) :: v(:,:,:,:) !< result

  !---------------------------------------------------------------------------
  ! local variables

  real(RWP) :: x ( size(A,1), size(A,2), size(A,2), 2 )
  real(RWP) :: y ( size(A,1), size(A,1), size(A,2), 2, 2 )
  real(RWP) :: At( size(A,2), size(A,1), 2 )
  integer   :: na1, na2
  integer   :: e, i, j, k, p, o

  !---------------------------------------------------------------------------
  ! initialization

  na1 = size(A,1)
  na2 = size(A,2)

  At(:,:,1) = transpose(A(:,:,1))
  At(:,:,2) = transpose(A(:,:,2))

  !---------------------------------------------------------------------------
  ! evaluation

  !$omp do private(e,o)
  do e = 1, size(u,4)

    ! x = IxIxA u^e ............................................................

    do k = 1, na2
    do j = 1, na2
    do i = 1, na1
      x(i,j,k,1) = 0
      x(i,j,k,2) = 0
      do p = 1, na2
        x(i,j,k,1) = x(i,j,k,1) + A(i,p,1) * u(p,j,k,e)
        x(i,j,k,2) = x(i,j,k,2) + A(i,p,2) * u(p,j,k,e)
      end do
    end do
    end do
    end do

    ! y = IxAxI x ..............................................................

    do k = 1, na2
    do j = 1, na1
    do i = 1, na1
      y(i,j,k,1:2,1:2) = 0
      do p = 1, na2
        y(i,j,k,1,1) = y(i,j,k,1,1) + x(i,p,k,1) * At(p,j,1)
        y(i,j,k,2,1) = y(i,j,k,2,1) + x(i,p,k,2) * At(p,j,1)
        y(i,j,k,1,2) = y(i,j,k,1,2) + x(i,p,k,1) * At(p,j,2)
        y(i,j,k,2,2) = y(i,j,k,2,2) + x(i,p,k,2) * At(p,j,2)
      end do
    end do
    end do
    end do

    ! v^e = AxIxI z2 ...........................................................

    ! offset of fine element index
    o = 8 * (e-1)

    do k = 1, na1
    do j = 1, na1
    do i = 1, na1
      v(i,j,k,o+1:o+8) = 0
      do p = 1, na2                                            ! position
        v(i,j,k,o+1) = v(i,j,k,o+1) + y(i,j,p,1,1) * At(p,k,1) !  1,1,1
        v(i,j,k,o+2) = v(i,j,k,o+2) + y(i,j,p,2,1) * At(p,k,1) !  2,1,1
        v(i,j,k,o+3) = v(i,j,k,o+3) + y(i,j,p,1,2) * At(p,k,1) !  1,2,1
        v(i,j,k,o+4) = v(i,j,k,o+4) + y(i,j,p,2,2) * At(p,k,1) !  2,2,1
        v(i,j,k,o+5) = v(i,j,k,o+5) + y(i,j,p,1,1) * At(p,k,2) !  1,1,2
        v(i,j,k,o+6) = v(i,j,k,o+6) + y(i,j,p,2,1) * At(p,k,2) !  2,1,2
        v(i,j,k,o+7) = v(i,j,k,o+7) + y(i,j,p,1,2) * At(p,k,2) !  1,2,2
        v(i,j,k,o+8) = v(i,j,k,o+8) + y(i,j,p,2,2) * At(p,k,2) !  2,2,2
      end do
    end do
    end do
    end do

  end do
  !$omp end do

end subroutine TPO_1To8_Gen_RWP
