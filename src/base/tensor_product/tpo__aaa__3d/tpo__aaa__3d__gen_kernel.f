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

!> summary:   3d generic AxAxA operator
!> author:    Jörg Stiller
!> date:      2020/05/07
!===============================================================================

!-----------------------------------------------------------------------------
!> Generic AxAxA operator

subroutine TPO_AAA_Gen_RWP(A, u, v)
  real(RWP), intent(in)  :: A(:,:)     !< 1D operator
  real(RWP), intent(in)  :: u(:,:,:,:) !< operand
  real(RWP), intent(out) :: v(:,:,:,:) !< result

  !---------------------------------------------------------------------------
  ! local variables

  real(RWP) :: At( size(A,2), size(A,1) )
  real(RWP) :: z2( size(A,2), size(A,1), size(A,1))
  real(RWP) :: z3( size(A,2), size(A,2), size(A,1))
  real(RWP) :: tmp
  integer   :: na1, na2, ne
  integer   :: e, i, j, k, p

  !---------------------------------------------------------------------------
  ! initialization

  na1 = size(A,1)
  na2 = size(A,2)
  ne  = size(u,4)

  At = transpose(A)

  !---------------------------------------------------------------------------
  ! evaluation

  !$acc data present(u,v) copyin(At)
  !$acc parallel
  !$acc loop gang worker private(z2,z3)

  !$omp do private(e)
  do e = 1, ne

    ! z3 = AxIxI u^e .........................................................

    !$acc loop collapse(3) vector
    do k = 1, na1
    do j = 1, na2
    !DIR$ SIMD
    do i = 1, na2
      tmp = 0
      do p = 1, na2
        tmp = tmp + At(p,k) * u(i,j,p,e)
      end do
      z3(i,j,k) = tmp
    end do
    end do
    end do

    ! z2 = IxAxI z3 ..........................................................

    !$acc loop collapse(3) vector
    do k = 1, na1
    do j = 1, na1
    !DIR$ SIMD
    do i = 1, na2
      tmp = 0
      do p = 1, na2
        tmp = tmp + At(p,j) * z3(i,p,k)
      end do
      z2(i,j,k) = tmp
    end do
    end do
    end do

    ! v^e = IxIxA z2 .........................................................

    !$acc loop collapse(3) vector
    do k = 1, na1
    do j = 1, na1
    !DIR$ SIMD
    do i = 1, na1
      tmp = 0
      do p = 1, na2
        tmp = tmp + At(p,i) * z2(p,j,k)
      end do
      v(i,j,k,e) = tmp
    end do
    end do
    end do

  end do
  !$omp end do

  !$acc end parallel
  !$acc end data

end subroutine TPO_AAA_Gen_RWP
