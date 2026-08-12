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

!-------------------------------------------------------------------------------
!> Parametrized 3d AxAxA kernel using hand-crafted suboperators

subroutine PROC(TPO_AxAxA_Hand__,_NA1_,_NA2_)(ne, A, u, v)
  integer,   intent(in)  :: ne                       !< num elements
  real(RWP), intent(in)  :: A(_NA1_,_NA2_)           !< 1D operator
  real(RWP), intent(in)  :: u(_NA2_,_NA2_,_NA2_,ne)  !< operand
  real(RWP), intent(out) :: v(_NA1_,_NA1_,_NA1_,ne)  !< result

  real(RWP), parameter :: alpha = 1
  real(RWP), parameter :: beta  = 0
  real(RWP) :: At(_NA2_,_NA1_)
  real(RWP) :: z2(_NA2_,_NA1_,_NA1_)
  real(RWP) :: z3(_NA2_,_NA2_,_NA1_)

  integer, parameter :: NA2_NA2 = _NA2_ * _NA2_
  integer :: e

  !-----------------------------------------------------------------------------
  ! initialization

  At = transpose(A)

  ! make sure that aux array does not contain NaNs
  z2 = 0
  z3 = 0

  !-----------------------------------------------------------------------------
  ! evaluation

  !$acc data present(u,v) copyin(At)
  !$acc parallel
  !$acc loop gang worker private(z2,z3)

  !$omp do
  do e = 1, ne
    v(:,:,:,e) = 0  ! avoid trouble with NaNs
    call PROC(CtxIxI__,_NA2_,_NA1_)(NA2_NA2, At, alpha, beta, u(:,:,:,e), z3)
    call PROC(IxBtxI__,_NA2_,_NA1_)(_NA2_, _NA1_, At, alpha, beta, z3, z2)
    call PROC(IxIxAt__,_NA2_,_NA1_)(_NA1_, _NA1_, At, alpha, beta, z2, v(:,:,:,e))
  end do

  !$acc end parallel
  !$acc end data

  !-----------------------------------------------------------------------------

end subroutine PROC(TPO_AxAxA_Hand__,_NA1_,_NA2_)
