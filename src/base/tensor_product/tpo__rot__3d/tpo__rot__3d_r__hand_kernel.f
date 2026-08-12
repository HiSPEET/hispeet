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
!> Parametrized 3d rotation kernel using hand-crafted suboperators (RCLI)

subroutine PROC(TPO_Rot_R_Hand__,_NP_)(ne, Ds, dx, u, v)

  use Constants, only: ZERO, ONE

  integer,   intent(in)  :: ne                     !< num elements
  real(RWP), intent(in)  :: Ds(_NP_,_NP_)          !< standard diff matrix
  real(RWP), intent(in)  :: dx(3)                  !< element extensions
  real(RWP), intent(in)  :: u(_NP_,_NP_,_NP_,ne,3) !< 3D vector field
  real(RWP), intent(out) :: v(_NP_,_NP_,_NP_,ne,3) !< element-wise rotation of u

  real(RWP) :: A(_NP_,_NP_)
  real(RWP) :: g(3)

  integer :: e

  !-----------------------------------------------------------------------------
  ! initialization

  A = transpose(Ds)

  ! metric coefficients
  g = 2 / dx

  !---------------------------------------------------------------------------
  ! evaluation

  !$acc data present(u,v) copyin(c,M,Lm)
  !$acc parallel
  !$acc loop gang worker private(M_u)

  !$omp do
  do e = 1, ne

    ! v2 = du1/dx3
    v(:,:,:,e,2) = 0  ! avoid trouble with NaNs
    call PROC(QtxIxI__,_NP_)(A, g(3), ZERO, u(:,:,:,e,1), v(:,:,:,e,2))
    ! v3 =-du1/dx2
    v(:,:,:,e,3) = 0  ! avoid trouble with NaNs
    call PROC(IxQtxI__,_NP_)(A,-g(2), ZERO, u(:,:,:,e,1), v(:,:,:,e,3))
    ! v3+= du2/dx1
    call PROC(IxIxQt__,_NP_)(A, g(1), ONE , u(:,:,:,e,2), v(:,:,:,e,3))
    ! v1 =-du2/dx3
    v(:,:,:,e,1) = 0  ! avoid trouble with NaNs
    call PROC(QtxIxI__,_NP_)(A,-g(3), ZERO, u(:,:,:,e,2), v(:,:,:,e,1))
    ! v1+= du3/dx2
    call PROC(IxQtxI__,_NP_)(A, g(2), ONE , u(:,:,:,e,3), v(:,:,:,e,1))
    ! v2+=-du3/dx1
    call PROC(IxIxQt__,_NP_)(A,-g(1), ONE , u(:,:,:,e,3), v(:,:,:,e,2))

  end do

  !$acc end parallel
  !$acc end data

  !-----------------------------------------------------------------------------

end subroutine PROC(TPO_Rot_R_Hand__,_NP_)
