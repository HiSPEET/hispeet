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
!> Parametrized 3d rotrot kernel using hand-crafted suboperators (RCLI)

subroutine PROC(TPO_RotRot_R_Hand__,_NP_)(ne, Ds, dx, u, v)

  use Constants, only: ZERO, ONE

  integer,   intent(in)  :: ne                     !< num elements
  real(RWP), intent(in)  :: Ds(_NP_,_NP_)          !< standard diff matrix
  real(RWP), intent(in)  :: dx(3)                  !< element extensions
  real(RWP), intent(in)  :: u(_NP_,_NP_,_NP_,ne,3) !< 3D vector field
  real(RWP), intent(out) :: v(_NP_,_NP_,_NP_,ne,3) !< v = rot(rot(u))

  real(RWP) :: A(_NP_,_NP_), DA(_NP_,_NP_)
  real(RWP) :: div_u(_NP_,_NP_,_NP_)
  real(RWP) :: g(3)

  integer :: e, i, j

  !-----------------------------------------------------------------------------
  ! initialization

  A  = transpose(Ds)

  do j = 1, _NP_
  do i = 1, _NP_
    DA(j,i) = sum(Ds(i,:) * Ds(:,j))
  end do
  end do

  ! metric coefficients
  g = 2 / dx

  ! make sure that aux array doe not contain NaNs
  div_u = 0

  !---------------------------------------------------------------------------
  ! evaluation

  !$acc data present(u,v) copyin(c,M,Lm)
  !$acc parallel
  !$acc loop gang worker private(M_u)

  !$omp do
  do e = 1, ne

    ! div(u) = du1/dx1
    call PROC(IxIxQt__,_NP_)( A, g(1)     , ZERO, u(:,:,:,e,1), div_u(:,:,:))
    ! div(u)+= du2/dx2
    call PROC(IxQtxI__,_NP_)( A, g(2)     , ONE , u(:,:,:,e,2), div_u(:,:,:))
    ! div(u)+= du3/dx3
    call PROC(QtxIxI__,_NP_)( A, g(3)     , ONE , u(:,:,:,e,3), div_u(:,:,:))

    ! v1 = d(div(u))/dx1
    v(:,:,:,e,1) = 0  ! avoid trouble with NaNs
    call PROC(IxIxQt__,_NP_)( A, g(1)     , ZERO, div_u(:,:,:), v(:,:,:,e,1))
    ! v1+=-d( du1/dx1 )/dx1
    call PROC(IxIxQt__,_NP_)(DA,-g(1)*g(1), ONE , u(:,:,:,e,1), v(:,:,:,e,1))
    ! v1+=-d( du1/dx2 )/dx2
    call PROC(IxQtxI__,_NP_)(DA,-g(2)*g(2), ONE , u(:,:,:,e,1), v(:,:,:,e,1))
    ! v1+=-d( du1/dx3 )/dx3
    call PROC(QtxIxI__,_NP_)(DA,-g(3)*g(3), ONE , u(:,:,:,e,1), v(:,:,:,e,1))

    ! v2 = d(div(u))/dx2
    v(:,:,:,e,2) = 0  ! avoid trouble with NaNs
    call PROC(IxQtxI__,_NP_)( A, g(2)     , ZERO, div_u(:,:,:), v(:,:,:,e,2))
    ! v2+=-d( du2/dx1 )/dx1
    call PROC(IxIxQt__,_NP_)(DA,-g(1)*g(1), ONE , u(:,:,:,e,2), v(:,:,:,e,2))
    ! v2+=-d( du2/dx2 )/dx2
    call PROC(IxQtxI__,_NP_)(DA,-g(2)*g(2), ONE , u(:,:,:,e,2), v(:,:,:,e,2))
    ! v2+=-d( du2/dx3 )/dx3
    call PROC(QtxIxI__,_NP_)(DA,-g(3)*g(3), ONE , u(:,:,:,e,2), v(:,:,:,e,2))

    ! v3 = d(div(u))/dx3
    v(:,:,:,e,3) = 0  ! avoid trouble with NaNs
    call PROC(QtxIxI__,_NP_)( A, g(3)     , ZERO, div_u(:,:,:), v(:,:,:,e,3))
    ! v3+=-d( du3/dx1 )/dx1
    call PROC(IxIxQt__,_NP_)(DA,-g(1)*g(1), ONE , u(:,:,:,e,3), v(:,:,:,e,3))
    ! v3+=-d( du3/dx2 )/dx2
    call PROC(IxQtxI__,_NP_)(DA,-g(2)*g(2), ONE , u(:,:,:,e,3), v(:,:,:,e,3))
    ! v3+=-d( du3/dx3 )/dx3
    call PROC(QtxIxI__,_NP_)(DA,-g(3)*g(3), ONE , u(:,:,:,e,3), v(:,:,:,e,3))

  end do

  !$acc end parallel
  !$acc end data

  !-----------------------------------------------------------------------------

end subroutine PROC(TPO_RotRot_R_Hand__,_NP_)
