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
!> Parametrized 3d isotropic spectral kernel using hand-crafted suboperators

subroutine PROC(TPO_Spectral_CI_Hand__,_NP_)(ne, S, Lambda, u, v)
  integer,   intent(in)  :: ne                     !< num elements
  real(RWP), intent(in)  :: S(_NP_,_NP_)           !< 1D eigenvectors
  real(RWP), intent(in)  :: Lambda(_NP_,_NP_,_NP_) !< 3D eigenvalues
  real(RWP), intent(in)  :: u(_NP_,_NP_,_NP_,ne)   !< operand
  real(RWP), intent(out) :: v(_NP_,_NP_,_NP_,ne)   !< result

  real(RWP), parameter :: alpha = 1
  real(RWP), parameter :: beta  = 0
  real(RWP) :: St(_NP_,_NP_)
  real(RWP) :: z (_NP_,_NP_,_NP_)
  integer   :: e

  !---------------------------------------------------------------------------
  ! initialization

  St = transpose(S)
  z  = 0

  !---------------------------------------------------------------------------
  ! evaluation

  !$acc data present(u,v) copyin(S,St,Lambda)
  !$acc parallel
  !$acc loop gang worker private(z)

  !$omp do private(e)
  do e = 1, ne
    v(:,:,:,e) = 0  ! avoid trouble with NaNs
    call PROC(QtxIxI__,_NP_)(S, alpha, beta, u(:,:,:,e), z)  ! z  = SᵀxIxI uᵉ
    call PROC(IxQtxI__,_NP_)(S, alpha, beta, z, v(:,:,:,e))  ! vᵉ = IxSᵀxI z
    call PROC(IxIxQt__,_NP_)(S, alpha, beta, v(:,:,:,e), z)  ! z  = IxIxSᵀ vᵉ
    z = Lambda * z
    call PROC(IxIxQt__,_NP_)(St, alpha, beta, z, v(:,:,:,e)) ! vᵉ = IxIxS z
    call PROC(IxQtxI__,_NP_)(St, alpha, beta, v(:,:,:,e), z) ! z  = IxSxI vᵉ
    call PROC(QtxIxI__,_NP_)(St, alpha, beta, z, v(:,:,:,e)) ! vᵉ = SxIxI z
  end do
  !$omp end do

  !$acc end parallel
  !$acc end data

end subroutine PROC(TPO_Spectral_CI_Hand__,_NP_)
