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
!> Parametrized 3d isotropic Schwarz kernel using hand-crafted suboperators

subroutine PROC(TPO_Schwarz_CI_Hand__,_NP_)(nc, nd, S, W, cfg, D_inv, f, u)
  integer,   intent(in)  :: nc                       !< num configurations
  integer,   intent(in)  :: nd                       !< num subdomains
  real(RWP), intent(in)  :: S(_NP_,_NP_,nc)          !< 1D eigenvectors
  real(RWP), intent(in)  :: W(_NP_,nc)               !< 1D weights
  integer,   intent(in)  :: cfg(:,:)                 !< subdomain configurations
  real(RWP), intent(in)  :: D_inv(_NP_,_NP_,_NP_,nd) !< inverse 3D eigenvalues
  real(RWP), intent(in)  :: f(_NP_,_NP_,_NP_,nd)     !< RHS
  real(RWP), intent(out) :: u(_NP_,_NP_,_NP_,nd)     !< solution

  real(RWP), parameter :: alpha = 1
  real(RWP), parameter :: beta  = 0

  real(RWP) :: WS_t ( _NP_, _NP_,  nc  )
  real(RWP) :: y    ( _NP_, _NP_, _NP_ )
  real(RWP) :: z    ( _NP_, _NP_, _NP_ )

  integer :: i, j, l
  integer :: c, c1, c2, c3

  !---------------------------------------------------------------------------
  ! initialization

  ! WS_t(:,:,c) = (WS)ᵀ(:,:,c)
  do c = 1, nc
  do j = 1, _NP_
  do i = 1, _NP_
    WS_t(j,i,c) = W(i,c) * S(i,j,c)
  end do
  end do
  end do

  ! assure that aux arrays do not contain NaNs
  y = 0
  z = 0

  !---------------------------------------------------------------------------
  ! evaluation

  !$acc data async present(cfg, D_inv, f, u) copyin(S, WS_t)
  !$acc parallel
  !$acc loop gang worker private(y,z)

  !$omp do private(l)
  subdomains: do l = 1, nd

    ! element-boundary configurations
    c1 = cfg(1,l)
    c2 = cfg(2,l)
    c3 = cfg(3,l)

    ! y = Sᵀ x I x I uᵉ
    call PROC(QtxIxI__,_NP_)(S(:,:,c3), alpha, beta, f(:,:,:,l), y)

    ! z = I x Sᵀ x I y
    call PROC(IxQtxI__,_NP_)(S(:,:,c2), alpha, beta, y, z)

    ! y = I x I x Sᵀ z
    call PROC(IxIxQt__,_NP_)(S(:,:,c1), alpha, beta, z, y)

    ! y = D⁻¹ y
    y = D_inv(:,:,:,l) * y

    ! z = I x I x WS y
    call PROC(IxIxQt__,_NP_)(WS_t(:,:,c1), alpha, beta, y, z)

    ! y  = I x WS x I z
    call PROC(IxQtxI__,_NP_)(WS_t(:,:,c2), alpha, beta, z, y)

    ! uᵉ = WS x I x I y
    u(:,:,:,l) = 0  ! avoid trouble with NaNs
    call PROC(QtxIxI__,_NP_)(WS_t(:,:,c3), alpha, beta, y, u(:,:,:,l))

  end do subdomains

  !$acc end parallel
  !$acc end data

end subroutine PROC(TPO_Schwarz_CI_Hand__,_NP_)
