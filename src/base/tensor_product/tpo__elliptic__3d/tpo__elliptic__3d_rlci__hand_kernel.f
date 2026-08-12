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
!> Parametrized 3d diffusion kernel using hand-crafted suboperators (RLCI)

subroutine PROC(TPO_Elliptic_RLCI_Hand__,_NP_) &
    (ne, Ms, Ls, dx, lambda, nu, u, v, Ds, ub, qb)

  integer,   intent(in)  :: ne                   !< num elements
  real(RWP), intent(in)  :: Ms(_NP_)             !< standard mass matrix
  real(RWP), intent(in)  :: Ls(_NP_,_NP_)        !< standard stiffness matrix
  real(RWP), intent(in)  :: dx(3)                !< element extensions
  real(RWP), intent(in)  :: lambda               !< Helmholtz parameter λ
  real(RWP), intent(in)  :: nu                   !< diffusivity nu
  real(RWP), intent(in)  :: u(_NP_,_NP_,_NP_,ne) !< operand
  real(RWP), intent(out) :: v(_NP_,_NP_,_NP_,ne) !< result

  real(RWP), optional, intent(in)    :: Ds(_NP_,_NP_)      !< 1D std diff matrix
  real(RWP), optional, intent(inout) :: ub(_NP_,_NP_,6,ne) !< elem bound values
  real(RWP), optional, intent(inout) :: qb(_NP_,_NP_,6,ne) !< elem bound fluxes

  !-----------------------------------------------------------------------------
  ! local variables

  real(RWP), parameter :: beta = 1
  real(RWP) :: M(_NP_,_NP_,_NP_), M_u(_NP_,_NP_,_NP_), Lm(_NP_,_NP_)
  real(RWP) :: Bs1(_NP_), Bs2(_NP_)
  real(RWP) :: c(3), tmp, tmp1, tmp2

  integer, parameter :: NP_E = _NP_**3
  integer :: e, i, j, k, p
  logical :: get_traces

  !-----------------------------------------------------------------------------
  ! initialization

  ! element mass matrix
  tmp = product(dx) / 8
  do k = 1, _NP_
  do j = 1, _NP_
  do i = 1, _NP_
    M(i,j,k) = tmp * Ms(k) * Ms(j) * Ms(i)
  end do
  end do
  end do

  ! mass-weighted stiffness matrix: Lm = Ms^-1 Ls = (Ls Ms^-1)^T
  do j = 1, _NP_
  do i = 1, _NP_
    Lm(i,j) = Ls(i,j) / Ms(i)
  end do
  end do

  ! coefficients
  c = 4 * nu / dx**2

  ! handling of boundary traces
  get_traces = present(Ds) .and. present(ub) .and. present(qb)

  if (get_traces) then
    Bs1 = -nu * Ds( 1  ,:)  ! n⋅νDˢ @ "left" boundary
    Bs2 =  nu * Ds(_NP_,:)  ! n⋅νDˢ @ "right" boundary
  end if

  !---------------------------------------------------------------------------
  ! evaluation

  !$acc data present(u,v) copyin(c,M,Lm,Bs1,Bs2)
  !$acc parallel
  !$acc loop gang worker private(M_u)

  !$omp do
  do e = 1, ne

    call SetOperands(lambda, M, u(:,:,:,e), M_u, v(:,:,:,e))
    call PROC(IxIxQt__,_NP_)(Lm, c(1), beta, M_u, v(:,:,:,e))
    call PROC(IxQtxI__,_NP_)(Lm, c(2), beta, M_u, v(:,:,:,e))
    call PROC(QtxIxI__,_NP_)(Lm, c(3), beta, M_u, v(:,:,:,e))

    ! extract/compute traces ...................................................

    if (get_traces) then

      ! ub = u, n⋅ν∇u  @ Γ₁ + Γ₂
      do k = 1, _NP_
      do j = 1, _NP_
        tmp1 = 0
        tmp2 = 0
        do p = 1, _NP_
          tmp1 = tmp1 + Bs1(p) * u(p,j,k,e)
          tmp2 = tmp2 + Bs2(p) * u(p,j,k,e)
        end do
        ub(j,k,1,e) = u( 1  ,j,k,e)
        ub(j,k,2,e) = u(_NP_,j,k,e)
        qb(j,k,1,e) = 2/dx(1) * tmp1
        qb(j,k,2,e) = 2/dx(1) * tmp2
      end do
      end do

      ! ub = u, n⋅ν∇u  @ Γ₃ + Γ₄
      do k = 1, _NP_
      do i = 1, _NP_
        tmp1 = 0
        tmp2 = 0
        do p = 1, _NP_
          tmp1 = tmp1 + Bs1(p) * u(i,p,k,e)
          tmp2 = tmp2 + Bs2(p) * u(i,p,k,e)
        end do
        ub(i,k,3,e) = u(i, 1  ,k,e)
        ub(i,k,4,e) = u(i,_NP_,k,e)
        qb(i,k,3,e) = 2/dx(2) * tmp1
        qb(i,k,4,e) = 2/dx(2) * tmp2
      end do
      end do

      ! ub = u, n⋅ν∇u  @ Γ₅ + Γ₆
      do j = 1, _NP_
      do i = 1, _NP_
        tmp1 = 0
        tmp2 = 0
        do p = 1, _NP_
          tmp1 = tmp1 + Bs1(p) * u(i,j,p,e)
          tmp2 = tmp2 + Bs2(p) * u(i,j,p,e)
        end do
        ub(i,j,5,e) = u(i,j, 1  ,e)
        ub(i,j,6,e) = u(i,j,_NP_,e)
        qb(i,j,5,e) = 2/dx(3) * tmp1
        qb(i,j,6,e) = 2/dx(3) * tmp2
      end do
      end do

    end if

  end do

  !$acc end parallel
  !$acc end data

contains

  !-----------------------------------------------------------------------------
  !> Projection to element mass matrix and computation of the Helmholtz term

  subroutine SetOperands(lambda, M, u, M_u, v)
    !$acc routine vector
    real(RWP), intent(in)  :: lambda
    real(RWP), intent(in)  :: M   (NP_E)
    real(RWP), intent(in)  :: u   (NP_E)
    real(RWP), intent(out) :: M_u (NP_E)
    real(RWP), intent(out) :: v   (NP_E)

    integer :: l

    !$acc loop vector
    !DIR$ SIMD
    do l = 1, NP_E
      M_u(l) = M(l) * u(l)
      v(l) = lambda * M_u(l)
    end do

  end subroutine SetOperands

  !-----------------------------------------------------------------------------

end subroutine PROC(TPO_Elliptic_RLCI_Hand__,_NP_)
