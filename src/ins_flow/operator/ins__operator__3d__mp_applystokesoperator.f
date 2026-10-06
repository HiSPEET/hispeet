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

!> summary:  Stokes operator for incompressible flow
!> author:   Joerg Stiller
!> date:     2025/05/17
!===============================================================================

submodule (INS__Operator__3D) MP_ApplyStokesOperator
  use TPO__Div__3D,  only: TPO_Div
  use TPO__Grad__3D, only: TPO_Grad
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Stokes operator for incompressible flow
  !>
  !> Skipping `bv` yields the homogeneous operator

  module subroutine ApplyStokesOperator(this, tau, mu, nu, bv, u, r)

    ! arguments ................................................................

    class(INS_Operator_3D), intent(in) :: this
    !< incompressible Navier-Stokes operator
    real(RNP), intent(in) :: tau
    !< τ, effective time step width
    real(RNP), contiguous, optional, intent(in) :: mu(:,:,:,:)
    !< μ, kinematic bulk viscosity
    real(RNP), contiguous, optional, intent(in) :: nu(:,:,:,:)
    !< ν, kinematic shear viscosity
    class(BoundaryVariable_3D), optional, intent(in) :: bv(:)
    !< boundary values
    real(RNP), contiguous, intent(in) :: u(:,:,:,:,:)
    !< u = [v, p], velocity and pressure
    real(RNP), contiguous, intent(out) :: r(:,:,:,:,:)
    !< result

    ! internal variables .......................................................

    real(RNP), allocatable, save :: mm(:,:,:,:)   ! diagonal mass matrix
    real(RNP), allocatable, save :: up(:,:,:,:,:) ! outer traces
    real(RNP), allocatable, save :: w (:,:,:,:,:) ! work

    real(RNP) :: lambda
    integer :: na, ne, np
    integer :: c, e

    ! initialization ...........................................................

    lambda = 1/tau

    np = size(u,1)
    ne = this % mesh % n_elem
    na = this % mesh % n_elem_active

    !$omp master
    allocate(mm (np, np, np, ne))
    allocate(up (np, np,  6, ne, 6), source = ZERO)
    allocate(w  (np, np, np, ne, 4), source = ZERO)
    !$omp end master
    !$omp barrier

    call this % sem_u % Get_DG_DiagonalMassMatrix(mm)

    !...........................................................................

    associate( r_m    => r (:,:,:,:,1:3) &
             , r_c    => r (:,:,:,:, 4 ) &
             , v      => u (:,:,:,:,1:3) &
             , p      => u (:,:,:,:, 4 ) &
             , vp     => up(:,:,:,:,1:3) &
             , sp     => up(:,:,:,:,4:6) &
             , pp     => up(:,:,:,:, 4 ) &
             , grad_p => w (:,:,:,:,1:3) )

      ! contributions ..........................................................

      ! diffusion
      call this % GetDiffusionTerm(mu, nu, bv, v, vp, sp, r_m)

      ! pressure gradient and velocity divergence
      call GetOuterTraces_3D(this % mesh, p, pp)
      call TPO_Grad(this % eop_u, this % sem_u, p, pp, grad_p)
      call TPO_Div( this % eop_u, this % sem_u, v, vp, r_c)

      ! complete result ........................................................

      !$omp do
      do e = 1, na
        do c = 1, 3
          r_m(:,:,:,e,c)                                                 &
             = mm(:,:,:,e) * (lambda * v(:,:,:,e,c) + grad_p(:,:,:,e,c)) &
             - r_m(:,:,:,e,c)
        end do
        r_c(:,:,:,e) = mm(:,:,:,e) * r_c(:,:,:,e)
      end do
      !$omp end do nowait

      !$omp do collapse(2)
      do c = 1, 4
      do e = na+1, ne
        r(:,:,:,e,c) = 0
      end do
      end do

    end associate

    !$omp master
    deallocate(mm, up, w)
    !$omp end master

  end subroutine ApplyStokesOperator

  !=============================================================================

end submodule MP_ApplyStokesOperator

