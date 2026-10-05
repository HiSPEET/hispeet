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

!> summary:  Incompressible Navier-Stokes DG-SEM diffusion residual (DV)
!> author:   Joerg Stiller
!> date:     2024/09/30
!===============================================================================

submodule(INS__Operator__3D) MP_GetDiffusionResidual_V
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Diffusion residual with constant viscosity
  !>
  !> Computes the DG-SEM residual of  of the viscous diffusion term including
  !> the implicit part of the discretized time derivative, i.e.,
  !>
  !>     r = F_d(v, vb, sb) - Mv/τ + Mf
  !>
  !> where `F_d` is the weak form of the diffusion term for the given velocity
  !> `v` and boundary values `vb`, `sb` obtained from the boundary variable
  !> `bv`, `M` is the diagonal mass matrix and `f` the nodal coefficients of
  !> the sources, which comprise the remaining coefficients of the momentum
  !> equation.

  module subroutine GetDiffusionResidual_V(this, tau, mu, nu, bv, f, v, r, form)

    class(INS_Operator_3D), intent(in) :: this
    !< incompressible Navier-Stokes operator

    real(RNP), intent(in) :: tau
    !< τ, effective time step width

    real(RNP), contiguous, intent(in) :: mu(:,:,:,:)
    !< kinematic bulk viscosity μ (np,np,np,ne)

    real(RNP), contiguous, intent(in) :: nu(:,:,:,:)
    !< kinematic shear viscosity ν (np,np,np,ne)

    class(BoundaryVariable_3D), intent(in) :: bv(:)
    !< boundary values
    !!   - Γᴰ :  [ v₁, v₂, v₃, - ]
    !!   - Γᴼ :  [ - , - , ∆p, p ]

    real(RNP), contiguous, intent(in) :: f(:,:,:,:,:)
    !< sources, f(np,np,np,ne,3)

    real(RNP), contiguous, intent(in) :: v(:,:,:,:,:)
    !< velocity, v(np,np,np,ne,3)

    real(RNP), contiguous, intent(out) :: r(:,:,:,:,:)
    !< residual, r(np,np,np,ne,3)

    integer, optional, intent(in) :: form
    !< form of `∇⋅τ`: 0/1/2 ↔︎ default/diffusion/rotational [0]

    ! local variables ..........................................................

    real(RNP), allocatable, save :: mm(:,:,:,:)   ! diagonal mass matrix M
    real(RNP), allocatable, save :: vp(:,:,:,:,:) ! velocity traces v⁺
    real(RNP), allocatable, save :: sp(:,:,:,:,:) ! viscous flux traces s⁺

    real(RNP) :: lambda
    integer   :: na, ne, np
    integer   :: d, e

    associate(mesh => this % sem_u % mesh)

      ! initialization .........................................................

      np = size(v,1)
      na = mesh % n_elem_active
      ne = mesh % n_elem

      !$omp master
      allocate( mm (np, np, np, ne) )
      allocate( vp (np, np,  6, ne, 3), source = ZERO )
      allocate( sp (np, np,  6, ne, 3), source = ZERO )
      !$omp end master
      !$omp barrier

      call this % sem_u % Get_DG_DiagonalMassMatrix(mm)

      lambda = 1 / tau

      ! compute residual .......................................................

      call this % GetDiffusionTerm_V(mu, nu, bv, r, v, vp, sp, form=form)

      !$omp do collapse(2)
      do e = 1, na
        do d = 1, 3
          r(:,:,:,e,d) = r(:,:,:,e,d) &
                       + mm(:,:,:,e) * (f(:,:,:,e,d) - lambda * v(:,:,:,e,d))
        end do
      end do
      !$omp end do nowait

      !$omp do collapse(2)
      do e = na+1, ne
        do d = 1, 3
          r(:,:,:,e,d) = 0
        end do
      end do

      ! cleanup ................................................................

      !$omp master
      deallocate(mm, vp, sp)
      !$omp end master

    end associate

  end subroutine GetDiffusionResidual_V

  !=============================================================================

end submodule MP_GetDiffusionResidual_V
