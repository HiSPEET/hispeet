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

!> summary:  Homogeneous incompressible Navier-Stokes DG-SEM diffusion operator
!> author:   Joerg Stiller
!> date:     2022/09/26
!===============================================================================

submodule(INS__Operator__3D) MP_ApplyDiffusionOperator_C
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Homogeneous diffusion operator with constant viscosity
  !>
  !> Computes the DG-SEM viscous diffusion operator including the implicit part
  !> of the discretized time derivative, i.e.,
  !>
  !>     r = Mv/τ - F_d(v, bv)
  !>
  !> where `F_d` is the weak form of the diffusion term for the given velocity
  !> `v` using the diagonal mass matrix `M`.
  !> Homogeneous boundary conditions are applied if `bv` is absent.

  module subroutine ApplyDiffusionOperator_C(this, tau, bv, v, r, form)

    class(INS_Operator_3D), intent(in) :: this
    !< incompressible Navier-Stokes operator

    real(RNP), intent(in) :: tau
    !< τ, effective time step width

    class(BoundaryVariable_3D), optional, intent(in) :: bv(:)
    !< boundary values
    !!   - Γᴰ :  [ v₁, v₂, v₃, - ]
    !!   - Γᴼ :  [ - , - , ∆p, p ]

    real(RNP), contiguous, intent(in) :: v(:,:,:,:,:)
    !< velocity, v(np,np,np,ne,3)

    real(RNP), contiguous, intent(out) :: r(:,:,:,:,:)
    !< result, r(np,np,np,ne,3)

    integer, optional, intent(in) :: form
    !< form of `∇⋅τ`: 0/1/2 ↔︎ default/diffusion/rotational [0]

    ! local variables ..........................................................

    real(RNP), allocatable, save :: mm(:,:,:,:)   ! diagonal mass matrix
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

      ! computation ............................................................

      call this % GetDiffusionTerm_C(bv, v, vp, sp, r, form=form)

      !$omp do collapse(2)
      do e = 1, na
        do d = 1, 3
          r(:,:,:,e,d) = lambda * mm(:,:,:,e) * v(:,:,:,e,d) - r(:,:,:,e,d)
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

  end subroutine ApplyDiffusionOperator_C

  !=============================================================================

end submodule MP_ApplyDiffusionOperator_C
