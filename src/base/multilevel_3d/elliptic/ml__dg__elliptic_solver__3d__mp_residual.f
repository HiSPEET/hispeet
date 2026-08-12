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

!> summary:  Residual of elliptic problems
!> author:   Joerg Stiller
!> date:     2020/11/20
!===============================================================================
submodule(ML__DG__Elliptic_Solver__3D) MP_Residual
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Residual with constant diffusivity
  !>
  !> Homogeneous boundary conditions are assumed if `bv` is absent.
  !> This case is used for implementing the correction scheme

  module subroutine Residual_C(this, l, bc, lambda, nu, f, bv, u, r)
    class(ML_DG_EllipticSolver_3D), intent(in) :: this
    integer,   intent(in) :: l
    character, intent(in) :: bc(:)
    real(RNP), intent(in) :: lambda
    real(RNP), intent(in) :: nu
    real(RNP), contiguous, intent(in) :: f(:,:,:,:)
    class(BoundaryVariable_3D), optional, intent(in) :: bv(:)
    real(RNP), contiguous, intent(in)  :: u(:,:,:,:)
    real(RNP), contiguous, intent(out) :: r(:,:,:,:)

    call  this % elliptic_op(l) % Residual(bc, lambda, nu, f, bv, u, r)

  end subroutine Residual_C

  !-----------------------------------------------------------------------------
  !> Residual with variable diffusivity
  !>
  !> Homogeneous boundary conditions are assumed if `bv` is absent.
  !> This case is used for implementing the correction scheme

  module subroutine Residual_V(this, l, bc, lambda, nu, f, bv, u, r)
    class(ML_DG_EllipticSolver_3D), intent(in) :: this
    integer,   intent(in) :: l
    character, intent(in) :: bc(:)
    real(RNP), intent(in) :: lambda
    real(RNP), contiguous, intent(in) :: nu(:,:,:,:)
    real(RNP), contiguous, intent(in) :: f(:,:,:,:)
    class(BoundaryVariable_3D), optional, intent(in) :: bv(:)
    real(RNP), contiguous, intent(in)  :: u(:,:,:,:)
    real(RNP), contiguous, intent(out) :: r(:,:,:,:)

    call  this % elliptic_op(l) % Residual(bc, lambda, nu, f, bv, u, r)

  end subroutine Residual_V

  !=============================================================================

end submodule MP_Residual
