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

!> summary:  Smoother for elliptic problems with constant diffusivity
!> author:   Joerg Stiller
!> date:     2020/11/20
!===============================================================================

submodule(ML__DG__Elliptic_Solver__3D) MP_Smoother
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Smoother with constant diffusivity

  module subroutine Smoother_C(this, l, bc, lambda, nu, u, f, bv, n_s)
    class(ML_DG_EllipticSolver_3D), intent(in) :: this
    integer,   intent(in) :: l
    character, intent(in) :: bc(:)
    real(RNP), intent(in) :: lambda
    real(RNP), intent(in) :: nu
    real(RNP), contiguous, intent(inout) :: u(:,:,:,:)
    real(RNP), contiguous, intent(in) :: f(:,:,:,:)
    class(BoundaryVariable_3D), optional, intent(in) :: bv(:)
    integer, intent(in) :: n_s

    call Smoother_X(this, l, bc, lambda, nu, null(), u, f, bv, n_s)

  end subroutine Smoother_C

  !-----------------------------------------------------------------------------
  !> Smoother with variable diffusivity

  module subroutine Smoother_V(this, l, bc, lambda, nu, u, f, bv, n_s)
    class(ML_DG_EllipticSolver_3D), intent(in) :: this
    integer,   intent(in) :: l
    character, intent(in) :: bc(:)
    real(RNP), intent(in) :: lambda
    real(RNP), contiguous, intent(in) :: nu(:,:,:,:)
    real(RNP), contiguous, intent(inout) :: u(:,:,:,:)
    real(RNP), contiguous, intent(in) :: f(:,:,:,:)
    class(BoundaryVariable_3D), optional, intent(in) :: bv(:)
    integer, intent(in) :: n_s

    call Smoother_X(this, l, bc, lambda, null(), nu, u, f, bv, n_s)

  end subroutine Smoother_V

  !-----------------------------------------------------------------------------
  !> Generic smoother with constant or variable diffusivity

  subroutine Smoother_X(this, l, bc, lambda, nu_c, nu_v, u, f, bv, n_s)
    class(ML_DG_EllipticSolver_3D), intent(in) :: this
    integer,   intent(in) :: l
    character, intent(in) :: bc(:)
    real(RNP), intent(in) :: lambda
    real(RNP), optional, intent(in) :: nu_c
    real(RNP), contiguous, optional, intent(in) :: nu_v(:,:,:,:)
    real(RNP), contiguous, intent(inout) :: u(:,:,:,:)
    real(RNP), contiguous, intent(in) :: f(:,:,:,:)
    class(BoundaryVariable_3D), optional, intent(in) :: bv(:)
    integer, intent(in) :: n_s

    if (n_s < 1) return

    select case(this % smooth_method)
    case(SOLVER_CG)
      ! Flexible CG
      call this % elliptic_op(l) % CG_Method_X &
             (bc, lambda, nu_c, nu_v, u, f, bv, n_s)
    case(SOLVER_WS)
      ! Weighted Additive Schwarz
      call this % elliptic_op(l) % Schwarz_Method_X &
             (bc, lambda, nu_c, nu_v, u, f, bv, n_s)
    case(SOLVER_SPCG)
      ! Schwarz-preconditioned flexible CG
      call this % elliptic_op(l) % SchwarzPCG_Method_X &
             (bc, lambda, nu_c, nu_v, u, f, bv, n_s)
    end select

  end subroutine Smoother_X

  !=============================================================================

end submodule MP_Smoother
