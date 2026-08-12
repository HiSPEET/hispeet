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

!> summary:  Coarse grid solver for elliptic problems
!> author:   Joerg Stiller
!> date:     2026/08/01
!===============================================================================

submodule(ML__DG__Elliptic_Solver__3D) MP_CoarseSolver
  use Boundary_Variable__3D
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Coarse grid solver with constant diffusivity

  module subroutine CoarseSolver_C(this, bc, lambda, nu, u, f, bv)
    class(ML_DG_EllipticSolver_3D), intent(in) :: this
    character, intent(in) :: bc(:)
    real(RNP), intent(in) :: lambda
    real(RNP), intent(in) :: nu
    real(RNP), contiguous, intent(inout) :: u(:,:,:,:)
    real(RNP), contiguous, intent(in) :: f(:,:,:,:)
    class(BoundaryVariable_3D), optional, intent(in) :: bv(:)

    call CoarseSolver_X(this, bc, lambda, nu_c = nu, u = u, f = f, bv = bv)

  end subroutine CoarseSolver_C

  !-----------------------------------------------------------------------------
  !> Coarse grid solver with variable diffusivity

  module subroutine CoarseSolver_V(this, bc, lambda, nu, u, f, bv)
    class(ML_DG_EllipticSolver_3D), intent(in) :: this
    character, intent(in) :: bc(:)
    real(RNP), intent(in) :: lambda
    real(RNP), contiguous, intent(in) :: nu(:,:,:,:)
    real(RNP), contiguous, intent(inout) :: u(:,:,:,:)
    real(RNP), contiguous, intent(in) :: f(:,:,:,:)
    class(BoundaryVariable_3D), optional, intent(in) :: bv(:)

    call CoarseSolver_X(this, bc, lambda, nu_v = nu, u = u, f = f, bv = bv)

  end subroutine CoarseSolver_V

  !-----------------------------------------------------------------------------
  !> Generic coarse grid solver with constant or variable diffusivity

  subroutine CoarseSolver_X(this, bc, lambda, nu_c, nu_v, u, f, bv)
    class(ML_DG_EllipticSolver_3D), intent(in) :: this
    character, intent(in) :: bc(:)
    real(RNP), intent(in) :: lambda
    real(RNP), optional, intent(in) :: nu_c
    real(RNP), contiguous, optional, intent(in) :: nu_v(:,:,:,:)
    real(RNP), contiguous, intent(inout) :: u(:,:,:,:)
    real(RNP), contiguous, intent(in) :: f(:,:,:,:)
    class(BoundaryVariable_3D), optional, intent(in) :: bv(:)

    select case(this % coarse_solver)
    case(SOLVER_CG)
      ! Flexible CG
      call this % elliptic_op(1) %                                    &
                      CG_Method_X( bc, lambda, nu_c, nu_v , u, f, bv  &
                                 , this%i_crs, this%r_crs, this%r_max )
    case(SOLVER_WS)
      ! Weighted Additive Schwarz
      call this % elliptic_op(1) %                                         &
                      Schwarz_Method_X( bc, lambda, nu_c, nu_v, u, f, bv   &
                                      , this%i_crs, this%r_crs, this%r_max )
    case(SOLVER_SPCG)
      ! Schwarz-preconditioned flexible CG
      call this % elliptic_op(1) %                                            &
                      SchwarzPCG_Method_X( bc, lambda, nu_c, nu_v, u, f, bv   &
                                         , this%i_crs, this%r_crs, this%r_max )
    end select

  end subroutine CoarseSolver_X

  !=============================================================================

end submodule MP_CoarseSolver
