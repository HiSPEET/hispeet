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

!> summary:  FAS-MG solver for elliptic problems
!> author:   Joerg Stiller
!> date:     2026/08/01
!===============================================================================

submodule(ML__DG__Elliptic_Solver__3D) MP_FAS_MG_Solver
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> FAS-MG solver for problems with constant diffusivity

  module subroutine FAS_MG_Solver_C(this, bc, lambda, nu, bv, f, u, ni, r_2)

    class(ML_DG_EllipticSolver_3D), intent(in) :: this

    character, intent(in) :: bc(:)  !< boundary conditions
    real(RNP), intent(in) :: lambda !< Helmholtz parameter
    real(RNP), intent(in) :: nu     !< diffusivity

    class(ML_BoundaryVariable_3D), intent(in)    :: bv !< boundary values
    class(ML_MeshVariable_3D),     intent(in)    :: f  !< RHS
    class(ML_MeshVariable_3D),     intent(inout) :: u  !< approx/final solution

    integer,   optional, intent(out) :: ni    !< num executed cycles
    real(RNP), optional, intent(out) :: r_2   !< Euclidean residual norm

    call FAS_MG_Solver_X(this, bc, lambda, nu, null(), bv, f, u, ni, r_2)

  end subroutine FAS_MG_Solver_C

  !-----------------------------------------------------------------------------
  !> FAS-MG solver for problems with variable diffusivity

  module subroutine FAS_MG_Solver_V(this, bc, lambda, nu, bv, f, u, ni, r_2 )

    class(ML_DG_EllipticSolver_3D), intent(in) :: this

    character, intent(in) :: bc(:)  !< boundary conditions
    real(RNP), intent(in) :: lambda !< Helmholtz parameter

    class(ML_MeshVariable_3D),     intent(in)    :: nu !< diffusivity
    class(ML_BoundaryVariable_3D), intent(in)    :: bv !< boundary values
    class(ML_MeshVariable_3D),     intent(in)    :: f  !< RHS
    class(ML_MeshVariable_3D),     intent(inout) :: u  !< approx/final solution

    integer,   optional, intent(out) :: ni    !< num executed cycles
    real(RNP), optional, intent(out) :: r_2   !< Euclidean residual norm

    call FAS_MG_Solver_X(this, bc, lambda, null(), nu, bv, f, u, ni, r_2)

  end subroutine FAS_MG_Solver_V

  !-----------------------------------------------------------------------------
  !> Generic FAS-MG solver for problems with constant or variable diffusivity
  !>
  !> Either `nu_0` or `nu_v` must be given.

  module subroutine FAS_MG_Solver_X( this, bc, lambda, nu_0, nu_v, bv, f, u &
                                   , ni, r_2 )

    class(ML_DG_EllipticSolver_3D), intent(in) :: this

    character, intent(in) :: bc(:)
      !< boundary conditions
    real(RNP), intent(in) :: lambda
      !< Helmholtz parameter
    real(RNP), optional, intent(in) :: nu_0
      !< constant diffusivity
    class(ML_MeshVariable_3D), optional, intent(in) :: nu_v
      !< variable diffusivity
    class(ML_BoundaryVariable_3D), intent(in) :: bv
      !< boundary values
    class(ML_MeshVariable_3D), intent(in) :: f
      !< RHS
    class(ML_MeshVariable_3D), intent(inout) :: u
      !< approx/final solution
    integer, optional, intent(out) :: ni
      !< num executed cycles
    real(RNP), optional, intent(out) :: r_2
      !< Euclidean residual norm

    ! internal variables .......................................................

    type(ML_MeshVariable_3D), allocatable, save :: r
    real(RNP), allocatable, save :: r0_2
    integer :: l_top

    character(len=:), allocatable :: prefix
    logical :: logging

    ! start logging ............................................................

    if (log_level > 0) then
      associate(proc => this % ml_op % sem(1) % mesh % proc)
        logging = proc == 0 .or. log_level > 1
        prefix  = LoggingPrefix('FAS_MG_Solver_X', proc)
      end associate
    else
      logging = .false.
    end if

    if (logging) then
      print '(2A)', prefix, 'start'
    end if

    ! initialization ...........................................................

    l_top = size(this % ml_op % sem)

    ! initial residual
    if (this%r_red > 0) then
      !$omp master
      allocate(r)
      call r % Init(this%ml_op, nc = 1, l_top = l_top)
      !$omp end master
      call this % FAS_MG_Residual_X(bc, lambda, nu_0, nu_v, bv, f, u, r)
      r0_2 = sqrt(ML_ScalarProduct_3D(r,r))
      !$omp master
      deallocate(r)
      !$omp end master
    end if

    ! start ....................................................................

    select case(this % start_method)
    case(START_CASC)
      call this % FAS_MG_Start_X(bc, lambda, nu_0, nu_v, bv, f, u, n_cyc = 0)
    case(START_FMG:)
      call this % FAS_MG_Start_X(bc, lambda, nu_0, nu_v, bv, f, u,         &
                                 n_cyc = this%start_method + 1 - START_FMG )
    end select

    ! V cycles .................................................................

    call this % FAS_MG_Cycle_X( bc, lambda, nu_0, nu_v, bv, f, u &
                              , r0_2 = r0_2, ni = ni, r_2 = r_2  )

    !$omp master
    if (allocated(r0_2)) deallocate(r0_2)
    !$omp end master

    ! exit logging .............................................................

    if (logging) then
      print '(2A)', prefix, 'exit'
    end if

  end subroutine FAS_MG_Solver_X

  !=============================================================================

end submodule MP_FAS_MG_Solver
