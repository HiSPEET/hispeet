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

!> summary:  3D FAS multigrid DG solver for elliptic problems
!> author:   Joerg Stiller
!> date:     2024/08/20
!===============================================================================

module ML__DG__Elliptic_Solver__3D
  use Kind_Parameters
  use Constants
  use XMPI
  use Logging_Levels
  use Array_Assignments
  use Array_Reductions
  use Boundary_Variable__3D
  use DG__Element_Operators__1D
  use DG__Elliptic_Operator__3D
  use DG__Schwarz_Operator__3D
  use Child_To_Parent_Projection__3D
  use Child_To_Parent_Restriction__3D
  use Parent_To_Child_Interpolation__3D
  use ML__Boundary_Variable__3D
  use ML__Mesh_Operators__3D
  use ML__Mesh_Variable__3D
  use ML__Array_Assignments__3D
  use ML__Array_Reductions__3D

  implicit none
  private

  public :: ML_DG_EllipticSolver_3D
  public :: ML_DG_EllipticOptions_3D

  !=============================================================================
  ! Parameters

  integer, parameter :: START_CASC  = 1 !< start with cascade
  integer, parameter :: START_FMG   = 2 !< start with full multigrid method

  integer, parameter :: SOLVER_CG   = 1 !< flexible conjugate gradient method
  integer, parameter :: SOLVER_WS   = 2 !< weighted Schwarz method
  integer, parameter :: SOLVER_SPCG = 3 !< Schwarz-preconditioned flexible CG

  !=============================================================================
  ! Types

  !-----------------------------------------------------------------------------
  !> Type providing 3D FAS multigrid DG solvers for elliptic problems

  type ML_DG_EllipticSolver_3D

    class(ML_MeshOperators_3D), pointer :: ml_op => null()
      !< multilevel mesh operators

    type(DG_EllipticOperator_3D), allocatable :: elliptic_op(:)
      !< elliptic operators for each level

    integer   :: start_method   !< starting method
    integer   :: smooth_method  !< smoothing method
    integer   :: coarse_solver  !< coarse grid solver
    character :: fc_projection  !< fine-to-coarse projection method
    character :: fc_restriction !< fine-to-coarse restriction method

    integer   :: i_crs   !< max number of coarse solver iterations
    integer   :: i_max   !< max number of multigrid iterations (cycles)
    integer   :: ns_0    !< num smoothing steps in starting cascade
    integer   :: ns_1    !< num pre-smoothing steps
    integer   :: ns_2    !< num post-smoothing steps
    integer   :: ns_c    !< num continuation smoothing steps
    integer   :: ns_f    !< num final smoothing steps on top level
    integer   :: vs_mode !< variable smoothing mode, 0/1/2: none/p/hp
    integer   :: vs_lmin !< variable smoothing minimum level
    integer   :: vs_lmax !< variable smoothing maximum level
    real(RNP) :: r_crs   !< min coarse residual reduction,  if > 0
    real(RNP) :: r_red   !< min residual reduction,  if > 0
    real(RNP) :: r_max   !< max admissible residual, if > 0

  contains

    procedure, public  :: Init_ML_DG_EllipticSolver_3D

    procedure, public  :: NumSmoothingSteps

    generic,   public  :: CS_MG_Solver => CS_MG_Solver_C, CS_MG_Solver_V
    procedure, private :: CS_MG_Solver_C, CS_MG_Solver_V, CS_MG_Solver_X

    generic,   public  :: CS_MGCG_Solver => CS_MGCG_Solver_C, CS_MGCG_Solver_V
    procedure, private :: CS_MGCG_Solver_C, CS_MGCG_Solver_V

    generic,   public  :: FAS_MG_Residual => FAS_MG_Residual_C, FAS_MG_Residual_V
    procedure, private :: FAS_MG_Residual_C, FAS_MG_Residual_V, FAS_MG_Residual_X

    procedure, private :: FAS_MG_Cycle_X
    procedure, private :: FAS_MG_Start_X

    generic,   public  :: FAS_MG_Solver => FAS_MG_Solver_C, FAS_MG_Solver_V
    procedure, private :: FAS_MG_Solver_C, FAS_MG_Solver_V, FAS_MG_Solver_X

    generic,   private :: Residual => Residual_C, Residual_V
    procedure, private :: Residual_C, Residual_V

    generic,   private :: CoarseSolver => CoarseSolver_C, CoarseSolver_V
    procedure, private :: CoarseSolver_C, CoarseSolver_V

    generic,   private :: Smoother => Smoother_C, Smoother_V
    procedure, private :: Smoother_C, Smoother_V

    generic,   private :: Monitoring => Monitoring_R, Monitoring_C, Monitoring_V
    procedure, private :: Monitoring_R, Monitoring_C, Monitoring_V

  end type ML_DG_EllipticSolver_3D

  ! constructor interface
  interface ML_DG_EllipticSolver_3D
    procedure New_ML_DG_EllipticSolver_3D
  end interface

  !-----------------------------------------------------------------------------
  !> 3D FAS multigrid solver options

  type ML_DG_EllipticOptions_3D

    real(RNP) :: penalty = -1                !< IP-DG penalty, keep for default
    type(DG_SchwarzOptions_3D) :: schwarz    !< Schwarz operator options

    integer   :: start_method  = START_FMG   !< starting method
    integer   :: smooth_method = SOLVER_WS   !< smoothing method
    integer   :: coarse_solver = SOLVER_SPCG !< coarse grid solver

    character :: fc_projection  = 'I' !< projection method {'I','P'}
    character :: fc_restriction = 'C' !< restriction method {'C','P'}
    character :: interior_bc    = ' ' !< coupling with frozen elements {' ','D'}

    integer   :: i_crs   =  1 !< max num coarse solver iterations
    integer   :: i_max   =  1 !< max num multigrid iterations (cycles)
    integer   :: ns_0    =  1 !< num smoothing steps in starting cascade
    integer   :: ns_1    =  1 !< num pre-smoothing steps
    integer   :: ns_2    =  1 !< num post-smoothing steps
    integer   :: ns_c    = -1 !< num continuation smoothing steps         [auto]
    integer   :: ns_f    = -1 !< num final smoothing steps on top level   [auto]
    integer   :: vs_mode =  0 !< variable smoothing mode, 0/1/2: none/p/hp
    integer   :: vs_lmin =  1 !< variable smoothing minimum level
    integer   :: vs_lmax = -1 !< variable smoothing maximum level        [l_top]
    real(RNP) :: r_crs   = -1 !< min coarse residual reduction, overrides `r_red`
    real(RNP) :: r_red   = -1 !< min residual reduction,  if > 0
    real(RNP) :: r_max   = -1 !< max admissible residual, if > 0

  contains

    procedure :: Bcast => Bcast_ML_DG_EllipticOptions_3D

  end type ML_DG_EllipticOptions_3D

  !=============================================================================
  ! ML_DG_EllipticSolver_3D: separate multi-level procedures

  interface

    !---------------------------------------------------------------------------
    !> CS-MG solver for problems with global refinement and constant diffusivity

    module subroutine CS_MG_Solver_C( this, bc, lambda, nu, bv, f, u &
                                    , i_max, l_top, ni, r_2 )
      class(ML_DG_EllipticSolver_3D), intent(in)    :: this
      character,                      intent(in)    :: bc(:)
      real(RNP),                      intent(in)    :: lambda
      real(RNP),                      intent(in)    :: nu
      class(ML_BoundaryVariable_3D),  intent(in)    :: bv
      class(ML_MeshVariable_3D),      intent(inout) :: f
      class(ML_MeshVariable_3D),      intent(inout) :: u
      integer,              optional, intent(in)    :: i_max
      integer,              optional, intent(in)    :: l_top
      integer,              optional, intent(out)   :: ni
      real(RNP),            optional, intent(out)   :: r_2
    end subroutine CS_MG_Solver_C

    !---------------------------------------------------------------------------
    !> CS-MG solver for problems with global refinement and variable diffusivity

    module subroutine CS_MG_Solver_V( this, bc, lambda, nu, bv, f, u &
                                    , i_max, l_top, ni, r_2          )
      class(ML_DG_EllipticSolver_3D), intent(in)    :: this
      character,                      intent(in)    :: bc(:)
      real(RNP),                      intent(in)    :: lambda
      class(ML_MeshVariable_3D),      intent(in)    :: nu
      class(ML_BoundaryVariable_3D),  intent(in)    :: bv
      class(ML_MeshVariable_3D),      intent(inout) :: f
      class(ML_MeshVariable_3D),      intent(inout) :: u
      integer,              optional, intent(in)    :: i_max
      integer,              optional, intent(in)    :: l_top
      integer,              optional, intent(out)   :: ni
      real(RNP),            optional, intent(out)   :: r_2
    end subroutine CS_MG_Solver_V

    !---------------------------------------------------------------------------
    !> Generic CS-MG solver for problems with constant or variable diffusivity

    module subroutine CS_MG_Solver_X( this, bc, lambda, nu_0, nu_v, bv, f, u &
                                    , i_max, l_top, ni, r_2                  )
      class(ML_DG_EllipticSolver_3D),                  intent(in)    :: this
      character,                                       intent(in)    :: bc(:)
      real(RNP),                                       intent(in)    :: lambda
      real(RNP),                             optional, intent(in)    :: nu_0
      class(ML_MeshVariable_3D),             optional, intent(in)    :: nu_v
      class(ML_BoundaryVariable_3D), target, optional, intent(in)    :: bv
      class(ML_MeshVariable_3D),                       intent(inout) :: f
      class(ML_MeshVariable_3D),                       intent(inout) :: u
      integer,                               optional, intent(in)    :: i_max
      integer,                               optional, intent(in)    :: l_top
      integer,                               optional, intent(out)   :: ni
      real(RNP),                             optional, intent(out)   :: r_2
    end subroutine CS_MG_Solver_X

    !---------------------------------------------------------------------------
    !> CS-MGCG solver for problems with global refinement and const diffusivity

    module subroutine CS_MGCG_Solver_C( this, bc, lambda, nu, bv, f, u &
                                      , i_max, l_top, ni, r_2          )
      class(ML_DG_EllipticSolver_3D), intent(in)    :: this
      character,                      intent(in)    :: bc(:)
      real(RNP),                      intent(in)    :: lambda
      real(RNP),                      intent(in)    :: nu
      class(ML_BoundaryVariable_3D),  intent(in)    :: bv
      class(ML_MeshVariable_3D),      intent(inout) :: f
      class(ML_MeshVariable_3D),      intent(inout) :: u
      integer,              optional, intent(in)    :: i_max
      integer,              optional, intent(in)    :: l_top
      integer,              optional, intent(out)   :: ni
      real(RNP),            optional, intent(out)   :: r_2
    end subroutine CS_MGCG_Solver_C

    !---------------------------------------------------------------------------
    !> CS-MG solver for problems with global refinement and variable diffusivity

    module subroutine CS_MGCG_Solver_V( this, bc, lambda, nu, bv, f, u &
                                    , i_max, l_top, ni, r_2          )
      class(ML_DG_EllipticSolver_3D), intent(in)    :: this
      character,                      intent(in)    :: bc(:)
      real(RNP),                      intent(in)    :: lambda
      class(ML_MeshVariable_3D),      intent(in)    :: nu
      class(ML_BoundaryVariable_3D),  intent(in)    :: bv
      class(ML_MeshVariable_3D),      intent(inout) :: f
      class(ML_MeshVariable_3D),      intent(inout) :: u
      integer,              optional, intent(in)    :: i_max
      integer,              optional, intent(in)    :: l_top
      integer,              optional, intent(out)   :: ni
      real(RNP),            optional, intent(out)   :: r_2
    end subroutine CS_MGCG_Solver_V

    !---------------------------------------------------------------------------
    !> FAS-MG residual with constant diffusivity

    module subroutine FAS_MG_Residual_C( this, bc, lambda, nu, bv, f &
                                       , u, r, l_top )
      class(ML_DG_EllipticSolver_3D), intent(in)    :: this
      character,                      intent(in)    :: bc(:)
      real(RNP),                      intent(in)    :: lambda
      real(RNP),                      intent(in)    :: nu
      class(ML_BoundaryVariable_3D),  intent(in)    :: bv
      class(ML_MeshVariable_3D),      intent(in)    :: f
      class(ML_MeshVariable_3D),      intent(in)    :: u
      class(ML_MeshVariable_3D),      intent(inout) :: r
      integer,              optional, intent(in)    :: l_top
    end subroutine FAS_MG_Residual_C

    !---------------------------------------------------------------------------
    !> FAS-MG residual with variable diffusivity

    module subroutine FAS_MG_Residual_V( this, bc, lambda, nu, bv, f &
                                       , u, r, l_top )
      class(ML_DG_EllipticSolver_3D), intent(in)    :: this
      character,                      intent(in)    :: bc(:)
      real(RNP),                      intent(in)    :: lambda
      class(ML_MeshVariable_3D),      intent(in)    :: nu
      class(ML_BoundaryVariable_3D),  intent(in)    :: bv
      class(ML_MeshVariable_3D),      intent(in)    :: f
      class(ML_MeshVariable_3D),      intent(in)    :: u
      class(ML_MeshVariable_3D),      intent(inout) :: r
      integer,              optional, intent(in)    :: l_top
    end subroutine FAS_MG_Residual_V

    !---------------------------------------------------------------------------
    !> Generic FAS-MG residual with constant or variable diffusivity

    module subroutine FAS_MG_Residual_X( this, bc, lambda, nu_0, nu_v, bv, f &
                                       , u, r, l_top )
      class(ML_DG_EllipticSolver_3D),                  intent(in)    :: this
      character,                                       intent(in)    :: bc(:)
      real(RNP),                                       intent(in)    :: lambda
      real(RNP),                             optional, intent(in)    :: nu_0
      class(ML_MeshVariable_3D),             optional, intent(in)    :: nu_v
      class(ML_BoundaryVariable_3D), target, optional, intent(in)    :: bv
      class(ML_MeshVariable_3D),             optional, intent(in)    :: f
      class(ML_MeshVariable_3D),                       intent(in)    :: u
      class(ML_MeshVariable_3D),                       intent(inout) :: r
      integer,                               optional, intent(in)    :: l_top
    end subroutine FAS_MG_Residual_X

    !---------------------------------------------------------------------------
    !> Generic FAS-MG V-cycle for problems with constant or variable diffusivity

    module subroutine FAS_MG_Cycle_X( this, bc, lambda, nu_0, nu_v, bv, f, u &
                                    , n_cyc, l_top, r0_2, ni, r_2 )
      class(ML_DG_EllipticSolver_3D),                  intent(in)    :: this
      character,                                       intent(in)    :: bc(:)
      real(RNP),                                       intent(in)    :: lambda
      real(RNP),                             optional, intent(in)    :: nu_0
      class(ML_MeshVariable_3D),             optional, intent(in)    :: nu_v
      class(ML_BoundaryVariable_3D), target, optional, intent(in)    :: bv
      class(ML_MeshVariable_3D),                       intent(in)    :: f
      class(ML_MeshVariable_3D),                       intent(inout) :: u
      integer,                               optional, intent(in)    :: n_cyc
      integer,                               optional, intent(in)    :: l_top
      real(RNP),                             optional, intent(in)    :: r0_2
      integer,                               optional, intent(out)   :: ni
      real(RNP),                             optional, intent(out)   :: r_2
    end subroutine FAS_MG_Cycle_X

    !---------------------------------------------------------------------------
    !> FAS-MG solver for problems with constant diffusivity

    module subroutine FAS_MG_Solver_C(this, bc, lambda, nu, bv, f, u, ni, r_2)
      class(ML_DG_EllipticSolver_3D), intent(in)    :: this
      character,                      intent(in)    :: bc(:)
      real(RNP),                      intent(in)    :: lambda
      real(RNP),                      intent(in)    :: nu
      class(ML_BoundaryVariable_3D),  intent(in)    :: bv
      class(ML_MeshVariable_3D),      intent(in)    :: f
      class(ML_MeshVariable_3D),      intent(inout) :: u
      integer,              optional, intent(out)   :: ni
      real(RNP),            optional, intent(out)   :: r_2
    end subroutine FAS_MG_Solver_C

    !---------------------------------------------------------------------------
    !> FAS-MG solver for problems with variable diffusivity

    module subroutine FAS_MG_Solver_V(this, bc, lambda, nu, bv, f, u, ni, r_2 )
      class(ML_DG_EllipticSolver_3D), intent(in)    :: this
      character,                      intent(in)    :: bc(:)
      real(RNP),                      intent(in)    :: lambda
      class(ML_MeshVariable_3D),      intent(in)    :: nu
      class(ML_BoundaryVariable_3D),  intent(in)    :: bv
      class(ML_MeshVariable_3D),      intent(in)    :: f
      class(ML_MeshVariable_3D),      intent(inout) :: u
      integer,              optional, intent(out)   :: ni
      real(RNP),            optional, intent(out)   :: r_2
    end subroutine FAS_MG_Solver_V

    !---------------------------------------------------------------------------
    !> Generic FAS-MG solver for problems with constant or variable diffusivity

    module subroutine FAS_MG_Solver_X( this, bc, lambda, nu_0, nu_v, bv, f, u &
                                     , ni, r_2 )
      class(ML_DG_EllipticSolver_3D),      intent(in)    :: this
      character,                           intent(in)    :: bc(:)
      real(RNP),                           intent(in)    :: lambda
      real(RNP),                 optional, intent(in)    :: nu_0
      class(ML_MeshVariable_3D), optional, intent(in)    :: nu_v
      class(ML_BoundaryVariable_3D),       intent(in)    :: bv
      class(ML_MeshVariable_3D),           intent(in)    :: f
      class(ML_MeshVariable_3D),           intent(inout) :: u
      integer,                   optional, intent(out)   :: ni
      real(RNP),                 optional, intent(out)   :: r_2
    end subroutine FAS_MG_Solver_X

    !---------------------------------------------------------------------------
    !> Cascade and FMG start for problems with constant or variable diffusivity

    module subroutine FAS_MG_Start_X( this, bc, lambda, nu_0, nu_v, bv, f, u &
                                    , n_cyc )
      class(ML_DG_EllipticSolver_3D),      intent(in)    :: this
      character,                           intent(in)    :: bc(:)
      real(RNP),                           intent(in)    :: lambda
      real(RNP),                 optional, intent(in)    :: nu_0
      class(ML_MeshVariable_3D), optional, intent(in)    :: nu_v
      class(ML_BoundaryVariable_3D),       intent(in)    :: bv
      class(ML_MeshVariable_3D),           intent(in)    :: f
      class(ML_MeshVariable_3D),           intent(inout) :: u
      integer,                             intent(in)    :: n_cyc
    end subroutine FAS_MG_Start_X

  end interface

  !=============================================================================
  ! ML_DG_EllipticSolver_3D:  separate single-level procedures

  interface

    !---------------------------------------------------------------------------
    !> Residual with constant diffusivity

    module subroutine Residual_C(this, l, bc, lambda, nu, f, bv, u, r)
      class(ML_DG_EllipticSolver_3D),       intent(in) :: this
      integer,                              intent(in)  :: l
      character,                            intent(in)  :: bc(:)
      real(RNP),                            intent(in)  :: lambda
      real(RNP),                            intent(in)  :: nu
      real(RNP), contiguous,                intent(in)  :: f(:,:,:,:)
      class(BoundaryVariable_3D), optional, intent(in)  :: bv(:)
      real(RNP), contiguous,                intent(in)  :: u(:,:,:,:)
      real(RNP), contiguous,                intent(out) :: r(:,:,:,:)
    end subroutine Residual_C

    !---------------------------------------------------------------------------
    !> Residual with variable diffusivity

    module subroutine Residual_V(this, l, bc, lambda, nu, f, bv, u, r)
      class(ML_DG_EllipticSolver_3D),       intent(in)  :: this
      integer,                              intent(in)  :: l
      character,                            intent(in)  :: bc(:)
      real(RNP),                            intent(in)  :: lambda
      real(RNP), contiguous,                intent(in)  :: nu(:,:,:,:)
      real(RNP), contiguous,                intent(in)  :: f(:,:,:,:)
      class(BoundaryVariable_3D), optional, intent(in)  :: bv(:)
      real(RNP), contiguous,                intent(in)  :: u(:,:,:,:)
      real(RNP), contiguous,                intent(out) :: r(:,:,:,:)
    end subroutine Residual_V

    !---------------------------------------------------------------------------
    !> Coarse grid solver with constant diffusivity

    module subroutine CoarseSolver_C(this, bc, lambda, nu, u, f, bv)
      class(ML_DG_EllipticSolver_3D),       intent(in)    :: this
      character,                            intent(in)    :: bc(:)
      real(RNP),                            intent(in)    :: lambda
      real(RNP),                            intent(in)    :: nu
      real(RNP), contiguous,                intent(inout) :: u(:,:,:,:)
      real(RNP), contiguous,                intent(in)    :: f(:,:,:,:)
      class(BoundaryVariable_3D), optional, intent(in)    :: bv(:)
    end subroutine CoarseSolver_C

    !---------------------------------------------------------------------------
    !> Coarse grid solver with variable diffusivity

    module subroutine CoarseSolver_V(this, bc, lambda, nu, u, f, bv)
      class(ML_DG_EllipticSolver_3D),       intent(in)    :: this
      character,                            intent(in)    :: bc(:)
      real(RNP),                            intent(in)    :: lambda
      real(RNP), contiguous,                intent(in)    :: nu(:,:,:,:)
      real(RNP), contiguous,                intent(inout) :: u(:,:,:,:)
      real(RNP), contiguous,                intent(in)    :: f(:,:,:,:)
      class(BoundaryVariable_3D), optional, intent(in)    :: bv(:)
    end subroutine CoarseSolver_V

    !---------------------------------------------------------------------------
    !> Smoother with constant diffusivity

    module subroutine Smoother_C(this, l, bc, lambda, nu, u, f, bv, n_s)
      class(ML_DG_EllipticSolver_3D),       intent(in)    :: this
      integer,                              intent(in)    :: l
      character,                            intent(in)    :: bc(:)
      real(RNP),                            intent(in)    :: lambda
      real(RNP),                            intent(in)    :: nu
      real(RNP), contiguous,                intent(inout) :: u(:,:,:,:)
      real(RNP), contiguous,                intent(in)    :: f(:,:,:,:)
      class(BoundaryVariable_3D), optional, intent(in)    :: bv(:)
      integer,                              intent(in)    :: n_s
    end subroutine Smoother_C

    !---------------------------------------------------------------------------
    !> Smoother with variable diffusivity

    module subroutine Smoother_V(this, l, bc, lambda, nu, u, f, bv, n_s)
      class(ML_DG_EllipticSolver_3D),       intent(in)    :: this
      integer,                              intent(in)    :: l
      character,                            intent(in)    :: bc(:)
      real(RNP),                            intent(in)    :: lambda
      real(RNP), contiguous,                intent(in)    :: nu(:,:,:,:)
      real(RNP), contiguous,                intent(inout) :: u(:,:,:,:)
      real(RNP), contiguous,                intent(in)    :: f(:,:,:,:)
      class(BoundaryVariable_3D), optional, intent(in)    :: bv(:)
      integer,                              intent(in)    :: n_s
    end subroutine Smoother_V

    !---------------------------------------------------------------------------
    !> Monitoring with given residual

    module subroutine Monitoring_R(this, l, step, r)
      class(ML_DG_EllipticSolver_3D), intent(in) :: this
      integer,               intent(in) :: l          !< level
      character(len=*),      intent(in) :: step       !< current step
      real(RNP), contiguous, intent(in) :: r(:,:,:,:) !< residual
    end subroutine Monitoring_R

    !---------------------------------------------------------------------------
    !> Monitoring with constant diffusivity

    module subroutine Monitoring_C(this, l, step, bc, lambda, nu, f, bv, u)
      class(ML_DG_EllipticSolver_3D),       intent(in) :: this
      integer,                              intent(in) :: l
      character(len=*),                     intent(in) :: step
      character,                            intent(in) :: bc(:)
      real(RNP),                            intent(in) :: lambda
      real(RNP),                            intent(in) :: nu
      real(RNP), contiguous,                intent(in) :: f(:,:,:,:)
      real(RNP), contiguous,                intent(in) :: u(:,:,:,:)
      class(BoundaryVariable_3D), optional, intent(in) :: bv(:)
    end subroutine Monitoring_C

    !---------------------------------------------------------------------------
    !> Monitoring with variable diffusivity

    module subroutine Monitoring_V(this, l, step, bc, lambda, nu, f, bv, u)
      class(ML_DG_EllipticSolver_3D),       intent(in) :: this
      integer,                              intent(in) :: l
      character(len=*),                     intent(in) :: step
      character,                            intent(in) :: bc(:)
      real(RNP),                            intent(in) :: lambda
      real(RNP), contiguous,                intent(in) :: nu(:,:,:,:)
      real(RNP), contiguous,                intent(in) :: f(:,:,:,:)
      real(RNP), contiguous,                intent(in) :: u(:,:,:,:)
      class(BoundaryVariable_3D), optional, intent(in) :: bv(:)
    end subroutine Monitoring_V

  end interface

contains

  !=============================================================================
  ! ML_DG_EllipticSolver_3D: constructor type bound procedures

  !-----------------------------------------------------------------------------
  !> Constructor of 3D FAS multigrid solver

  function New_ML_DG_EllipticSolver_3D(ml_op, opt) result(this)
    class(ML_MeshOperators_3D), target, intent(in) :: ml_op
    class(ML_DG_EllipticOptions_3D),    intent(in) :: opt

    type(ML_DG_EllipticSolver_3D) :: this

    call Init_ML_DG_EllipticSolver_3D(this, ml_op, opt)

  end function New_ML_DG_EllipticSolver_3D

  !-----------------------------------------------------------------------------
  !> Initialization of 3D FAS multigrid solver

  subroutine Init_ML_DG_EllipticSolver_3D(this, ml_op, opt)
    class(ML_DG_EllipticSolver_3D),     intent(inout) :: this
    class(ML_MeshOperators_3D), target, intent(in)    :: ml_op
    class(ML_DG_EllipticOptions_3D),    intent(in)    :: opt

    integer :: l, l_top

    l_top = size(ml_op % sem)

    allocate(this % elliptic_op(l_top))

    do l = 1, l_top

      this % elliptic_op(l) = DG_EllipticOperator_3D( ml_op % sem(l)    &
                                                    , opt % schwarz     &
                                                    , opt % penalty     &
                                                    , opt % interior_bc )
    end do

    this % ml_op => ml_op

    this % start_method   = opt % start_method
    this % smooth_method  = opt % smooth_method
    this % coarse_solver  = opt % coarse_solver
    this % fc_projection  = opt % fc_projection
    this % fc_restriction = opt % fc_restriction

    this % i_crs = opt % i_crs
    this % i_max = opt % i_max
    this % ns_0  = opt % ns_0
    this % ns_1  = opt % ns_1
    this % ns_2  = opt % ns_2

    if (opt % ns_c < 0) then
      this % ns_c = max(this%ns_1, this%ns_2)
    else
      this % ns_c = opt % ns_c
    end if

    if (opt % ns_f < 0) then
      this % ns_f = this % ns_2
    else
      this % ns_f = opt % ns_f
    end if

    this % vs_mode = opt % vs_mode
    this % vs_lmin = opt % vs_lmin
    if (opt % vs_lmax > 0) then
      this % vs_lmax = opt % vs_lmax
    else
      this % vs_lmax = l_top
    end if

    this % r_crs = opt % r_crs
    this % r_red = opt % r_red
    this % r_max = opt % r_max

    if (this % r_crs < 0) then
      this % r_crs = this % r_red
    end if

  end subroutine Init_ML_DG_EllipticSolver_3D

  !-----------------------------------------------------------------------------
  !> Number of pre-smoothing steps on level l accounting for possible variation

  pure integer function NumSmoothingSteps(this, level, stage) result(n)
    class(ML_DG_EllipticSolver_3D), intent(in)  :: this
    integer, intent(in) :: level !< current level
    integer, intent(in) :: stage !< 0/1/2: cascade/pre/post-smoothing

    integer :: l, m, p, q

    select case(stage)
    case(0)
      n = this % ns_0
    case(1)
      n = this % ns_1
    case(2)
      n = this % ns_2
    case default
      error stop 'NumSmoothingSteps: invalid stage'
    end select

    if (this % vs_mode == 0 .or. this % vs_lmax <= level) return

    l = max(this % vs_lmin, level)
    m = this % vs_lmax

    select case(this % vs_mode)
    case(1)
      ! polynomial coarsening
      n = n * this%elliptic_op(m)%eop%po / this%elliptic_op(l)%eop%po
    case(2)
      ! arbitrary coarsening
      n = n * 2 ** (m-l)
    end select

  end function NumSmoothingSteps

  !=============================================================================
  ! ML_DG_EllipticOptions_3D: constructor type bound procedures

  !-----------------------------------------------------------------------------
  !> Broadcasting multilevel mesh options

  subroutine Bcast_ML_DG_EllipticOptions_3D(this, root, comm)
    class(ML_DG_EllipticOptions_3D), intent(inout) :: this !< options
    integer       , intent(in) :: root !< rank root process
    type(MPI_Comm), intent(in) :: comm !< MPI communicator

    call XMPI_Bcast(this % penalty        , root, comm)
    call XMPI_Bcast(this % start_method   , root, comm)
    call XMPI_Bcast(this % smooth_method  , root, comm)
    call XMPI_Bcast(this % coarse_solver  , root, comm)
    call XMPI_Bcast(this % fc_projection  , root, comm)
    call XMPI_Bcast(this % fc_restriction , root, comm)
    call XMPI_Bcast(this % i_crs          , root, comm)
    call XMPI_Bcast(this % i_max          , root, comm)
    call XMPI_Bcast(this % ns_0           , root, comm)
    call XMPI_Bcast(this % ns_1           , root, comm)
    call XMPI_Bcast(this % ns_2           , root, comm)
    call XMPI_Bcast(this % ns_c           , root, comm)
    call XMPI_Bcast(this % ns_f           , root, comm)
    call XMPI_Bcast(this % vs_mode        , root, comm)
    call XMPI_Bcast(this % vs_lmin        , root, comm)
    call XMPI_Bcast(this % vs_lmax        , root, comm)
    call XMPI_Bcast(this % r_crs          , root, comm)
    call XMPI_Bcast(this % r_red          , root, comm)
    call XMPI_Bcast(this % r_max          , root, comm)

    call this % schwarz % Bcast(root, comm)

  end subroutine Bcast_ML_DG_EllipticOptions_3D

  !=============================================================================

end module ML__DG__Elliptic_Solver__3D
