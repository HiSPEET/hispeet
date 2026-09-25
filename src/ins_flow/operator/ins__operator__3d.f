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

!> summary:  Incompressible Navier-Stokes DG-SEM operator
!> author:   Joerg Stiller
!> date:     2021/12/30
!===============================================================================

module INS__Operator__3D
  use Kind_Parameters
  use Constants
  use Execution_Control
  use Array_Assignments
  use XMPI

  use TPO__INS_Convection__3D

  use Standard_Element_Operators__1D
  use Embedded_Interpolation_Operator__1D
  use Projection_Operator__1D
  use DG__Element_Operators__1D
  use DG__Elliptic_Operator__3D
  use DG__Schwarz_Operator__3D

  use Mesh__3D
  use Trace_Operators__3D
  use Spectral_Element_Mesh__3D
  use Boundary_Variable__3D

  use ML__Mesh_Variable__3D
  use ML__Boundary_Variable__3D
  use ML__DG__Elliptic_Solver__3D

  use INS__Problem__3D
  use INS__SGS__Model__3D

  implicit none
  private

  public :: INS_Operator_3D
  public :: INS_OperatorOptions_3D

  !-----------------------------------------------------------------------------
  !> DG-SEM mesh and operators for incompressible Navier-Stokes problems

  type INS_Operator_3D

    class(INS_Problem_3D), pointer :: problem => null() !< flow problem

    integer :: level !< rank in multilevel hierarchy (0 if none)

    character(len=4) :: convection_term  !< form of the convection term
    character(len=4) :: pressure_solver  !< pressure solver
    character(len=4) :: diffusion_solver !< diffusion solver
    character(len=1) :: interior_bc      !< coupling to frozen elements

    integer   :: div_stab  !< grad-div stabilization method
    real(RNP) :: mu_0      !< bulk viscosity,  μ = ζ/ρ (if const)
    real(RNP) :: nu_0      !< shear viscosity, ν = η/ρ (if const)
    real(RNP) :: delta_out !< δ parameter of outflow conditions

    ! element operators
    type(DG_ElementOperators_1D)      :: eop_u !< DG operators for u \ p
    type(DG_ElementOperators_1D)      :: eop_p !< DG operators for p
    type(StandardElementOperators_1D) :: sop_q !< quadrature ops for convection

    ! projection and interpolation operators
    type(ProjectionOperator_1D)            :: pop_up !< u to p L² projection
    type(EmbeddedInterpolationOperator_1D) :: iop_up !< u to p interpolation
    type(EmbeddedInterpolationOperator_1D) :: iop_pu !< p to u interpolation
    type(EmbeddedInterpolationOperator_1D) :: iop_uq !< u to q interpolation

    ! mesh and metrics
    type(Mesh_3D), pointer                    :: mesh  !< local mesh partition
    type(SpectralElementMesh_3D), pointer     :: sem_u !< mesh + metrics for u
    type(SpectralElementMesh_3D), pointer     :: sem_p !< mesh + metrics for p
    type(SpectralElementMesh_3D), allocatable :: sem_q !< quadrature operators

    ! operators and solvers for elliptic subsystems
    type(DG_EllipticOperator_3D) :: elliptic_p !< elliptic operator for p
    type(DG_SchwarzOperator_3D)  :: schwarz_u  !< Schwarz operators for u
    type(ML_DG_EllipticSolver_3D), pointer :: ml_solver_p ! ML pessure solver

    type(INS_SGS_Model_3D) :: sgs_model !< subgrid-scale model

    ! iterative solver settings
    integer   :: i_max_p   !< max num p-iterations in projection solver
    integer   :: i_max_v   !< max num v-iterations in projection solver
    integer   :: k_max     !< max num Krylov iterations
    integer   :: k_pre_p   !< max num p-iterations in Krylov preconditioner
    integer   :: k_pre_v   !< max num v-iterations in Krylov preconditioner
    real(RNP) :: r_red     !< min residual reduction, if > 0
    real(RNP) :: r_max     !< max residual to reach,  if > 0

  contains

    generic   :: Init => Init_INS_Operator_3D
    procedure :: Init_INS_Operator_3D

    procedure :: HasVariableViscosity

    procedure :: ApplyEssentialBC
    procedure :: ApplyNaturalBC

    procedure :: ApplyDiffusionOperator
    procedure :: ApplyDiffusionOperator_C
    procedure :: ApplyDiffusionOperator_V
    procedure :: ApplyStokesOperator

    procedure :: DiffusionSolver

    procedure :: GetBackflowPenalty
    procedure :: GetConvectionTerm

    procedure :: GetDiffusionResidual
    procedure :: GetDiffusionResidual_C
    procedure :: GetDiffusionResidual_V

    procedure :: GetDiffusionTerm
    procedure :: GetDiffusionTerm_C
    procedure :: GetDiffusionTerm_V

    procedure :: GetPressureBoundaryValues
    procedure :: GetStokesResidual
    procedure :: GetVariableViscosity

    procedure :: GetViscousBoundaryStress
    procedure :: GetViscousBoundaryStress_C
    procedure :: GetViscousBoundaryStress_V

    procedure :: PressureSolver
    procedure :: StokesSolver

  end type INS_Operator_3D

  ! constructor interface
  interface INS_Operator_3D
    procedure New_INS_Operator_3D
  end interface

  !-----------------------------------------------------------------------------
  !> Options for INS_Operator_3D initialization

  type INS_OperatorOptions_3D

    character(4) :: convection_term  = 'flux' !< {'flux','skew','conv'}
    character(4) :: pressure_solver  = 'SPCG' !< {'AS','CG','SPCG','MG','MGCG'}
    character(4) :: diffusion_solver = 'DPCG' !< {'DPCG','SPCG'}
    character(1) :: interior_bc      = ' '    !< {' ','D','M'}

    logical   :: dealiasing = .false. !< F: no dealiasing, T: 3/2 rule

    real(RNP) :: penalty_u = -1 !< penalty for u-solver, -1: auto
    real(RNP) :: penalty_p = -1 !< penalty for p-solver, -1: auto
    integer   :: div_stab  =  1 !< grad-div stabilization method
      !!  - `1`  constant, μ₁ = μ₀
      !!  - `2`  constant, μ₂ = μ₀ max(ν, vh), using reference values for ν and v
      !!  - `3`  variable, μ₃ = μ₀ max(ν, vh), using element averages for ν and v
      !!  - `4`  variable, μ₄ = μ₂ averaged across faces
      !!  - `5`  variable, μ₅ = μ₃ averaged across faces

    real(RNP) :: mu_0      =     0 !< bulk viscosity coefficient μ₀
    real(RNP) :: delta_out = 1e-02 !< outflow parameter

    integer   :: i_max_p   =  1000 !< max num p-iterations in projection
    integer   :: i_max_v   =   200 !< max num v-iterations in projection
    integer   :: k_max     =     0 !< max num Krylov iterations
    integer   :: k_pre_p   =    10 !< max num p-iterations in Krylov precon
    integer   :: k_pre_v   =    10 !< max num v-iterations in Krylov precon
    real(RNP) :: r_red     = 1e-08 !< min residual reduction, if > 0
    real(RNP) :: r_max     = 1e-12 !< max residual to reach,  if > 0

    type(DG_SchwarzOptions_3D) :: schwarz_u !< Schwarz options for u-solver
    type(DG_SchwarzOptions_3D) :: schwarz_p !< Schwarz options for p-solver
    type(INS_SGS_Options_3D)   :: sgs_model !< subgrid-scale model

  contains
    procedure :: Bcast => Bcast_INS_OperatorOptions_3D
  end type INS_OperatorOptions_3D

  !=============================================================================
  ! Module procedures (alphabetically)

  interface

    !---------------------------------------------------------------------------
    !> Homogeneous diffusion operator with constant viscosity

    module subroutine ApplyDiffusionOperator_C(this, tau, bv, v, r, form)
      class(INS_Operator_3D),               intent(in)  :: this
      real(RNP),                            intent(in)  :: tau
      class(BoundaryVariable_3D), optional, intent(in)  :: bv(:)
      real(RNP), contiguous,                intent(in)  :: v(:,:,:,:,:)
      real(RNP), contiguous,                intent(out) :: r(:,:,:,:,:)
      integer,                    optional, intent(in)  :: form
    end subroutine ApplyDiffusionOperator_C

    !---------------------------------------------------------------------------
    !> Homogeneous diffusion operator with variable viscosity

    module subroutine ApplyDiffusionOperator_V(this, tau, mu, nu, bv, v, r, form)
      class(INS_Operator_3D),               intent(in)  :: this
      real(RNP),                            intent(in)  :: tau
      real(RNP), contiguous,                intent(in)  :: mu(:,:,:,:)
      real(RNP), contiguous,                intent(in)  :: nu(:,:,:,:)
      class(BoundaryVariable_3D), optional, intent(in)  :: bv(:)
      real(RNP), contiguous,                intent(in)  :: v(:,:,:,:,:)
      real(RNP), contiguous,                intent(out) :: r(:,:,:,:,:)
      integer,                    optional, intent(in)  :: form
    end subroutine ApplyDiffusionOperator_V

    !---------------------------------------------------------------------------
    !> Application of velocity boundary conditions to trace variables

    module subroutine ApplyEssentialBC(this, bv_u, vm, vp)
      class(INS_Operator_3D),               intent(in)    :: this
      class(BoundaryVariable_3D), optional, intent(in)    :: bv_u(:)
      real(RNP), contiguous,                intent(in)    :: vm(:,:,:,:,:)
      real(RNP), contiguous,                intent(inout) :: vp(:,:,:,:,:)
    end subroutine ApplyEssentialBC

    !---------------------------------------------------------------------------
    !> Application of natural boundary conditions to velocity trace variables

    module subroutine ApplyNaturalBC(this, bv_s, sm, sp)
      class(INS_Operator_3D),               intent(in)    :: this
      class(BoundaryVariable_3D), optional, intent(in)    :: bv_s(:)
      real(RNP), contiguous,                intent(in)    :: sm(:,:,:,:,:)
      real(RNP), contiguous,                intent(inout) :: sp(:,:,:,:,:)
    end subroutine ApplyNaturalBC

    !---------------------------------------------------------------------------
    !> Application of the Stokes operator

    module subroutine ApplyStokesOperator(this, tau, bv, mu, nu, u, r)
      class(INS_Operator_3D),               intent(in)  :: this
      real(RNP),                            intent(in)  :: tau
      class(BoundaryVariable_3D), optional, intent(in)  :: bv(:)
      real(RNP), contiguous,      optional, intent(in)  :: mu(:,:,:,:)
      real(RNP), contiguous,      optional, intent(in)  :: nu(:,:,:,:)
      real(RNP), contiguous,                intent(in)  :: u(:,:,:,:,:)
      real(RNP), contiguous,                intent(out) :: r(:,:,:,:,:)
    end subroutine ApplyStokesOperator

    !---------------------------------------------------------------------------
    !> Diffusion solver

    module subroutine DiffusionSolver( this, tau, mu, nu, f, bv, v, i_max &
                                     , precon, ni )
      class(INS_Operator_3D),          intent(in)    :: this
      real(RNP),                       intent(in)    :: tau
      real(RNP), contiguous, optional, intent(in)    :: mu(:,:,:,:)
      real(RNP), contiguous, optional, intent(in)    :: nu(:,:,:,:)
      real(RNP), contiguous,           intent(in)    :: f(:,:,:,:,:)
      class(BoundaryVariable_3D),      intent(in)    :: bv(:)
      real(RNP), contiguous,           intent(inout) :: v(:,:,:,:,:)
      integer,               optional, intent(in)    :: i_max
      logical,               optional, intent(in)    :: precon
      integer,               optional, intent(out)   :: ni
    end subroutine DiffusionSolver

    !---------------------------------------------------------------------------
    !> Backflow pressure penalty

    module subroutine GetBackflowPenalty(this, problem, b, v, dp)
      class(INS_Operator_3D), intent(in)  :: this
      class(INS_Problem_3D),  intent(in)  :: problem
      integer,                intent(in)  :: b
      real(RNP),  contiguous, intent(in)  :: v(:,:,:,:,:)
      real(RNP),  contiguous, intent(out) :: dp(:,:,:)
    end subroutine GetBackflowPenalty

    !---------------------------------------------------------------------------
    !> Diffusion residual with constant viscosity

    module subroutine GetDiffusionResidual_C(this, tau, f, bv, v, r, form)
      class(INS_Operator_3D),     intent(in)  :: this
      real(RNP),                  intent(in)  :: tau
      real(RNP), contiguous,      intent(in)  :: f(:,:,:,:,:)
      class(BoundaryVariable_3D), intent(in)  :: bv(:)
      real(RNP), contiguous,      intent(in)  :: v(:,:,:,:,:)
      real(RNP), contiguous,      intent(out) :: r(:,:,:,:,:)
      integer,          optional, intent(in)  :: form
    end subroutine GetDiffusionResidual_C

    !---------------------------------------------------------------------------
    !> Diffusion residual with variable viscosity

    module subroutine GetDiffusionResidual_V &
        (this, tau, mu, nu, f, bv, v, r, form)
      class(INS_Operator_3D),     intent(in)  :: this
      real(RNP),                  intent(in)  :: tau
      real(RNP), contiguous,      intent(in)  :: mu(:,:,:,:)
      real(RNP), contiguous,      intent(in)  :: nu(:,:,:,:)
      real(RNP), contiguous,      intent(in)  :: f(:,:,:,:,:)
      class(BoundaryVariable_3D), intent(in)  :: bv(:)
      real(RNP), contiguous,      intent(in)  :: v(:,:,:,:,:)
      real(RNP), contiguous,      intent(out) :: r(:,:,:,:,:)
      integer,          optional, intent(in)  :: form
    end subroutine GetDiffusionResidual_V

    !---------------------------------------------------------------------------
    !> Diffusion term with constant viscosity on irregular (deformed) mesh

    module subroutine GetDiffusionTerm_C(this, v, vp, sp, f_d, bv_v, xout, form)
      class(INS_Operator_3D),               intent(in)  :: this
      real(RNP),                contiguous, intent(in)  :: v(:,:,:,:,:)
      real(RNP),                contiguous, intent(out) :: vp(:,:,:,:,:)
      real(RNP),                contiguous, intent(out) :: sp(:,:,:,:,:)
      real(RNP),                contiguous, intent(out) :: f_d(:,:,:,:,:)
      class(BoundaryVariable_3D), optional, intent(in)  :: bv_v(:)
      logical,                    optional, intent(in)  :: xout
      integer,                    optional, intent(in)  :: form
    end subroutine GetDiffusionTerm_C

    !---------------------------------------------------------------------------
    !> Diffusion term with variable viscosity on irregular (deformed) mesh

    module subroutine GetDiffusionTerm_V &
        (this, mu, nu, v, vp, sp, f_d, bv_v, xout, form)
      class(INS_Operator_3D),               intent(in)  :: this
      real(RNP), contiguous,                intent(in)  :: mu(:,:,:,:)
      real(RNP), contiguous,                intent(in)  :: nu(:,:,:,:)
      real(RNP), contiguous,                intent(in)  :: v(:,:,:,:,:)
      real(RNP), contiguous,                intent(out) :: vp(:,:,:,:,:)
      real(RNP), contiguous,                intent(out) :: sp(:,:,:,:,:)
      real(RNP), contiguous,                intent(out) :: f_d(:,:,:,:,:)
      class(BoundaryVariable_3D), optional, intent(in)  :: bv_v(:)
      logical,                    optional, intent(in)  :: xout
      integer,                    optional, intent(in)  :: form
    end subroutine GetDiffusionTerm_V

    !---------------------------------------------------------------------------
    !> Provision of pressure boundary values

    module subroutine GetPressureBoundaryValues(this, tau, v, bv_u, bv_p, bv_q)
      class(INS_Operator_3D),               intent(in)    :: this
      real(RNP),                            intent(in)    :: tau
      real(RNP), contiguous,                intent(in)    :: v(:,:,:,:,:)
      class(BoundaryVariable_3D),           intent(in)    :: bv_u(:)
      class(BoundaryVariable_3D),           intent(inout) :: bv_p(:)
      class(BoundaryVariable_3D), optional, intent(inout) :: bv_q(:)
    end subroutine GetPressureBoundaryValues

    !---------------------------------------------------------------------------
    !> Stokes residual for incompressible flow

    module subroutine GetStokesResidual(this, tau, f, bv, mu, nu, u, r)
      class(INS_Operator_3D),          intent(in)  :: this
      real(RNP),                       intent(in)  :: tau
      real(RNP), contiguous,           intent(in)  :: f(:,:,:,:,:)
      class(BoundaryVariable_3D),      intent(in)  :: bv(:)
      real(RNP), contiguous, optional, intent(in)  :: mu(:,:,:,:)
      real(RNP), contiguous, optional, intent(in)  :: nu(:,:,:,:)
      real(RNP), contiguous,           intent(in)  :: u(:,:,:,:,:)
      real(RNP), contiguous,           intent(out) :: r(:,:,:,:,:)
    end subroutine GetStokesResidual

    !---------------------------------------------------------------------------
    !> Variable viscosity coefficients

    module subroutine GetVariableViscosity(this, t, u, mu, nu)
      class(INS_Operator_3D), intent(in)  :: this
      real(RNP),              intent(in)  :: t
      real(RNP),  contiguous, intent(in)  :: u(:,:,:,:,:)
      real(RNP),  contiguous, intent(out) :: mu(:,:,:,:)
      real(RNP),  contiguous, intent(out) :: nu(:,:,:,:)
    end subroutine GetVariableViscosity

    !---------------------------------------------------------------------------
    !> Viscous stress vector on a boundary (C)

    module subroutine GetViscousBoundaryStress_C &
        (this, b, v, sb, bv, xout, form)
      class(INS_Operator_3D),               intent(in)  :: this
      integer,                              intent(in)  :: b
      real(RNP), contiguous,                intent(in)  :: v(:,:,:,:,:)
      real(RNP), contiguous,                intent(out) :: sb(:,:,:,:)
      class(BoundaryVariable_3D), optional, intent(in)  :: bv(:)
      logical,                    optional, intent(in)  :: xout
      integer,                    optional, intent(in)  :: form
    end subroutine GetViscousBoundaryStress_C

    !---------------------------------------------------------------------------
    !> Viscous stress vector on a boundary (V)

    module subroutine GetViscousBoundaryStress_V &
        (this, b, mu, nu, v, sb, bv, xout, form)
      class(INS_Operator_3D),               intent(in)  :: this
      integer,                              intent(in)  :: b
      real(RNP), contiguous,                intent(in)  :: mu(:,:,:,:)
      real(RNP), contiguous,                intent(in)  :: nu(:,:,:,:)
      real(RNP), contiguous,                intent(in)  :: v(:,:,:,:,:)
      real(RNP), contiguous,                intent(out) :: sb(:,:,:,:)
      class(BoundaryVariable_3D), optional, intent(in)  :: bv(:)
      logical,                    optional, intent(in)  :: xout
      integer,                    optional, intent(in)  :: form
    end subroutine GetViscousBoundaryStress_V

    !---------------------------------------------------------------------------
    !> Projection-based pressure solver

    module subroutine PressureSolver(this, tau, bv_u, v, f, p, precon, ni)
      class(INS_Operator_3D),     intent(in)    :: this
      real(RNP),                  intent(in)    :: tau
      class(BoundaryVariable_3D), intent(in)    :: bv_u(:)
      real(RNP), contiguous,      intent(in)    :: v(:,:,:,:,:)
      real(RNP), contiguous,      intent(in)    :: f(:,:,:,:)
      real(RNP), contiguous,      intent(inout) :: p(:,:,:,:)
      logical,          optional, intent(in)    :: precon
      integer,          optional, intent(out)   :: ni
    end subroutine PressureSolver

    !---------------------------------------------------------------------------
    !> Projection-diffusion step for incompressible flow

    module subroutine StokesProjection &
        (this, tau, f, bv, mu, nu, u, f_d0, precon)
      class(INS_Operator_3D),          intent(in)    :: this
      real(RNP),                       intent(in)    :: tau
      real(RNP), contiguous,           intent(in)    :: f(:,:,:,:,:)
      class(BoundaryVariable_3D),      intent(in)    :: bv(:)
      real(RNP), contiguous, optional, intent(in)    :: mu(:,:,:,:)
      real(RNP), contiguous, optional, intent(in)    :: nu(:,:,:,:)
      real(RNP), contiguous,           intent(inout) :: u(:,:,:,:,:)
      real(RNP), contiguous, optional, intent(in)    :: f_d0(:,:,:,:,:)
      logical,               optional, intent(in)    :: precon
    end subroutine StokesProjection

    !---------------------------------------------------------------------------
    !> FGMRES for Stokes part with projection-diffusion preconditioner

    module subroutine StokesFGMRES(this, tau, f, bv, mu, nu, u, f_d0)
      class(INS_Operator_3D),          intent(in)    :: this
      real(RNP),                       intent(in)    :: tau
      real(RNP), contiguous,           intent(in)    :: f(:,:,:,:,:)
      class(BoundaryVariable_3D),      intent(in)    :: bv(:)
      real(RNP), contiguous, optional, intent(in)    :: mu(:,:,:,:)
      real(RNP), contiguous, optional, intent(in)    :: nu(:,:,:,:)
      real(RNP), contiguous,           intent(inout) :: u(:,:,:,:,:)
      real(RNP), contiguous, optional, intent(in)    :: f_d0(:,:,:,:,:)
    end subroutine StokesFGMRES

  end interface

contains

  !=============================================================================
  ! Type-bound procedures of INS_Operator_3D

  !-----------------------------------------------------------------------------
  !> Constructor of INS_Operator_3D

  function New_INS_Operator_3D(opt, problem, sem_u, sem_p, ml_solver_p, level) &
        result(this)

    class(INS_OperatorOptions_3D), intent(in) :: opt
      !< INS operator options
    class(INS_Problem_3D), intent(in) :: problem
      !< INS flow problem
    type(SpectralElementMesh_3D), target, intent(in) :: sem_u
      !< spectral element mesh and operators for u
    type(SpectralElementMesh_3D),  optional, target, intent(in) :: sem_p
      !< spectral element mesh and operators for p, if different from `sem_u`
    type(ML_DG_EllipticSolver_3D), optional, target, intent(in) :: ml_solver_p
      !< multilevel pressure solver [none]
    integer, optional, intent(in) :: level
      !< rank in multilevel hierarchy (0 if none) 0[]
    type(INS_Operator_3D) :: this

    call Init_INS_Operator_3D( this, opt, problem, sem_u, sem_p &
                             , ml_solver_p, level               )

  end function New_INS_Operator_3D

  !-----------------------------------------------------------------------------
  !> Initialization of INS_Operator_3D

  subroutine Init_INS_Operator_3D &
               (this, opt, problem, sem_u, sem_p, ml_solver_p, level)

    class(INS_Operator_3D), intent(inout) :: this
      !< new INS operator
    class(INS_OperatorOptions_3D), intent(in) :: opt
      !< INS operator options
    class(INS_Problem_3D), target, intent(in) :: problem
      !< INS flow problem
    type(SpectralElementMesh_3D), target, intent(in) :: sem_u
      !< spectral element mesh and operators for u
    type(SpectralElementMesh_3D),  optional, target, intent(in) :: sem_p
      !< spectral element mesh and operators for p, if different from `sem_u`
    type(ML_DG_EllipticSolver_3D), optional, target, intent(in) :: ml_solver_p
      !< multilevel pressure solver [none]
    integer, optional, intent(in) :: level
      !< rank in multilevel hierarchy (0 if none) 0[]

    ! problem ..................................................................

    this % problem => problem

    ! parameters ...............................................................

    if (present(level)) then
      this % level = level
    else
      this % level = 0
    end if

    select case(opt % convection_term)
    case('flux','skew','conv')
      this % convection_term = opt % convection_term
    case default
      call Error( 'Init_INS_Operator_3D'                                     &
                , 'ivalid convection term "'//trim(opt%convection_term)//'"' &
                , 'INS__Operator__3D'                                        )
    end select
    select case(opt % pressure_solver)
    case('AS','CG','SPCG','MG','MGCG')
      this % pressure_solver = opt % pressure_solver
    case default
      call Error( 'Init_INS_Operator_3D'                                     &
                , 'ivalid pressure solver "'//trim(opt%pressure_solver)//'"' &
                , 'INS__Operator__3D'                                        )
    end select

    select case(opt % diffusion_solver)
    case('DPCG','SPCG')
      this % diffusion_solver = opt % diffusion_solver
    case default
      call Error( 'Init_INS_Operator_3D'                                       &
                , 'ivalid diffusion solver "'//trim(opt%diffusion_solver)//'"' &
                , 'INS__Operator__3D'                                          )
    end select

    this % div_stab  = opt % div_stab
    this % mu_0      = opt % mu_0
    this % nu_0      = problem % nu_ref
    this % delta_out = opt % delta_out

    ! element operators ........................................................

    this % eop_u = DG_ElementOperators_1D(sem_u % std_op, opt % penalty_u)
    this % eop_p = DG_ElementOperators_1D(sem_p % std_op, opt % penalty_p)

    if (this%eop_u%nodes /= 'L' .or. this%eop_p%nodes /= 'L') then
      call Error( 'Init_INS_Operator_3D'               &
                , 'Lobatto nodes required for u and p' &
                , 'INS__Operator__3D'                  )
    end if

    if (opt % dealiasing) then
      ! use 3/2 rule for dealiasing
      this % sop_q = StandardElementOperators_1D &
                         (po = ceiling(1.5 * this%eop_u%po), no_vdm = .true.)
    else
      ! use velocity Lobatto points for for convection
      this % sop_q = sem_u % std_op
    end if

    ! projection and interpolation operators ...................................

    this % pop_up = ProjectionOperator_1D(this%eop_p, this%eop_u%x, nodes='L')
    this % iop_pu = EmbeddedInterpolationOperator_1D( this%eop_p, this%eop_u%x )
    this % iop_up = EmbeddedInterpolationOperator_1D( this%eop_u, this%eop_p%x )
    this % iop_uq = EmbeddedInterpolationOperator_1D( this%eop_u, this%sop_q%x )

    ! mesh and metrics .........................................................

    this % mesh  => sem_u % mesh
    this % sem_u => sem_u

    if (present(sem_p)) then
      this % sem_p => sem_p
    else
      this % sem_p => sem_u
    end if

    if (opt % dealiasing) then
      this % sem_q = SpectralElementMesh_3D(this % mesh, this % sop_q % po)
    else
      this % sem_q = this % sem_u
    end if

    ! operators and solvers for elliptic subsystems ............................

    ! elliptic operator for pressure
    this % elliptic_p = DG_EllipticOperator_3D( sem_p             &
                                              , opt % schwarz_p   &
                                              , opt % penalty_p   &
                                              , opt % interior_bc )

    ! Schwarz operators for viscous diffusion
    this % schwarz_u = DG_SchwarzOperator_3D( opt  % schwarz_u &
                                            , this % eop_u     &
                                            , this % mesh      )

    ! multilevel pressure solver
    if (present(ml_solver_p)) then
      this % ml_solver_p => ml_solver_p
    else
      this % ml_solver_p => null()
    end if

    ! subgrid-scale model
    this % sgs_model = INS_SGS_Model_3D(opt % sgs_model)

    ! iterative solver settings
    this % i_max_p = opt % i_max_p
    this % i_max_v = opt % i_max_v
    this % k_max   = opt % k_max
    this % k_pre_p = opt % k_pre_p
    this % k_pre_v = opt % k_pre_v
    this % r_red   = opt % r_red
    this % r_max   = opt % r_max

  end subroutine Init_INS_Operator_3D

  !-----------------------------------------------------------------------------
  !> Function for enquiring wether viscosity coefficients are variable

  logical function HasVariableViscosity(this) result(hvv)
    class(INS_Operator_3D), intent(in) :: this

    ! variable bulk viscosity
    hvv = this % div_stab > 1

    ! problem dependent shear viscosity
    hvv = hvv .or. this % problem % HasVariableProperties()

    ! variation due to SGS model
    hvv = hvv .or. this % sgs_model % model > 0

  end function HasVariableViscosity

  !-----------------------------------------------------------------------------
  !> Diffusion term with constant viscosity

  subroutine GetConvectionTerm(this, v, vp, f_c, frozen)
    class(INS_Operator_3D), intent(in) :: this
    real(RNP), contiguous, intent(in)  :: v(:,:,:,:,:)   !< velocity
    real(RNP), contiguous, intent(in)  :: vp(:,:,:,:,:)  !< outer velocity v⁺
    real(RNP), contiguous, intent(out) :: f_c(:,:,:,:,:) !< convection term

    logical, optional, intent(in) :: frozen !< T/F in/exclude frozen elements [F]

    character, allocatable, save :: bc_elem(:,:)

    character :: bc
    integer   :: b, e, i, j, k, ne, e_max

    ! range ....................................................................

    ne = this % mesh % n_elem

    e_max = ne
    if (present(frozen)) then
      if (frozen) then
        e_max = this % mesh % n_elem_active
      end if
    end if

    ! element face boundary conditions .........................................

    !$omp master
    allocate(bc_elem(6,this%mesh%n_elem), source = ' ')
    do b = 1, this % mesh % n_bound
      bc = this % problem % bc_v(b)
      select case(bc)
      case('D')
        associate(bface => this % mesh % boundary(b) % face)
          do k = 1, this % mesh % boundary(b) % n_face
            bc_elem(bface(k)%element_face, bface(k)%element_id) = bc
          end do
        end associate
      end select
    end do
    !$omp end master
    !$omp barrier

    ! evaluation ...............................................................

    call TPO_INS_Convection( nv    = this % eop_u % po + 1       &
                           , nq    = this % sop_q % po + 1       &
                           , ne    = ne                          &
                           , bc    = bc_elem                     &
                           , D_v   = this % eop_u  % D           &
                           , I_vq  = this % iop_uq % A           &
                           , w_q   = this % sop_q  % w           &
                           , Jd_q  = this % sem_q % metrics % Jd &
                           , Ji_q  = this % sem_q % metrics % Ji &
                           , a_q   = this % sem_q % metrics % a  &
                           , n_q   = this % sem_q % metrics % n  &
                           , v     = v                           &
                           , vp    = vp                          &
                           , f_c   = f_c                         &
                           , form  = this % convection_term      &
                           , e_max = e_max                       )

    ! Coriolis force in rotating frame .........................................

    if (any(this % problem % omega /= 0)) then

      associate( omega => this % problem % omega      &
               , po    => this % eop_u % po           &
               , w     => this % eop_u % w            &
               , Jd    => this % sem_u % metrics % Jd )

        block
          real(RNP) :: Ms(0:po,0:po,0:po)
          real(RNP) :: ce(0:po,0:po,0:po)

          do k = 0, po
          do j = 0, po
          do i = 0, po
            Ms(i,j,k) = w(i) * w(j) * w(k)
          end do
          end do
          end do

          !$omp do
          do e = 1, e_max

            ce = -2 * Ms * Jd(:,:,:,e)

            f_c(:,:,:,e,1) = f_c(:,:,:,e,1) + ce * ( omega(2) * v(:,:,:,e,3) &
                                                   - omega(3) * v(:,:,:,e,2) )
            f_c(:,:,:,e,2) = f_c(:,:,:,e,2) + ce * ( omega(3) * v(:,:,:,e,1) &
                                                   - omega(1) * v(:,:,:,e,3) )
            f_c(:,:,:,e,3) = f_c(:,:,:,e,3) + ce * ( omega(1) * v(:,:,:,e,2) &
                                                   - omega(2) * v(:,:,:,e,1) )
          end do

        end block
      end associate
    end if

    ! finalization .............................................................

    !$omp master
    deallocate(bc_elem)
    !$omp end master

  end subroutine GetConvectionTerm

  !-----------------------------------------------------------------------------
  !> Diffusion operator with constant or variable viscosity

  subroutine ApplyDiffusionOperator(this, tau, mu, nu, bv, v, r, form)
    class(INS_Operator_3D),               intent(in)  :: this
    real(RNP),                            intent(in)  :: tau
    real(RNP), contiguous,      optional, intent(in)  :: mu(:,:,:,:)
    real(RNP), contiguous,      optional, intent(in)  :: nu(:,:,:,:)
    class(BoundaryVariable_3D), optional, intent(in)  :: bv(:)
    real(RNP), contiguous,                intent(in)  :: v(:,:,:,:,:)
    real(RNP), contiguous,                intent(out) :: r(:,:,:,:,:)
    integer,                    optional, intent(in)  :: form

    if (.not. this % HasVariableViscosity()) then
      call this % ApplyDiffusionOperator_C(tau, bv, v, r, form)

    else if (present(mu) .and. present(nu)) then
      call this % ApplyDiffusionOperator_V(tau, mu, nu, bv, v, r, form)

    else
      call Error( 'ApplyDiffusionOperator'                     &
                , 'mu and nu required with variable viscosity' &
                , 'INS_Operator_3D'                            )
    end if

  end subroutine ApplyDiffusionOperator

  !-----------------------------------------------------------------------------
  !> Diffusion residual with constant or variable viscosity

  subroutine GetDiffusionResidual(this, tau, mu, nu, f, bv, v, r, form)
    class(INS_Operator_3D),          intent(in)  :: this
    real(RNP),                       intent(in)  :: tau
    real(RNP), contiguous, optional, intent(in)  :: mu(:,:,:,:)
    real(RNP), contiguous, optional, intent(in)  :: nu(:,:,:,:)
    real(RNP), contiguous,           intent(in)  :: f(:,:,:,:,:)
    class(BoundaryVariable_3D),      intent(in)  :: bv(:)
    real(RNP), contiguous,           intent(in)  :: v(:,:,:,:,:)
    real(RNP), contiguous,           intent(out) :: r(:,:,:,:,:)
    integer,               optional, intent(in)  :: form

    if (.not. this % HasVariableViscosity()) then
      call this % GetDiffusionResidual_C(tau, f, bv, v, r, form)

    else if (present(mu) .and. present(nu)) then
      call this % GetDiffusionResidual_V(tau, mu, nu, f, bv, v, r, form)

    else
      call Error( 'GetDiffusionResidual'                     &
                , 'mu and nu required with variable viscosity' &
                , 'INS_Operator_3D'                            )
    end if

  end subroutine GetDiffusionResidual

  !-----------------------------------------------------------------------------
  !> Diffusion term with constant or variable viscosity

  subroutine GetDiffusionTerm(this, mu, nu, v, vp, sp, f_d, bv, xout, form)
    class(INS_Operator_3D),               intent(in)  :: this
    real(RNP), contiguous,      optional, intent(in)  :: mu(:,:,:,:)
    real(RNP), contiguous,      optional, intent(in)  :: nu(:,:,:,:)
    real(RNP), contiguous,                intent(in)  :: v(:,:,:,:,:)
    real(RNP), contiguous,                intent(out) :: vp(:,:,:,:,:)
    real(RNP), contiguous,                intent(out) :: sp(:,:,:,:,:)
    real(RNP), contiguous,                intent(out) :: f_d(:,:,:,:,:)
    class(BoundaryVariable_3D), optional, intent(in)  :: bv(:)
    logical,                    optional, intent(in)  :: xout
    integer,                    optional, intent(in)  :: form

    if (.not. this % HasVariableViscosity()) then
      call this % GetDiffusionTerm_C(v, vp, sp, f_d, bv, xout, form)

    else if (present(mu) .and. present(nu)) then
      call this % GetDiffusionTerm_V(mu, nu, v, vp, sp, f_d, bv, xout, form)

    else
      call Error( 'GetDiffusionTerm'                           &
                , 'mu and nu required with variable viscosity' &
                , 'INS_Operator_3D'                            )
    end if

  end subroutine GetDiffusionTerm

  !-----------------------------------------------------------------------------
  !> Viscous stress vector with constant or variable viscosity on a boundary

  subroutine GetViscousBoundaryStress(this, b, mu, nu, v, sb, bv, xout, form)
    class(INS_Operator_3D),               intent(in)  :: this
    integer,                              intent(in)  :: b
    real(RNP), contiguous,      optional, intent(in)  :: mu(:,:,:,:)
    real(RNP), contiguous,      optional, intent(in)  :: nu(:,:,:,:)
    real(RNP), contiguous,                intent(in)  :: v(:,:,:,:,:)
    real(RNP), contiguous,                intent(out) :: sb(:,:,:,:)
    class(BoundaryVariable_3D), optional, intent(in)  :: bv(:)
    logical,                    optional, intent(in)  :: xout
    integer,                    optional, intent(in)  :: form

    if (.not. this % HasVariableViscosity()) then
      call this % GetViscousBoundaryStress_C(b, v, sb, bv, xout, form)

    else if (present(mu) .and. present(nu)) then
      call this % GetViscousBoundaryStress_V(b, mu, nu, v, sb, bv, xout, form)

    else
      call Error( 'GetViscousBoundaryStress'                   &
                , 'mu and nu required with variable viscosity' &
                , 'INS_Operator_3D'                            )
    end if

  end subroutine GetViscousBoundaryStress

  !-----------------------------------------------------------------------------
  !> Solver for Stokes part using either projection step or FGMRES
  !>
  !> If `f_d0` is present, the velocity is initialized as `v = τ(f + f_d0)`.
  !> Otherwise, the given `u = [v,p]` is used as the initial approximation.

  subroutine StokesSolver(this, tau, f, bv, mu, nu, u, f_d0)
    class(INS_Operator_3D), intent(in) :: this
    !< incompressible Navier-Stokes time integrator
    real(RNP), intent(in) :: tau
    !< τ, effective time step width
    real(RNP), contiguous, intent(in) :: f(:,:,:,:,:)
    !< unweighted RHS: f = v₀/τ + f_c + f_s + ...
    class(BoundaryVariable_3D), intent(in) :: bv(:)
    !< boundary values
    !!   - Γᴰ :  [ v₁, v₂, v₃, - ]
    !!   - Γᴼ :  [ - , - , ∆p, p ]
    real(RNP), contiguous, optional, intent(inout) :: mu(:,:,:,:)
    !< μ, kinematic bulk viscosity
    real(RNP), contiguous, optional, intent(inout) :: nu(:,:,:,:)
    !< ν, kinematic shear viscosity
    real(RNP), contiguous, intent(inout) :: u(:,:,:,:,:)
    !< u = [v, p], velocity and pressure
    real(RNP), contiguous, optional, intent(in) :: f_d0(:,:,:,:,:)
    !< approximate diffusion term

    if (.not. this % HasVariableViscosity()) then

      if (this % k_max > 0) then
        call StokesFGMRES(this, tau, f, bv, null(), null(), u, f_d0)
      else
        call StokesProjection(this, tau, f, bv, null(), null(), u, f_d0)
      end if

    else if (present(mu) .and. present(nu)) then

      if (this % k_max > 0) then
        call StokesFGMRES(this, tau, f, bv, mu, nu, u, f_d0)
      else
        call StokesProjection(this, tau, f, bv, mu, nu, u, f_d0)
      end if

    else
      call Error( 'StokesSolver'                               &
                , 'mu and nu required with variable viscosity' &
                , 'INS_Operator_3D'                            )
    end if

  end subroutine StokesSolver

  !=============================================================================
  ! Type-bound procedures of INS_OperatorOptions_3D

  !-----------------------------------------------------------------------------
  !> MPI_Bcast for objects of type INS_OperatorOptions_3D

   subroutine Bcast_INS_OperatorOptions_3D(this, root, comm)
    class(INS_OperatorOptions_3D), intent(inout) :: this
    integer,        intent(in) :: root !< rank of broadcast root
    type(MPI_Comm), intent(in) :: comm !< MPI communicator

    call XMPI_Bcast(this % convection_term , root, comm)
    call XMPI_Bcast(this % pressure_solver , root, comm)
    call XMPI_Bcast(this % diffusion_solver, root, comm)
    call XMPI_Bcast(this % interior_bc     , root, comm)
    call XMPI_Bcast(this % dealiasing      , root, comm)
    call XMPI_Bcast(this % penalty_u       , root, comm)
    call XMPI_Bcast(this % penalty_p       , root, comm)
    call XMPI_Bcast(this % div_stab        , root, comm)
    call XMPI_Bcast(this % mu_0            , root, comm)
    call XMPI_Bcast(this % delta_out       , root, comm)
    call XMPI_Bcast(this % i_max_p         , root, comm)
    call XMPI_Bcast(this % i_max_v         , root, comm)
    call XMPI_Bcast(this % k_max           , root, comm)
    call XMPI_Bcast(this % k_pre_p         , root, comm)
    call XMPI_Bcast(this % k_pre_v         , root, comm)
    call XMPI_Bcast(this % r_red           , root, comm)
    call XMPI_Bcast(this % r_max           , root, comm)

    call this % schwarz_p % Bcast(root, comm)
    call this % schwarz_u % Bcast(root, comm)
    call this % sgs_model % Bcast(root, comm)

  end subroutine Bcast_INS_OperatorOptions_3D

  !=============================================================================

end module INS__Operator__3D
