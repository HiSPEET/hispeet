!> summary:  Incompressible Navier-Stokes DG-SEM operator
!> author:   Joerg Stiller
!> date:     2021/12/30
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> @todo
!>   - revision of boundary variables
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

  implicit none
  private

  public :: INS_Operator_3D
  public :: INS_OperatorOptions_3D

  !-----------------------------------------------------------------------------
  !> DG-SEM mesh and operators for incompressible Navier-Stokes problems

  type INS_Operator_3D

    class(INS_Problem_3D), pointer :: problem => null() !< flow problem

    integer :: level !< rank in multilevel hierarchy (0 if none)
    character(len=4) :: pressure_solver  !< pressure solver
    character(len=4) :: diffusion_solver !< diffusion solver

    real(RNP) :: mu_0      !< bulk viscosity,  μ = ζ/ρ
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

    procedure :: ApplyEssentialBC
    procedure :: ApplyNaturalBC

    procedure :: DiffusionSolver

    procedure :: GetBackflowPenalty
    procedure :: GetConvectionTerm

    procedure :: ApplyDiffusionOperator
    procedure :: ApplyDiffusionOperator_C
    procedure :: ApplyDiffusionOperator_V

    procedure :: GetDiffusionResidual
    procedure :: GetDiffusionResidual_C
    procedure :: GetDiffusionResidual_V

    procedure :: GetDiffusionTerm
    procedure :: GetDiffusionTerm_C
    procedure :: GetDiffusionTerm_V

    procedure :: GetStokesResidual

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

    character(4) :: pressure_solver  = 'SPCG'  !< {'AS','CG','SPCG','MG','MGCG'}
    character(4) :: diffusion_solver = 'DPCG'  !< {'DPCG','SPCG'}
    logical      :: dealiasing       = .false. !< F: no dealiasing, T: 3/2 rule

    real(RNP)    :: penalty_p =    -1 !< penalty for p-solver, -1: auto
    real(RNP)    :: penalty_u =    -1 !< penalty for u-solver, -1: auto
    real(RNP)    :: mu_0      =     0 !< bulk viscosity, μ = ζ/ρ
    real(RNP)    :: delta_out =  0.01 !< outflow parameter
    integer      :: i_max_p   =  1000 !< max num p-iterations in projection
    integer      :: i_max_v   =   200 !< max num v-iterations in projection
    integer      :: k_max     =     0 !< max num Krylov iterations
    integer      :: k_pre_p   =    10 !< max num p-iterations in Krylov precon
    integer      :: k_pre_v   =    10 !< max num v-iterations in Krylov precon
    real(RNP)    :: r_red     = 1e-08 !< min residual reduction, if > 0
    real(RNP)    :: r_max     = 1e-12 !< max residual to reach,  if > 0

    type(DG_SchwarzOptions_3D) :: schwarz_u !< Schwarz options for u-solver
    type(DG_SchwarzOptions_3D) :: schwarz_p !< Schwarz options for p-solver

  contains
    procedure :: Bcast => Bcast_INS_OperatorOptions_3D
  end type INS_OperatorOptions_3D

  !=============================================================================
  ! Module procedures (alphabetically)

  interface

    !---------------------------------------------------------------------------
    !> Homogeneous diffusion operator with constant viscosity

    module subroutine ApplyDiffusionOperator_C(this, tau, v, r, form)
      class(INS_Operator_3D), intent(in)  :: this
      real(RNP),              intent(in)  :: tau
      real(RNP),  contiguous, intent(in)  :: v(:,:,:,:,:)
      real(RNP),  contiguous, intent(out) :: r(:,:,:,:,:)
      integer,      optional, intent(in)  :: form
    end subroutine ApplyDiffusionOperator_C

    !---------------------------------------------------------------------------
    !> Homogeneous diffusion operator with variable viscosity

    module subroutine ApplyDiffusionOperator_V(this, tau, mu, nu, v, r, form)
      class(INS_Operator_3D), intent(in)  :: this
      real(RNP),              intent(in)  :: tau
      real(RNP), contiguous,  intent(in)  :: mu(:,:,:,:)
      real(RNP), contiguous,  intent(in)  :: nu(:,:,:,:)
      real(RNP), contiguous,  intent(in)  :: v(:,:,:,:,:)
      real(RNP), contiguous,  intent(out) :: r(:,:,:,:,:)
      integer,     optional,  intent(in)  :: form
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
      class(INS_Operator_3D),     intent(in)    :: this
      class(BoundaryVariable_3D), intent(in)    :: bv_s(:)
      real(RNP), contiguous,      intent(in)    :: sm(:,:,:,:,:)
      real(RNP), contiguous,      intent(inout) :: sp(:,:,:,:,:)
    end subroutine ApplyNaturalBC

    !---------------------------------------------------------------------------
    !> Diffusion solver

    module subroutine DiffusionSolver(this, tau, mu, nu, f, bv_u, v, precon, ni)
      class(INS_Operator_3D),          intent(in)    :: this
      real(RNP),                       intent(in)    :: tau
      real(RNP), contiguous, optional, intent(in)    :: mu(:,:,:,:)
      real(RNP), contiguous, optional, intent(in)    :: nu(:,:,:,:)
      real(RNP), contiguous,           intent(in)    :: f(:,:,:,:,:)
      class(BoundaryVariable_3D),      intent(in)    :: bv_u(:)
      real(RNP), contiguous,           intent(inout) :: v(:,:,:,:,:)
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

    module subroutine GetDiffusionResidual_C(this, tau, f, bv_u, v, r, form)
      class(INS_Operator_3D),     intent(in)  :: this
      real(RNP),                  intent(in)  :: tau
      real(RNP), contiguous,      intent(in)  :: f(:,:,:,:,:)
      class(BoundaryVariable_3D), intent(in)  :: bv_u(:)
      real(RNP), contiguous,      intent(in)  :: v(:,:,:,:,:)
      real(RNP), contiguous,      intent(out) :: r(:,:,:,:,:)
      integer,          optional, intent(in)  :: form
    end subroutine GetDiffusionResidual_C

    !---------------------------------------------------------------------------
    !> Diffusion residual with variable viscosity

    module subroutine GetDiffusionResidual_V &
        (this, tau, mu, nu, f, bv_u, v, r, form)
      class(INS_Operator_3D),     intent(in)  :: this
      real(RNP),                  intent(in)  :: tau
      real(RNP), contiguous,      intent(in)  :: mu(:,:,:,:)
      real(RNP), contiguous,      intent(in)  :: nu(:,:,:,:)
      real(RNP), contiguous,      intent(in)  :: f(:,:,:,:,:)
      class(BoundaryVariable_3D), intent(in)  :: bv_u(:)
      real(RNP), contiguous,      intent(in)  :: v(:,:,:,:,:)
      real(RNP), contiguous,      intent(out) :: r(:,:,:,:,:)
      integer,          optional, intent(in)  :: form
    end subroutine GetDiffusionResidual_V

    !---------------------------------------------------------------------------
    !> Diffusion term with constant viscosity on irregular (deformed) mesh

    module subroutine GetDiffusionTerm_C(this, v, vp, sp, F_d, bv_u, xout, form)
      class(INS_Operator_3D),               intent(in)  :: this
      real(RNP),                contiguous, intent(in)  :: v(:,:,:,:,:)
      real(RNP),                contiguous, intent(out) :: vp(:,:,:,:,:)
      real(RNP),                contiguous, intent(out) :: sp(:,:,:,:,:)
      real(RNP),                contiguous, intent(out) :: F_d(:,:,:,:,:)
      class(BoundaryVariable_3D), optional, intent(in)  :: bv_u(:)
      logical,                    optional, intent(in)  :: xout
      integer,                    optional, intent(in)  :: form
    end subroutine GetDiffusionTerm_C

    !---------------------------------------------------------------------------
    !> Diffusion term with variable viscosity on irregular (deformed) mesh

    module subroutine GetDiffusionTerm_V &
        (this, mu, nu, v, vp, sp, F_d, bv_u, xout, form)
      class(INS_Operator_3D),               intent(in)  :: this
      real(RNP), contiguous,                intent(in)  :: mu(:,:,:,:)
      real(RNP), contiguous,                intent(in)  :: nu(:,:,:,:)
      real(RNP), contiguous,                intent(in)  :: v(:,:,:,:,:)
      real(RNP), contiguous,                intent(out) :: vp(:,:,:,:,:)
      real(RNP), contiguous,                intent(out) :: sp(:,:,:,:,:)
      real(RNP), contiguous,                intent(out) :: F_d(:,:,:,:,:)
      class(BoundaryVariable_3D), optional, intent(in)  :: bv_u(:)
      logical,                    optional, intent(in)  :: xout
      integer,                    optional, intent(in)  :: form
    end subroutine GetDiffusionTerm_V

    !---------------------------------------------------------------------------
    !> Stokes residual for incompressible flow

    module subroutine GetStokesResidual(this, tau, f, bv_u, mu, nu, u, r)
      class(INS_Operator_3D),          intent(in)    :: this
      real(RNP),                       intent(in)    :: tau
      real(RNP), contiguous,           intent(in)    :: f(:,:,:,:,:)
      class(BoundaryVariable_3D),      intent(in)    :: bv_u(:)
      real(RNP), contiguous, optional, intent(inout) :: mu(:,:,:,:)
      real(RNP), contiguous, optional, intent(inout) :: nu(:,:,:,:)
      real(RNP), contiguous,           intent(in)    :: u(:,:,:,:,:)
      real(RNP), contiguous,           intent(out)   :: r(:,:,:,:,:)
    end subroutine GetStokesResidual

    !---------------------------------------------------------------------------
    !> Viscous stress vector on a boundary (C)

    module subroutine GetViscousBoundaryStress_C &
        (this, b, v, sb, bv_u, xout, form)
      class(INS_Operator_3D),               intent(in)  :: this
      integer,                              intent(in)  :: b
      real(RNP), contiguous,                intent(in)  :: v(:,:,:,:,:)
      real(RNP), contiguous,                intent(out) :: sb(:,:,:,:)
      class(BoundaryVariable_3D), optional, intent(in)  :: bv_u(:)
      logical,                    optional, intent(in)  :: xout
      integer,                    optional, intent(in)  :: form
    end subroutine GetViscousBoundaryStress_C

    !---------------------------------------------------------------------------
    !> Viscous stress vector on a boundary (V)

    module subroutine GetViscousBoundaryStress_V &
        (this, b, mu, nu, v, sb, bv_u, xout, form)
      class(INS_Operator_3D),               intent(in)  :: this
      integer,                              intent(in)  :: b
      real(RNP), contiguous,                intent(in)  :: mu(:,:,:,:)
      real(RNP), contiguous,                intent(in)  :: nu(:,:,:,:)
      real(RNP), contiguous,                intent(in)  :: v(:,:,:,:,:)
      real(RNP), contiguous,                intent(out) :: sb(:,:,:,:)
      class(BoundaryVariable_3D), optional, intent(in)  :: bv_u(:)
      logical,                    optional, intent(in)  :: xout
      integer,                    optional, intent(in)  :: form
    end subroutine GetViscousBoundaryStress_V

    !---------------------------------------------------------------------------
    !> Projection-based pressure solver

    module subroutine PressureSolver(this, tau, bv_u, v, f, p, precon, ni)
      class(INS_Operator_3D),     intent(in)    :: this
      real(RNP),                  intent(in)    :: tau
      class(BoundaryVariable_3D), intent(inout) :: bv_u(:)
      real(RNP), contiguous,      intent(in)    :: v(:,:,:,:,:)
      real(RNP), contiguous,      intent(in)    :: f(:,:,:,:)
      real(RNP), contiguous,      intent(inout) :: p(:,:,:,:)
      logical,          optional, intent(in)    :: precon
      integer,          optional, intent(out)   :: ni
    end subroutine PressureSolver

    !---------------------------------------------------------------------------
    !> Projection-diffusion step for incompressible flow

    module subroutine StokesProjection &
        (this, tau, t, v_0, F_c, F_d, Q, bv_u, mu, nu, u, precon)
      class(INS_Operator_3D),          intent(in)    :: this
      real(RNP),                       intent(in)    :: tau
      real(RNP),                       intent(in)    :: t
      real(RNP), contiguous,           intent(in)    :: v_0(:,:,:,:,:)
      real(RNP), contiguous,           intent(in)    :: F_c(:,:,:,:,:)
      real(RNP), contiguous,           intent(in)    :: F_d(:,:,:,:,:)
      real(RNP), contiguous,           intent(in)    :: Q(:,:,:,:,:)
      class(BoundaryVariable_3D),      intent(in)    :: bv_u(:)
      real(RNP), contiguous, optional, intent(inout) :: mu(:,:,:,:)
      real(RNP), contiguous, optional, intent(inout) :: nu(:,:,:,:)
      real(RNP), contiguous,           intent(inout) :: u(:,:,:,:,:)
      logical,               optional, intent(in)    :: precon
    end subroutine StokesProjection

    !---------------------------------------------------------------------------
    !> FGMRES for Stokes part with projection-diffusion preconditioner

    module subroutine StokesFGMRES &
        (this, tau, t, v_0, F_c, F_d, Q, bv_u, mu, nu, u)
      class(INS_Operator_3D),          intent(in)    :: this
      real(RNP),                       intent(in)    :: tau
      real(RNP),                       intent(in)    :: t
      real(RNP), contiguous,           intent(in)    :: v_0(:,:,:,:,:)
      real(RNP), contiguous,           intent(in)    :: F_c(:,:,:,:,:)
      real(RNP), contiguous,           intent(in)    :: F_d(:,:,:,:,:)
      real(RNP), contiguous,           intent(in)    :: Q(:,:,:,:,:)
      class(BoundaryVariable_3D),      intent(in)    :: bv_u(:)
      real(RNP), contiguous, optional, intent(inout) :: mu(:,:,:,:)
      real(RNP), contiguous, optional, intent(inout) :: nu(:,:,:,:)
      real(RNP), contiguous,           intent(inout) :: u(:,:,:,:,:)
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

    this % pressure_solver  = opt % pressure_solver
    this % diffusion_solver = opt % diffusion_solver
    this % mu_0             = opt % mu_0
    this % delta_out        = opt % delta_out

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
    this % elliptic_p = DG_EllipticOperator_3D( sem_p           &
                                              , opt % schwarz_p &
                                              , opt % penalty_p )

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

    ! iterative solver settings
    this % i_max_p = opt % i_max_p
    this % i_max_v = opt % i_max_v
    this % k_max   = opt % k_max
    this % k_pre_p = opt % k_pre_p
    this % k_pre_v = opt % k_pre_v

  end subroutine Init_INS_Operator_3D

  !-----------------------------------------------------------------------------
  !> Diffusion term with constant viscosity

  subroutine GetConvectionTerm(this, v, vp, F_c)
    class(INS_Operator_3D), intent(in) :: this
    real(RNP), contiguous, intent(in)  :: v(:,:,:,:,:)   !< velocity
    real(RNP), contiguous, intent(in)  :: vp(:,:,:,:,:)  !< outer velocity v⁺
    real(RNP), contiguous, intent(out) :: F_c(:,:,:,:,:) !< convection term

!   if (this % mesh % regular) then
!     not implemented yet
!   else
      call TPO_INS_Convection( nv   = this % eop_u % po + 1       &
                             , nq   = this % sop_q % po + 1       &
                             , ne   = this % mesh % n_elem        &
                             , D_v  = this % eop_u  % D           &
                             , I_vq = this % iop_uq % A           &
                             , w_q  = this % sop_q  % w           &
                             , Jd_q = this % sem_q % metrics % Jd &
                             , Ji_q = this % sem_q % metrics % Ji &
                             , a_q  = this % sem_q % metrics % a  &
                             , n_q  = this % sem_q % metrics % n  &
                             , v    = v                           &
                             , vp   = vp                          &
                             , F_c  = F_c                         )
!   end if

  end subroutine GetConvectionTerm

  !-----------------------------------------------------------------------------
  !> Homogeneous diffusion operator with constant or variable viscosity

  subroutine ApplyDiffusionOperator(this, tau, mu, nu, v, r, form)
    class(INS_Operator_3D),          intent(in)  :: this
    real(RNP),                       intent(in)  :: tau
    real(RNP), contiguous, optional, intent(in)  :: mu(:,:,:,:)
    real(RNP), contiguous, optional, intent(in)  :: nu(:,:,:,:)
    real(RNP), contiguous,           intent(in)  :: v(:,:,:,:,:)
    real(RNP), contiguous,           intent(out) :: r(:,:,:,:,:)
    integer,               optional, intent(in)  :: form

    if (present(mu) .and. present(nu)) then
      call this % ApplyDiffusionOperator_V(tau, mu, nu, v, r, form)
    else
      call this % ApplyDiffusionOperator_C(tau, v, r, form)
    end if

  end subroutine ApplyDiffusionOperator

  !-----------------------------------------------------------------------------
  !> Diffusion residual with constant or variable viscosity

  subroutine GetDiffusionResidual(this, tau, mu, nu, f, bv_u, v, r, form)
    class(INS_Operator_3D),          intent(in)  :: this
    real(RNP),                       intent(in)  :: tau
    real(RNP), contiguous, optional, intent(in)  :: mu(:,:,:,:)
    real(RNP), contiguous, optional, intent(in)  :: nu(:,:,:,:)
    real(RNP), contiguous,           intent(in)  :: f(:,:,:,:,:)
    class(BoundaryVariable_3D),      intent(in)  :: bv_u(:)
    real(RNP), contiguous,           intent(in)  :: v(:,:,:,:,:)
    real(RNP), contiguous,           intent(out) :: r(:,:,:,:,:)
    integer,               optional, intent(in)  :: form

    if (present(mu) .and. present(nu)) then
      call this % GetDiffusionResidual_V(tau, mu, nu, f, bv_u, v, r, form)
    else
      call this % GetDiffusionResidual_C(tau, f, bv_u, v, r, form)
    end if

  end subroutine GetDiffusionResidual

  !-----------------------------------------------------------------------------
  !> Diffusion term with constant or variable viscosity

  subroutine GetDiffusionTerm(this, mu, nu, v, vp, sp, F_d, bv_u, xout, form)
    class(INS_Operator_3D),               intent(in)  :: this
    real(RNP), contiguous,      optional, intent(in)  :: mu(:,:,:,:)
    real(RNP), contiguous,      optional, intent(in)  :: nu(:,:,:,:)
    real(RNP), contiguous,                intent(in)  :: v(:,:,:,:,:)
    real(RNP), contiguous,                intent(out) :: vp(:,:,:,:,:)
    real(RNP), contiguous,                intent(out) :: sp(:,:,:,:,:)
    real(RNP), contiguous,                intent(out) :: F_d(:,:,:,:,:)
    class(BoundaryVariable_3D), optional, intent(in)  :: bv_u(:)
    logical,                    optional, intent(in)  :: xout
    integer,                    optional, intent(in)  :: form

    if (present(mu) .and. present(nu)) then
      call this % GetDiffusionTerm_V(mu, nu, v, vp, sp, F_d, bv_u, xout, form)
    else
      call this % GetDiffusionTerm_C(v, vp, sp, F_d, bv_u, xout, form)
    end if

  end subroutine GetDiffusionTerm

  !-----------------------------------------------------------------------------
  !> Viscous stress vector with constant or variable viscosity on a boundary

  subroutine GetViscousBoundaryStress(this, b, mu, nu, v, sb, bv_u, xout, form)
    class(INS_Operator_3D),               intent(in)  :: this
    integer,                              intent(in)  :: b
    real(RNP), contiguous,      optional, intent(in)  :: mu(:,:,:,:)
    real(RNP), contiguous,      optional, intent(in)  :: nu(:,:,:,:)
    real(RNP), contiguous,                intent(in)  :: v(:,:,:,:,:)
    real(RNP), contiguous,                intent(out) :: sb(:,:,:,:)
    class(BoundaryVariable_3D), optional, intent(in)  :: bv_u(:)
    logical,                    optional, intent(in)  :: xout
    integer,                    optional, intent(in)  :: form

    if (present(mu) .and. present(nu)) then
      call this % GetViscousBoundaryStress_V(b, mu, nu, v, sb, bv_u, xout, form)
    else
      call this % GetViscousBoundaryStress_C(b, v, sb, bv_u, xout, form)
    end if

  end subroutine GetViscousBoundaryStress

  !-----------------------------------------------------------------------------
  !>

  subroutine StokesSolver(this, tau, t, v_0, F_c, F_d, Q, bv_u, mu, nu, u)
    class(INS_Operator_3D), intent(in) :: this
    !< incompressible Navier-Stokes time integrator
    real(RNP), intent(in) :: tau
    !< τ, effective time step width
    real(RNP), intent(in) :: t
    !< t, final time
    real(RNP), contiguous, intent(in) :: v_0(:,:,:,:,:)
    !< v₀, effective initial value of velocity
    real(RNP), contiguous, intent(in) :: F_c(:,:,:,:,:)
    !< convective term at final time t
    real(RNP), contiguous, intent(in) :: F_d(:,:,:,:,:)
    !< diffusion term at final time t
    real(RNP), contiguous, intent(in) :: Q(:,:,:,:,:)
    !< sources at time t and further known terms
    class(BoundaryVariable_3D), intent(in) :: bv_u(:)
    !< boundary values at final time t
    !!   - Γᴰ :  [ v₁, v₂, v₃, - , -  ]
    !!   - Γᴼ :  [ - , - , - , p , ∆p ]
    real(RNP), contiguous, optional, intent(inout) :: mu(:,:,:,:)
    !< μ, kinematic bulk viscosity
    real(RNP), contiguous, optional, intent(inout) :: nu(:,:,:,:)
    !< ν, kinematic shear viscosity
    real(RNP), contiguous, intent(inout) :: u(:,:,:,:,:)
    !< u = [v, p], velocity and pressure at final time u

    if (this % k_max > 0) then
      call StokesFGMRES(this, tau, t, v_0, F_c, F_d, Q, bv_u, mu, nu, u)
    else
      call StokesProjection(this, tau, t, v_0, F_c, F_d, Q, bv_u, mu, nu, u)
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

    call XMPI_Bcast(this % pressure_solver , root, comm)
    call XMPI_Bcast(this % diffusion_solver, root, comm)
    call XMPI_Bcast(this % dealiasing      , root, comm)
    call XMPI_Bcast(this % penalty_p       , root, comm)
    call XMPI_Bcast(this % penalty_u       , root, comm)
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

  end subroutine Bcast_INS_OperatorOptions_3D

  !=============================================================================

end module INS__Operator__3D
