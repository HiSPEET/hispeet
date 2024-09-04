!> summary:  3D FAS multigrid DG solver for elliptic problems
!> author:   Joerg Stiller
!> date:     2024/08/20
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module ML__DG__Elliptic_Solver__3D
  use Kind_Parameters
  use Constants
  use Array_Assignments
  use DG__Element_Operators__1D
  use DG__Elliptic_Operator__3D
  use DG__Schwarz_Operator__3D
  use Parent_To_Child_Interpolation__3D
  use ML__Boundary_Variable__3D
  use ML__Mesh_Operators__3D
  use ML__Mesh_Variable__3D

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

    integer   :: start_method  !< starting method
    integer   :: smooth_method !< smoothing method
    integer   :: coarse_solver !< coarse grid solver

    integer   :: i_crs  !< max number of coarse solver iterations
    integer   :: i_max  !< max number of multigrid iterations (cycles)
    integer   :: ns_1   !< number of pre-smoothing steps
    integer   :: ns_2   !< number of post-smoothing steps
    integer   :: ns_c   !< number of continuation smoothing steps
    real(RNP) :: r_red  !< min residual reduction,  if > 0
    real(RNP) :: r_max  !< max admissible residual, if > 0

  contains

    procedure, public  :: Init_ML_DG_EllipticSolver_3D

    generic,   public  :: MG_Solver => MG_Solver_C, MG_Solver_V
    procedure, private :: MG_Solver_C, MG_Solver_V

  ! generic, public :: MK_Solver => MK_Solver_C, MK_Solver_V

    generic,   private :: CoarseSolver => CoarseSolver_C, CoarseSolver_V
    procedure, private :: CoarseSolver_C, CoarseSolver_V

    generic,   private :: Smoother => Smoother_C, Smoother_V
    procedure, private :: Smoother_C, Smoother_V

  end type ML_DG_EllipticSolver_3D

  ! constructor interface
  interface ML_DG_EllipticSolver_3D
    procedure New_ML_DG_EllipticSolver_3D
  end interface

  !-----------------------------------------------------------------------------
  !> 3D FAS multigrid solver options

  type ML_DG_EllipticOptions_3D

    real(RNP)                  :: penalty = 2 !< penalty parameter > 1
    type(DG_SchwarzOptions_3D) :: schwarz     !< Schwarz operator options
    character, allocatable     :: bc(:)       !< boundary conditions

    integer   :: start_method  = START_FMG    !< starting method
    integer   :: smooth_method = SOLVER_WS    !< smoothing method
    integer   :: coarse_solver = SOLVER_SPCG  !< coarse grid solver

    integer   :: i_crs =  1 !< max number of coarse solver iterations
    integer   :: i_max =  1 !< max number of multigrid iterations (cycles)
    integer   :: ns_1  =  1 !< number of pre-smoothing steps
    integer   :: ns_2  =  1 !< number of post-smoothing steps
    integer   :: ns_c  =  1 !< number of continuation smoothing steps
    real(RNP) :: r_red = -1 !< min residual reduction
    real(RNP) :: r_max = -1 !< max admissible residual

  end type ML_DG_EllipticOptions_3D

  !=============================================================================
  ! External module procedures

  interface

    !---------------------------------------------------------------------------
    !> Coarse grid solver with constant diffusivity

    module subroutine CoarseSolver_C(this, lambda, nu, u, f, bv)
      use Boundary_Variable__3D
      class(ML_DG_EllipticSolver_3D), intent(in) :: this
      real(RNP), intent(in) :: lambda
      real(RNP), intent(in) :: nu
      real(RNP), contiguous, intent(inout) :: u(:,:,:,:)
      real(RNP), contiguous, intent(in) :: f(:,:,:,:)
      class(BoundaryVariable_3D), intent(in) :: bv(:)
    end subroutine CoarseSolver_C

    !---------------------------------------------------------------------------
    !> Coarse grid solver with variable diffusivity

    module subroutine CoarseSolver_V(this, lambda, nu, u, f, bv)
      use Boundary_Variable__3D
      class(ML_DG_EllipticSolver_3D), intent(in) :: this
      real(RNP), intent(in) :: lambda
      real(RNP), contiguous, intent(in) :: nu(:,:,:,:)
      real(RNP), contiguous, intent(inout) :: u(:,:,:,:)
      real(RNP), contiguous, intent(in) :: f(:,:,:,:)
      class(BoundaryVariable_3D), intent(in) :: bv(:)
    end subroutine CoarseSolver_V

    !---------------------------------------------------------------------------
    !> Smoother with constant diffusivity

    module subroutine Smoother_C(this, lambda, nu, u, f, bv, n_s)
      use Boundary_Variable__3D
      class(ML_DG_EllipticSolver_3D), intent(in) :: this
      real(RNP), intent(in) :: lambda
      real(RNP), intent(in) :: nu
      real(RNP), contiguous, intent(inout) :: u(:,:,:,:)
      real(RNP), contiguous, intent(in) :: f(:,:,:,:)
      class(BoundaryVariable_3D), intent(in) :: bv(:)
      integer, intent(in) :: n_s
    end subroutine Smoother_C

    !---------------------------------------------------------------------------
    !> Smoother with variable diffusivity

    module subroutine Smoother_V(this, lambda, nu, u, f, bv, n_s)
      use Boundary_Variable__3D
      class(ML_DG_EllipticSolver_3D), intent(in) :: this
      real(RNP), intent(in) :: lambda
      real(RNP), contiguous, intent(in) :: nu(:,:,:,:)
      real(RNP), contiguous, intent(inout) :: u(:,:,:,:)
      real(RNP), contiguous, intent(in) :: f(:,:,:,:)
      class(BoundaryVariable_3D), intent(in) :: bv(:)
      integer, intent(in) :: n_s
    end subroutine Smoother_V

    !---------------------------------------------------------------------------
    !> FAS-MG solver for problems with constant diffusivity

    module subroutine MG_Solver_C(this, lambda, nu, u, f, bv, n_i, r_2)
      class(ML_DG_EllipticSolver_3D), intent(in) :: this
      real(RNP), intent(in) :: lambda
      real(RNP), intent(in) :: nu
      class(ML_MeshVariable_3D), intent(inout) :: u
      class(ML_MeshVariable_3D), intent(inout) :: f
      class(ML_BoundaryVariable_3D), intent(in) :: bv
      integer, optional, intent(out) :: n_i
      real(RNP), optional, intent(out) :: r_2(:)
    end subroutine MG_Solver_C

    !---------------------------------------------------------------------------
    !> FAS-MG solver for problems with variable diffusivity

    module subroutine MG_Solver_V(this, lambda, nu, u, f, bv, n_i, r_2)
      class(ML_DG_EllipticSolver_3D), intent(in) :: this
      real(RNP), intent(in) :: lambda
      class(ML_MeshVariable_3D), intent(in) :: nu
      class(ML_MeshVariable_3D), intent(inout) :: u
      class(ML_MeshVariable_3D), intent(inout) :: f
      class(ML_BoundaryVariable_3D), intent(in) :: bv
      integer, optional, intent(out) :: n_i
      real(RNP), optional, intent(out) :: r_2(:)
    end subroutine MG_Solver_V

  end interface

contains

  !-----------------------------------------------------------------------------
  !> Constructor of  3D FAS multigrid solver

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

    type(DG_ElementOptions_1D) :: dg_opt
    integer :: l, l_top

    l_top = size(ml_op % sem)

    allocate(this % elliptic_op(l_top))

    do l = 1, l_top

      dg_opt % po      = ml_op % sem(l) % std_op % po
      dg_opt % basis   = ml_op % sem(l) % std_op % basis
      dg_opt % penalty = opt % penalty

      this % elliptic_op(l) = &
                 DG_EllipticOperator_3D( sem         = ml_op % sem(l) &
                                       , dg_opt      = dg_opt         &
                                       , schwarz_opt = opt % schwarz  &
                                       , bc          = opt % bc       )

    end do

    this % ml_op => ml_op

    this % start_method  = opt % start_method
    this % smooth_method = opt % smooth_method
    this % coarse_solver = opt % coarse_solver

    this % i_crs  = opt % i_crs
    this % i_max  = opt % i_max
    this % ns_1   = opt % ns_1
    this % ns_2   = opt % ns_2
    this % ns_c   = opt % ns_c
    this % r_red  = opt % r_red
    this % r_max  = opt % r_max

  end subroutine Init_ML_DG_EllipticSolver_3D

  !=============================================================================

end module ML__DG__Elliptic_Solver__3D
