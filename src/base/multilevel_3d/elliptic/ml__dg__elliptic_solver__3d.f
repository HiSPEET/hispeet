!> summary:  3D FAS multigrid DG solver for elliptic problems
!> author:   Joerg Stiller
!> date:     2024/08/20
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module ML__DG__Elliptic_Solver__3D
  use Kind_Parameters
  use DG__Element_Operators__1D
  use DG__Elliptic_Operator__3D
  use DG__Schwarz_Operator__3D
  use ML__Mesh_Operators__3D

  implicit none
  private

  public :: ML_DG_EllipticSolver_3D
  public :: ML_DG_EllipticOptions_3D

  !-----------------------------------------------------------------------------
  !> Type providing 3D FAS multigrid DG solvers for elliptic problems

  type ML_DG_EllipticSolver_3D

    class(ML_MeshOperators_3D), pointer :: ml_op => null()
      !< multilevel mesh operators

    type(DG_EllipticOperator_3D), allocatable :: elliptic_op(:)
      !< elliptic operators for each level

    character(len=3) :: smooth_method !< smoothing method
    character(len=3) :: coarse_method !< coarse grid solver

    integer, allocatable :: n_s1(:) !< number of pre-smoothing steps
    integer, allocatable :: n_s2(:) !< number of post-smoothing steps

    integer   :: i_crs   !< max number of coarse solver iterations
    integer   :: i_max   !< max number of multigrid iterations (cycles)
    real(RNP) :: r_red   !< min residual reduction
    real(RNP) :: r_max   !< max admissible residual
    real(RNP) :: dr_min  !< termination threshold for Δr

  contains

    procedure :: Init_ML_DG_EllipticSolver_3D
  ! procedure :: MG_Solver !< FAS multigrid solver
  ! procedure :: MK_Solver !< FAS-accelerated multilevel Krylov solver

  end type ML_DG_EllipticSolver_3D

  ! constructor interface
  interface ML_DG_EllipticSolver_3D
    procedure New_ML_DG_EllipticSolver_3D
  end interface

  !-----------------------------------------------------------------------------
  !> 3D FAS multigrid solver options
  !>
  !> Smoothing and coarse grid solution methods
  !>   - `'WAS'`  weighted additive Schwarz
  !>   - `'FCG'`  flexible conjugate gradient method
  !>   - `'SCG'`  Schwarz-preconditioned flexible conjugate gradient method

  type ML_DG_EllipticOptions_3D

    real(RNP)                  :: penalty = 2 !< penalty parameter > 1
    type(DG_SchwarzOptions_3D) :: schwarz     !< Schwarz operator options
    character, allocatable     :: bc(:)       !< boundary conditions

    character(len=3) :: smooth_method = 'WAS' !< smoothing method
    character(len=3) :: coarse_method = 'SCG' !< coarse grid solver

    integer   :: n_s1_top =  1 !< num pre-smoothing  steps on top level
    integer   :: n_s2_top =  1 !< num post-smoothing steps on top level
    integer   :: n_s_mult =  1 !< variable smoothing multiplier
    integer   :: i_crs    =  1 !< max number of coarse solver iterations
    integer   :: i_max    =  1 !< max number of multigrid iterations (cycles)
    real(RNP) :: r_red    = -1 !< min residual reduction
    real(RNP) :: r_max    = -1 !< max admissible residual
    real(RNP) :: dr_min   = -1 !< termination threshold for Δr

  end type ML_DG_EllipticOptions_3D

contains

  !-----------------------------------------------------------------------------
  !> Constructor of  3D FAS multigrid solver

  function New_ML_DG_EllipticSolver_3D(ml_op, opt) result(this)
    class(ML_MeshOperators_3D), target, intent(in)    :: ml_op
    class(ML_DG_EllipticOptions_3D),    intent(in)    :: opt

    type(ML_DG_EllipticSolver_3D) :: this

    call Init_ML_DG_EllipticSolver_3D(this, ml_op, opt)

  end function New_ML_DG_EllipticSolver_3D

  !-----------------------------------------------------------------------------
  !> Initialization of 3D FAS multigrid solver

  subroutine Init_ML_DG_EllipticSolver_3D(this, ml_op, opt)
    class(ML_DG_EllipticSolver_3D),     intent(inout) :: this
    class(ML_MeshOperators_3D), target, intent(in)    :: ml_op
    class(ML_DG_EllipticOptions_3D),    intent(in)    :: opt

    type(DG_EllipticOperator_3D), allocatable :: elliptic_op(:)
    integer, allocatable :: n_s1(:), n_s2(:)

    type(DG_ElementOptions_1D) :: dg_opt
    integer :: l, l_top

    l_top = size(ml_op % sem)

    allocate(elliptic_op(l_top), n_s1(l_top), n_s2(l_top))

    do l = 1, l_top

      dg_opt % po      = ml_op % sem(l) % std_op % po
      dg_opt % basis   = ml_op % sem(l) % std_op % basis
      dg_opt % penalty = opt % penalty

      elliptic_op(l) = DG_EllipticOperator_3D( sem         = ml_op % sem(l) &
                                             , dg_opt      = dg_opt         &
                                             , schwarz_opt = opt % schwarz  &
                                             , bc          = opt % bc       )

      if (l == 1 .or. opt % n_s_mult == 1) then
        n_s1(l) = opt%n_s1_top
        n_s2(l) = opt%n_s2_top
      else
        n_s1(l) = n_s1(l-1) * opt%n_s_mult
        n_s2(l) = n_s2(l-1) * opt%n_s_mult
      end if

    end do

    this % ml_op => ml_op

    call move_alloc( elliptic_op, this % elliptic_op )
    call move_alloc( n_s1       , this % n_s1        )
    call move_alloc( n_s2       , this % n_s2        )

    this % smooth_method = opt % smooth_method
    this % coarse_method = opt % coarse_method

    this % i_crs  = opt % i_crs
    this % i_max  = opt % i_max
    this % r_red  = opt % r_red
    this % r_max  = opt % r_max
    this % dr_min = opt % dr_min

  end subroutine Init_ML_DG_EllipticSolver_3D

  !=============================================================================

end module ML__DG__Elliptic_Solver__3D
