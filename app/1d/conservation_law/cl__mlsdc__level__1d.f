!> summary:  MLSDC level of 1D conservation laws
!> author:   Joerg Stiller
!> date:     2023/08/22
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module CL__MLSDC__Level__1D
  use Kind_Parameters
  use Constants
  use Array_Assignments
  use Execution_Control
  use Coarse_To_Fine_Interpolation__1D
  use Fine_To_Coarse_Projection__1D
  use CL__Problem__1D
  use CL__Operator__1D
  use CL__Time_Integrator__1D
  use CL__SDC__Method__1D
  use CL__SDC__Method__ISD1__1D

  implicit none
  private

  !-----------------------------------------------------------------------------
  !> Type for keeping one level of the MLSDC data structure

  type, public :: CL_MLSDC_Level_1D

    logical :: is_root  !< T if root (bottom) level mesh
    logical :: is_top   !< T if top level mesh

    integer :: n_space  !< number of elements in space
    integer :: p_space  !< polynomial degree in space
    integer :: n_time   !< number of steps in time
    integer :: p_time   !< polynomial degree in time

    ! discretization and solvers ...............................................

    class(CL_Problem_1D),    pointer     :: cl_problem   !< same for all levels
    class(CL_Operator_1D),   allocatable :: cl_operator
    class(CL_SDC_Method_1D), allocatable :: cl_sdc

    ! transfer operators .......................................................

    type(CoarseToFineInterpolation_1D) :: iop_cf_x !< C-F interpolation in x
    type(CoarseToFineInterpolation_1D) :: iop_cf_t !< C-F interpolation in t

    type(FineToCoarseProjection_1D)    :: pop_fc_x !< F-C projection in x
    type(FineToCoarseProjection_1D)    :: pop_fc_t !< F-C projection in t

  contains

    procedure :: GetResidual    !< compute multi-step collocation residual
    procedure :: ApplyPredictor !< application of SDC predictor
    procedure :: ApplyCorrector !< application of SDC corrector
    procedure :: Interpolate_CF !< solution interpolation to next finer   level
    procedure :: Project_FC     !< solution interpolation to next coarser level
    procedure :: Restrict_FC    !< residual restriction from next finer   level

  end type CL_MLSDC_Level_1D

  ! constructor
  interface CL_MLSDC_Level_1D
    module procedure New_CL_MLSDC_Level_1D
  end interface

  !-----------------------------------------------------------------------------
  !> Options for CL_MLSDC_Level_1D initialization

  type, public :: CL_MLSDC_Level_Options_1D

    integer :: n_space = -1 !< number of elements in space
    integer :: p_space = -1 !< polynomial degree in space
    integer :: n_time  = -1 !< number of steps in time
    integer :: p_time  = -1 !< polynomial degree in time

    type(CL_Operator_Options_1D) :: cl_operator

    class(CL_TimeIntegrator_Options_1D), allocatable :: cl_pre
    class(CL_SDC_Options_1D),            allocatable :: cl_sdc

    type(CoarseToFineInterpolationOptions_1D) :: iop_cf_x
    type(CoarseToFineInterpolationOptions_1D) :: iop_cf_t

    type(FineToCoarseProjectionOptions_1D) :: pop_fc_x
    type(FineToCoarseProjectionOptions_1D) :: pop_fc_t

  end type CL_MLSDC_Level_Options_1D

  !=============================================================================
  ! module procedures

  interface

    !---------------------------------------------------------------------------
    !> Application of the SDC predictor

    module subroutine ApplyPredictor(this, dt, t_0, u)
      class(CL_MLSDC_Level_1D), intent(in) :: this
      real(RNP), intent(in)    :: dt             !< size of the time slice
      real(RNP), intent(in)    :: t_0            !< start time of the slice
      real(RNP), intent(inout) :: u(0:,:,:,0:,:) !< approximate solution
    end subroutine ApplyPredictor

    !---------------------------------------------------------------------------
    !> Application of the SDC corrector

    module subroutine ApplyCorrector(this, dt, t_0, G, u, n_sweep)
      class(CL_MLSDC_Level_1D), intent(in) :: this
      real(RNP), intent(in)    :: dt             !< size of the time slice
      real(RNP), intent(in)    :: t_0            !< start time of the slice
      real(RNP), intent(in)    :: G(0:,:,:,0:,:) !< FAS defect correction
      real(RNP), intent(inout) :: u(0:,:,:,0:,:) !< approximate solution
      integer,   intent(in)    :: n_sweep        !< number of sweeps
    end subroutine ApplyCorrector

    !---------------------------------------------------------------------------
    !> Residual of the collocation method

    module subroutine GetResidual(this, dt, t_0, u, r)
      class(CL_MLSDC_Level_1D), intent(in) :: this !< MLSDC level
      real(RNP), intent(in)  :: dt             !< size of the time slice
      real(RNP), intent(in)  :: t_0            !< start time of the slice
      real(RNP), intent(in)  :: u(0:,:,:,0:,:) !< approximate solution
      real(RNP), intent(out) :: r(0:,:,:,0:,:) !< residual
    end subroutine GetResidual

    !---------------------------------------------------------------------------
    !> Coarse-to-fine space-time interpolation of solution-like variables

    module subroutine Interpolate_CF(this, u_c, u_f, complete)
      class(CL_MLSDC_Level_1D), intent(in) :: this !< coarse level
      real(RNP), intent(in)    :: u_c(0:,:,:,0:,:) !< coarse solution variable
      real(RNP), intent(inout) :: u_f(0:,:,:,0:,:) !< fine solution variable
      logical,   intent(in)    :: complete         !< T/F for all/refined elements
    end subroutine Interpolate_CF

    !---------------------------------------------------------------------------
    !> Fine-to-coarse space-time projection of solution-like variables

    module subroutine Project_FC(this, u_f, u_c)
      class(CL_MLSDC_Level_1D), intent(in) :: this !< coarse level
      real(RNP), intent(in)    :: u_f(0:,:,:,0:,:) !< fine solution variable
      real(RNP), intent(inout) :: u_c(0:,:,:,0:,:) !< coarse solution variable
    end subroutine Project_FC

    !---------------------------------------------------------------------------
    !> Fine-to-coarse restriction of residual-like variables

    module subroutine Restrict_FC(this, r_f, r_c)
      class(CL_MLSDC_Level_1D), intent(in) :: this !< coarse level
      real(RNP), intent(in)    :: r_f(0:,:,:,0:,:) !< fine residual variable
      real(RNP), intent(inout) :: r_c(0:,:,:,0:,:) !< coarse residual variable
    end subroutine Restrict_FC

  end interface

contains

  !-----------------------------------------------------------------------------
  !> Returns a new CL_MLSDC_Level_1D object

  function New_CL_MLSDC_Level_1D(opt, cl_problem) result(this)
    class(CL_MLSDC_Level_Options_1D), intent(in) :: opt
    class(CL_Problem_1D), target,     intent(in) :: cl_problem
    type(CL_MLSDC_Level_1D) :: this

    call Init_CL_MLSDC_Level_1D(this, opt, cl_problem)

  end function New_CL_MLSDC_Level_1D

  !-----------------------------------------------------------------------------
  !> Inititialization of an MLSDC level

  subroutine Init_CL_MLSDC_Level_1D(this, opt, cl_problem)
    class(CL_MLSDC_Level_1D),         intent(inout) :: this
    class(CL_MLSDC_Level_Options_1D), intent(in)    :: opt
    class(CL_Problem_1D), target,     intent(in)    :: cl_problem

    ! safeguard ................................................................

    if (.not. allocated(opt % cl_pre)) then
      call Error( 'Init_CL_MLSDC_Level_1D', 'opt % cl_pre not allocated')
    end if

    if (.not. allocated(opt % cl_sdc)) then
      call Error( 'Init_CL_MLSDC_Level_1D', 'opt % cl_sdc not allocated')
    end if

    ! parameters ...............................................................

    this % is_root = opt % pop_fc_x % po_c <= 0 .or. opt % pop_fc_t % po_c <= 0
    this % is_top  = opt % iop_cf_x % po_f <= 0 .or. opt % iop_cf_t % po_f <= 0

    this % n_space = opt % n_space
    this % p_space = opt % p_space
    this % n_time  = opt % n_time
    this % p_time  = opt % p_time

    ! discretization and solvers ...............................................

    this % cl_problem   => cl_problem
    this % cl_operator  =  CL_Operator_1D(opt % cl_operator)

    select type(sdc_opt => opt % cl_sdc)
    type is(CL_SDC_Options_ISD1_1D)
      this % cl_sdc = CL_SDC_Method_ISD1_1D(opt % cl_pre, sdc_opt)
    end select

    ! transfer operators .......................................................

    if (.not. this % is_top) then
      this % iop_cf_x = CoarseToFineInterpolation_1D(opt % iop_cf_x)
      this % iop_cf_t = CoarseToFineInterpolation_1D(opt % iop_cf_t)
    end if

    if (.not. this % is_root) then
      this % pop_fc_x = FineToCoarseProjection_1D(opt % pop_fc_x)
      this % pop_fc_t = FineToCoarseProjection_1D(opt % pop_fc_t)
    end if

  end subroutine Init_CL_MLSDC_Level_1D

  !=============================================================================

end module CL__MLSDC__Level__1D
