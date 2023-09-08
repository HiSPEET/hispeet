!> summary:  MultiLevel Spectral Deferred Correction of 1D conservation laws
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
  use Fine_To_Coarse_Interpolation__1D
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

    logical :: is_root      !< T if root (bottom) level mesh
    logical :: is_top       !< T if top level mesh
    integer :: n_step       !< number of time steps in one space-time slice
    integer :: x_refinement !< spatial refinement:  {0,1,2} = {none,p,h}
    integer :: t_refinement !< temporal refinement: {0,1,2} = {none,p,h}

    ! discretization and solvers ...............................................

    class(CL_Problem_1D),    pointer     :: cl_problem   !< same for all levels
    class(CL_Operator_1D),   allocatable :: cl_operator
    class(CL_SDC_Method_1D), allocatable :: cl_sdc


    ! interpolation operators ..................................................

    type(CoarseToFineInterpolation_1D) :: iop_cf_x !< C-F interpolation in x
    type(CoarseToFineInterpolation_1D) :: iop_cf_t !< C-F interpolation in t

    type(FineToCoarseInterpolation_1D) :: iop_fc_x !< F-C interpolation in x
    type(FineToCoarseInterpolation_1D) :: iop_fc_t !< F-C interpolation in t

    ! data .....................................................................

    real(RNP), allocatable :: u(:,:,:,:,:)
    !< approximate solution dimensioned as `u(0:po,1:ne,1:nc,0:ns,1:nt)`, where
    !!  - 'po         '  is the polynomial order (degree) in space
    !!  - 'ne = n_elem'  is the number of elements in space
    !!  - 'nc = n_comp'  is the number of solution components
    !!  - 'ns = n_sub '  is the number of subintervals in each time step
    !!  - 'nt = n_step'  is the number of time steps, as defined above

    real(RNP), allocatable :: v(:,:,:,:,:) !< restricted solution or correction
    real(RNP), allocatable :: g(:,:,:,:,:) !< FAS part of right hand side

  contains

    procedure :: GetResidual    !< compute multi-step collocation residual
    procedure :: ApplyPredictor !< application of SDC predictor
    procedure :: ApplyCorrector !< application of SDC corrector
    procedure :: Interpolate_CF !< solution interpolation to next finer   level
    procedure :: Interpolate_FC !< solution interpolation to next coarser level
    procedure :: Restrict_FC    !< residual restriction from next finer   level

  end type CL_MLSDC_Level_1D

  ! constructor
  interface CL_MLSDC_Level_1D
    module procedure New_CL_MLSDC_Level_1D
  end interface

  !-----------------------------------------------------------------------------
  !> Options for CL_MLSDC_Level_1D initialization

  type CL_MLSDC_Level_Options_1D

    logical :: is_root = .true. !< T if root (bottom) level mesh
    logical :: is_top  = .true. !< T if top level mesh
    integer :: n_step  = 1      !< number of time steps in one space-time slice
    integer :: x_refinement = 0 !< spatial refinement:  {0,1,2} = {none,p,h}
    integer :: t_refinement = 0 !< temporal refinement: {0,1,2} = {none,p,h}

    type(CL_Operator_Options_1D) :: cl_operator

    class(CL_TimeIntegrator_Options_1D), allocatable :: cl_predictor
    class(CL_SDC_Options_1D),            allocatable :: cl_sdc

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
    !> Fine-to-coarse space-time interpolation of solution-like variables

    module subroutine Interpolate_FC(this, u_f, u_c)
      class(CL_MLSDC_Level_1D), intent(in) :: this !< coarse level
      real(RNP), intent(in)    :: u_f(0:,:,:,0:,:) !< fine solution variable
      real(RNP), intent(inout) :: u_c(0:,:,:,0:,:) !< coarse solution variable
    end subroutine Interpolate_FC

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

    if (.not. allocated(opt % cl_predictor)) then
      call Error( 'Init_CL_MLSDC_Level_1D', 'opt % cl_predictor not allocated')
    end if

    if (.not. allocated(opt % cl_sdc)) then
      call Error( 'Init_CL_MLSDC_Level_1D', 'opt % cl_sdc not allocated')
    end if

    if (allocated(this % u)) deallocate(this % u)
    if (allocated(this % v)) deallocate(this % v)
    if (allocated(this % g)) deallocate(this % g)

    ! parameters ...............................................................

    this % is_root       =  opt % is_root
    this % is_top        =  opt % is_top
    this % n_step        =  opt % n_step
    this % x_refinement  =  opt % x_refinement
    this % t_refinement  =  opt % t_refinement

    ! problem and operators ....................................................

    this % cl_problem    => cl_problem
    this % cl_operator   =  CL_Operator_1D(opt % cl_operator)

    select type(sdc_opt => opt % cl_sdc)
    type is(CL_SDC_Options_ISD1_1D)
      this % cl_sdc = CL_SDC_Method_ISD1_1D(opt % cl_predictor, sdc_opt)
    end select

  end subroutine Init_CL_MLSDC_Level_1D

  !=============================================================================

end module CL__MLSDC__Level__1D
