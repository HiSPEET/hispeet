!> summary:  MultiLevel Spectral Deferred Correction of 1D conservation laws
!> author:   Joerg Stiller
!> date:     2023/08/22
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module CL__MLSDC__Level__1D
  use Kind_Parameters
  use Constants
  use Execution_Control
  use Coarse_To_Fine_Interpolation__1D
  use Fine_To_Coarse_Interpolation__1D
  use CL__Problem__1D
  use CL__Operator__1D
  use CL__SDC__Method__1D

  implicit none
  private

  !-----------------------------------------------------------------------------
  !> Type for keeping one level of the MLSDC data structure

  type, public :: CL_MLSDC_Level_1D

    logical :: is_root = .true. !< T if root (bottom) level mesh
    logical :: is_top  = .true. !< T if top level mesh
    integer :: x_refinement = 0 !< spatial refinement:  {0,1,2} = {none,p,h}
    integer :: t_refinement = 0 !< temporal refinement: {0,1,2} = {none,p,h}

    ! discretization and solvers ...............................................

    class(CL_Problem_1D),    pointer     :: cl_problem   !< same for all levels
    class(CL_Operator_1D),   allocatable :: cl_operator
    class(CL_SDC_Method_1D), allocatable :: cl_sdc

    integer :: nt = 1  !< number of time steps in one space-time slice

    ! interpolation operators ..................................................

    type(CoarseToFineInterpolation_1D) :: iop_cf_x !< C-F interpolation in x
    type(CoarseToFineInterpolation_1D) :: iop_cf_t !< C-F interpolation in t

    type(FineToCoarseInterpolation_1D) :: iop_fc_x !< F-C interpolation in x
    type(FineToCoarseInterpolation_1D) :: iop_fc_t !< F-C interpolation in t

    ! data .....................................................................

    real(RNP), allocatable :: u(:,:,:,:,:)
    !< approximate solution dimensioned as `u(0:po,1:ne,1:nc,0:ns,1:nt)`, where
    !!  - 'po'  is the polynomial order (degree) in space
    !!  - 'ne'  is the number of elements in space
    !!  - 'nc'  is the number of solution components
    !!  - 'ns'  is the number of subintervals in each time step (degree in time)
    !!  - 'nt'  is the number of time steps, as defined above

    real(RNP), allocatable :: g(:,:,:,:,:) !< FAS part of right hand side
    real(RNP), allocatable :: v(:,:,:,:,:) !< restricted solution or correction

  contains

    procedure :: GetResidual    !< compute multi-step collocation residual
!   procedure :: ApplyPredictor !< application of SDC predictor
!   procedure :: ApplyCorrector !< application of SDC corrector
    procedure :: Interpolate_CF !< solution interpolation to next finer   level
    procedure :: Interpolate_FC !< solution interpolation to next coarser level
!   procedure :: Restrict_FC    !< residual restriction   to next coarser level

  end type CL_MLSDC_Level_1D

  !=============================================================================
  ! module procedures

  interface

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

    !-----------------------------------------------------------------------------
    !> Fine-to-coarse space-time interpolation of solution-like variables

    module subroutine Interpolate_FC(this, u_f, u_c)
      class(CL_MLSDC_Level_1D), intent(in) :: this !< coarse level
      real(RNP), intent(in)    :: u_f(0:,:,:,0:,:) !< coarse solution variable
      real(RNP), intent(inout) :: u_c(0:,:,:,0:,:) !< fine solution variable
    end subroutine

  end interface

  !=============================================================================

end module CL__MLSDC__Level__1D
