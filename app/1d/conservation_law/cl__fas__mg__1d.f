  !-----------------------------------------------------------------------------
  !>
  !>  Notes
  !>    - In the simplest case the smoother is just `sdc % CorrectorStep`

  type CL_FAS_MG

    logical :: is_root = .true.   !< T if root (bottom) level mesh
    logical :: is_top  = .true.   !< T if top level mesh

    ! discretization and solvers ...............................................

    class(CL_Problem_1D),    pointer     :: cl_problem   !< same for all levels
    class(CL_Operator_1D),   allocatable :: cl_operator
    class(CL_SDC_Method_1D), allocatable :: cl_sdc

    integer :: nt = 1  !< number of time steps

    ! interpolation operators ..................................................

    real(RNP), allocatable :: icf_x(:,:,:)
    !< coarse-to-fine interpolation operator in space
    !!   - `size(icf_x,3) = 0`  no refinement
    !!   - `size(icf_x,3) = 1`  p refinement
    !!   - `size(icf_x,3) = 2`  h or hp refinement

    real(RNP), allocatable  :: ifc_x(:,:,:)
    !< fine-to-coarse interpolation operator in space
    !!   - `size(icf_x,3) = 0`  no coarsening
    !!   - `size(icf_x,3) = 1`  p coarsening
    !!   - `size(icf_x,3) = 2`  h or hp coarsening

    real(RNP), allocatable :: icf_t(:,:,:)
    !< coarse-to-fine interpolation operator in time

    real(RNP), allocatable :: ifc_t(:,:,:)
    !< fine-to-coarse interpolation operator in time

    ! data .....................................................................

    real(RNP), allocatable :: u(:,:,:,:,:)
    !< approximate solution dimensioned as `u(0:po,1:ne,1:nc,0:ns,1:nt)`, where
    !!  - 'po'  is the polynomial order (degree) in space
    !!  - 'ne'  is the number of elements in space
    !!  - 'nc'  is the number of solution components
    !!  - 'ns'  is the number of subintervals in each time step (degree in time)
    !!  - 'nt'  is the number of time steps

    ! possibly some more/auxiliary data

  contains

    procedure :: Interpolate_F !< space-time interpolation to next finer level
    procedure :: Interpolate_C !< space-time interpolation to next coarser level
    procedure :: Restrict_C    !< restriction to next coarser level

  end type CL_FAS_MG
