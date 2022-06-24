module Element_Distribution_Data__3D
  use Kind_Parameters
  use XMPI
  use Mesh_Element__3D

  implicit none
  private

  !-----------------------------------------------------------------------------
  !>

  type ElementDistributionData_3D

    ! MPI
    integer            :: proc = -1   !< MPI rank of source/destination
    type(MPI_Comm)     :: comm        !< MPI communicator
    type(MPI_Request)  :: req_dim(3)  !< MPI request for dimensions
    type(MPI_Request)  :: req_dat(6)  !< MPI request for data arrays

    ! dimensions
    integer :: n_element  = -1        !< number of elements
    integer :: n_neighbor = -1        !< number of neighbors
    integer :: n_points   = -1        !< number of element points

    ! static element components
    type(MeshElement_3D), allocatable :: element(:)      !< static components

    ! concatenated element neighbor components
    integer     , allocatable :: neighbor_id(:)          !< neighbor%id
    integer     , allocatable :: neighbor_part(:)        !< neighbor%part
    integer(IXS), allocatable :: neighbor_component(:)   !< neighbor%component
    integer(IXS), allocatable :: neighbor_orientation(:) !< neighbor%orientation

    ! concatenated element points
    real(RNP), allocatable :: x_e(:,:)                   !< geometry%x_e

  contains

    procedure :: Send_Dimensions
    procedure :: Send_Data
    procedure :: Send_Finish
    procedure :: Recv_Dimensions
    procedure :: Recv_Data
    procedure :: Recv_Finish

  end type ElementDistributionData_3D

  !-----------------------------------------------------------------------------
  !> MPI datatype for elements

  type(MPI_Datatype) :: MPI_Element = MPI_DATATYPE_NULL

contains

  !-----------------------------------------------------------------------------
  !>

  !=============================================================================

end module Element_Distribution_Data__3D
