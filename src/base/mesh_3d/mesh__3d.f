!> summary:  3D mesh partition type
!> author:   Joerg Stiller
!> date:     2020/11/13
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Mesh__3D
  use Kind_Parameters, only: IXL, IXS, RNP
  use XMPI
  use Mesh_Face__3D
  use Mesh_Element__3D
  use Mesh_Boundary__3D
  use Mesh_Link__3D
  use Mesh_Element_Indexing__3D
  implicit none
  private

  public :: Mesh_3D
  public :: MeshAttributes_3D

  !-----------------------------------------------------------------------------
  !> 3D mesh partition type
  !>
  !> ### Ghost elements
  !>
  !> The ghost elements stored in `ghost(1:n_ghost)` ...

  type Mesh_3D

    ! attributes ...............................................................

    integer :: n_bound = 0           !< number of domain boundaries
    integer :: n_parts = 0           !< number of non-empty partitions

    logical :: structured = .false.  !< T if mapping to structured mesh exists
    logical :: regular    = .false.  !< T if equidistant Cartesian

    logical :: is_root    = .true.   !< T if root (bottom) level mesh
    logical :: is_top     = .true.   !< T if top level mesh

    ! MPI ......................................................................

    type(MPI_Comm) :: comm_world     !< communicator covering all processes
    type(MPI_Comm) :: comm_parts     !< communicator covering active partitions

    integer :: proc = -1             !< process    =  ID in comm_world ≥ 0
    integer :: part = -1             !< partition  =  ID in comm_parts ≥ 0

    integer, allocatable :: proc_part(:) !< map from comm_parts to comm_world

    ! dimensions ...............................................................

    integer :: n_vert  = 0           !< number of mesh vertices
    integer :: n_edge  = 0           !< number of mesh edges
    integer :: n_face  = 0           !< number of mesh faces
    integer :: n_elem  = 0           !< number of mesh elements
    integer :: n_ghost = 0           !< number of ghost elements
    integer :: n_link  = 0           !< number of mesh links
    integer :: p_geom  = 0           !< max polynomial order of element geometry

    integer :: max_vert_val = 0      !< maximum vertex valency
    integer :: max_edge_val = 0      !< maximum edge valency

    ! structured mesh properties ...............................................

    integer :: n_elem_1 = 0          !< number of elements in direction 1
    integer :: n_elem_2 = 0          !< number of elements in direction 2
    integer :: n_elem_3 = 0          !< number of elements in direction 3

    ! regular mesh properties ..................................................

    integer   :: n_face_1 =  0       !< number of faces normal to direction 1
    integer   :: n_face_2 =  0       !< number of faces normal to direction 2
    integer   :: n_face_3 =  0       !< number of faces normal to direction 3
    real(RNP) :: dx(3)    = -1       !< mesh element spacing in directions 1:3

    ! mesh components and links ................................................

    type(MeshFace_3D)    , allocatable :: face(:)     !< mesh faces
    type(MeshElement_3D) , allocatable :: element(:)  !< mesh elements
    type(MeshElement_3D) , allocatable :: ghost(:)    !< ghost elements
    type(MeshBoundary_3D), allocatable :: boundary(:) !< mesh boundaries
    type(MeshLink_3D)    , allocatable :: link(:)     !< mesh links

  contains

    procedure :: Init_Mesh_3D

    procedure :: GetCuboids
    procedure :: GetPoints
    procedure :: GetPointValency

    ! automatic identification and generation of components
    procedure :: BuildCommunicator
    procedure :: BuildFaces
    procedure :: BuildGhosts
    procedure :: BuildLinks
    procedure :: BuildCuboids
    procedure :: IdentifyEdges
    procedure :: IdentifyRanks

    ! import/export
    procedure :: ImportGenericMesh

  end type Mesh_3D

  ! constructor
  interface Mesh_3D
    module procedure EmptyMeshPartition
  end interface

  interface

    !===========================================================================
    ! Service routines

    !---------------------------------------------------------------------------
    !> Returns cuboid points to all elements of a partition

    module subroutine GetCuboids(mesh, x)
      class(Mesh_3D),         intent(in)  :: mesh         !< mesh parition
      real(RNP), allocatable, intent(out) :: x(:,:,:,:,:) !< mesh points
    end subroutine GetCuboids

    !---------------------------------------------------------------------------
    !> Generates Gauss-Lobatto or Gauss points to all elements of a partition

    module subroutine GetPoints(mesh, po, basis, x)
      class(Mesh_3D), intent(in)  :: mesh  !< mesh parition
      integer,                 intent(in)  :: po    !< polynomial order
      character,     optional, intent(in)  :: basis !< 'G' or 'L' ['L']
      real(RNP),  allocatable, intent(out) :: x(:,:,:,:,:) !< mesh points
    end subroutine GetPoints

    !---------------------------------------------------------------------------
    !> Computes the valency of mesh points

    module subroutine GetPointValency(mesh, v)
      class(Mesh_3D), intent(in)  :: mesh          !< mesh parition
      integer,        intent(out) :: v(0:,0:,0:,:) !< point valency
    end subroutine GetPointValency

    !===========================================================================
    ! Routines for constructing mesh components

    !---------------------------------------------------------------------------
    !> Build communicator between active partitions

    module subroutine BuildCommunicator(mesh)
      class(Mesh_3D), intent(inout) :: mesh !< mesh partition
    end subroutine BuildCommunicator

    !---------------------------------------------------------------------------
    !> Generation of approximate cuboids

    module subroutine BuildCuboids(mesh)
      class(Mesh_3D), intent(inout)  :: mesh !< mesh partition
    end subroutine BuildCuboids

    !---------------------------------------------------------------------------
    !> Generation of ghost elements

    module subroutine BuildGhosts(mesh)
      class(Mesh_3D), intent(inout) :: mesh !< local partition
    end subroutine BuildGhosts

    !---------------------------------------------------------------------------
    !> Generation of mesh faces

    module subroutine BuildFaces(mesh)
      class(Mesh_3D), intent(inout) :: mesh !< mesh partition
    end subroutine BuildFaces

    !---------------------------------------------------------------------------
    !> Generation of mesh links from global element neighbor information

    module subroutine BuildLinks(mesh)
      class(Mesh_3D), intent(inout) :: mesh !< local partition
    end subroutine BuildLinks

    !---------------------------------------------------------------------------
    !> Identification of mesh edges

    module subroutine IdentifyEdges(mesh)
      class(Mesh_3D), intent(inout) :: mesh !< mesh partition
    end subroutine IdentifyEdges

    !---------------------------------------------------------------------------
    !> Identification of primary element components

    module subroutine IdentifyRanks(mesh)
      class(Mesh_3D), intent(inout) :: mesh !< mesh partition
    end subroutine IdentifyRanks

    !---------------------------------------------------------------------------
    !> Import a generic 3D mesh

    module subroutine ImportGenericMesh(mesh, generic_mesh, comm)
      use Generic_Mesh__3D
      class(Mesh_3D),        intent(out) :: mesh         !< mesh partition
      class(GenericMesh_3D), intent(in)  :: generic_mesh !< generic mesh
      type(MPI_Comm),        intent(in)  :: comm         !< "world" communicator
    end subroutine ImportGenericMesh

  end interface

  !-----------------------------------------------------------------------------
  !> Type for collecting and transmitting the mesh attributes

  type MeshAttributes_3D

    integer :: n_bound    = 0        !< number of domain boundaries
    integer :: n_parts    = 0        !< number of non-empty partitions
    logical :: structured = .false.  !< T if mapping to structured mesh exists
    logical :: regular    = .false.  !< T if equidistant Cartesian
    logical :: is_root    = .true.   !< T if root (bottom) level mesh
    logical :: is_top     = .true.   !< T if top level mesh

    type(MeshBoundaryAttributes_3D), allocatable :: boundary(:)

  contains
    procedure :: Bcast => Bcast_MeshAttributes_3D
  end type MeshAttributes_3D

  ! constructor
  interface MeshAttributes_3D
    module procedure ExtractMeshAttributes
  end interface

contains

  !=============================================================================
  ! Mesh_3D constructor and type-bound procedures

  !-----------------------------------------------------------------------------
  !> Creates an empty 3D mesh partition

  function EmptyMeshPartition(attrib, comm) result(this)
    class(MeshAttributes_3D), intent(in) :: attrib !< mesh attributes
    type(MPI_Comm),           intent(in) :: comm   !< "world" communicator
    type(Mesh_3D) :: this

    call Init_Mesh_3D(this, attrib, comm)

  end function EmptyMeshPartition

  !-----------------------------------------------------------------------------
  !> Basic initialization of a 3D mesh partition

  subroutine Init_Mesh_3D(this, attrib, comm)
    class(Mesh_3D),           intent(inout) :: this
    class(MeshAttributes_3D), intent(in)    :: attrib !< mesh attributes
    type(MPI_Comm),           intent(in)    :: comm   !< "world" communicator

    integer :: b

    this % n_bound    = attrib % n_bound
    this % n_parts    = attrib % n_parts
    this % structured = attrib % structured
    this % regular    = attrib % regular
    this % is_root    = attrib % is_root
    this % is_top     = attrib % is_top

    allocate(this % boundary( this%n_bound ))
    do b = 1, this % n_bound
      this % boundary(b) = MeshBoundary_3D(attrib % boundary(b))
    end do

    this % comm_world = comm

    call MPI_Comm_rank(comm, this % proc)

  end subroutine Init_Mesh_3D

  !=============================================================================
  ! MeshAttributes_3D constructor and type-bound procedures

  !-----------------------------------------------------------------------------
  !> Creates attributes by extraction from existing mesh partition

  function ExtractMeshAttributes(mesh) result(attrib)
    class(Mesh_3D), intent(in) :: mesh
    type(MeshAttributes_3D) :: attrib

    attrib % n_bound    = mesh % n_bound
    attrib % n_parts    = mesh % n_parts
    attrib % structured = mesh % structured
    attrib % regular    = mesh % regular
    attrib % is_root    = mesh % is_root
    attrib % is_top     = mesh % is_top

    attrib % boundary = MeshBoundaryAttributes_3D(mesh % boundary)

  end function ExtractMeshAttributes

  !-----------------------------------------------------------------------------
  !> Broadcast mesh attributes

  subroutine Bcast_MeshAttributes_3D(this, root, comm)
    class(MeshAttributes_3D), intent(inout) :: this
    integer       , intent(in) :: root !< MPI root process
    type(MPI_Comm), intent(in) :: comm !< MPI communicator

    integer :: attrib_int(2), b
    logical :: attrib_log(4)

    ! mesh attributes ..........................................................

    attrib_int(1) = this % n_bound
    attrib_int(2) = this % n_parts

    attrib_log(1) = this % structured
    attrib_log(2) = this % regular
    attrib_log(3) = this % is_root
    attrib_log(4) = this % is_top

    call XMPI_Bcast(attrib_int, root, comm)
    call XMPI_Bcast(attrib_log, root, comm)

    this % n_bound    = attrib_int(1)
    this % n_parts    = attrib_int(2)

    this % structured = attrib_log(1)
    this % regular    = attrib_log(2)
    this % is_root    = attrib_log(3)
    this % is_top     = attrib_log(4)

    ! boundary attributes ......................................................

    if (allocated(this % boundary)) then
      if (size(this % boundary) /= this % n_bound) then
        deallocate(this % boundary)
      end if
    end if

    if (.not. allocated(this % boundary)) then
      allocate(this % boundary( this%n_bound ))
    end if

    do b = 1, this % n_bound
      call this % boundary(b) % Bcast(root, comm)
    end do

  end subroutine Bcast_MeshAttributes_3D

  !=============================================================================

end module Mesh__3D
