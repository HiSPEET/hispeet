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
  use Mesh_Element_Indexing__3D
  use Mesh_Boundary__3D
  use Mesh_Link__3D
  use Mesh_Map_To_Child__3D
  use Mesh_Map_To_Parent__3D
  implicit none
  private

  public :: Mesh_3D
  public :: MeshAttributes_3D

  public :: Get_H5T_MeshAttributes_3D

  !-----------------------------------------------------------------------------
  !> 3D mesh partition type
  !>
  !> ### Frozen elements
  !>
  !> Frozen elements are elements whose values are defined by interpolation from
  !> the next lower grid level. To easily exclude them from certain operations,
  !> they are paced at the end of the element list. In most other aspects they
  !> are treated as ordinary mesh elements.
  !>
  !> The main purpose of frozen elements is to enclose refinement zones and to
  !> provide boundary conditions for the latter. Due to their external position,
  !> they may have no neighbors at certain faces, edges or vertices.
  !>
  !> ### Ghost elements
  !>
  !> The ghost elements stored in `ghost(1:n_ghost)` ...

  type Mesh_3D

    ! attributes ...............................................................

    integer   :: n_bound    = 0        !< number of domain boundaries
    integer   :: n_parts    = 0        !< number of non-empty partitions

    logical   :: structured = .false.  !< T if mapping to structured mesh exists
    logical   :: regular    = .false.  !< T if equidistant Cartesian
    real(RNP) :: dx(3)      =  0       !< regular mesh spacing in directions 1:3

    logical   :: is_root    = .true.   !< T if root (bottom) level mesh
    logical   :: is_top     = .true.   !< T if top level mesh
    character :: refinement = ''       !< 'c' clone or 's' subdivision

    ! MPI ......................................................................

    type(MPI_Comm) :: comm_world = MPI_COMM_NULL !< global communicator
    type(MPI_Comm) :: comm_parts = MPI_COMM_NULL !< active parts communicator

    integer :: proc = -1 !< process    =  ID in comm_world ≥ 0
    integer :: part = -1 !< partition  =  ID in comm_parts ≥ 0

    integer, allocatable :: proc_part(:) !< map from partition to process IDs

    ! dimensions ...............................................................

    integer :: n_vert        = 0  !< number of mesh vertices
    integer :: n_edge        = 0  !< number of mesh edges
    integer :: n_face        = 0  !< number of mesh faces
    integer :: n_elem        = 0  !< number of mesh elements
    integer :: n_elem_active = 0  !< number of active elements
    integer :: n_elem_frozen = 0  !< number of frozen elements
    integer :: n_cluster     = 0  !< number of sibling element clusters
    integer :: n_ghost       = 0  !< number of ghost elements
    integer :: n_link        = 0  !< number of mesh links
    integer :: n_child       = 0  !< number of child partitions
    integer :: n_parent      = 0  !< number of parent partitions

    integer :: p_geom        = 0  !< max polynomial order of element geometry

    integer :: max_vert_val  = 0  !< maximum vertex valency
    integer :: max_edge_val  = 0  !< maximum edge valency

    ! structured mesh only
    integer :: n_elem_1      = 0  !< number of elements in direction 1
    integer :: n_elem_2      = 0  !< number of elements in direction 2
    integer :: n_elem_3      = 0  !< number of elements in direction 3

    ! mesh components and links ................................................

    type(MeshFace_3D)       , allocatable :: face(:)       !< mesh faces
    type(MeshElement_3D)    , allocatable :: element(:)    !< mesh elements
    type(MeshElement_3D)    , allocatable :: ghost(:)      !< ghost elements
    type(MeshBoundary_3D)   , allocatable :: boundary(:)   !< mesh boundaries
    type(MeshLink_3D)       , allocatable :: link(:)       !< mesh links
    type(MeshMapToChild_3D) , allocatable :: map_child(:)  !< map to children
    type(MeshMapToParent_3D), allocatable :: map_parent(:) !< map to parents

  contains

    procedure :: Init_Mesh_3D
    procedure :: Delete_Mesh_3D

    procedure :: GetCuboids
    procedure :: GetPoints
    procedure :: GetPointValency

    ! automatic identification and generation of components
    procedure :: BuildCommunicator
    procedure :: BuildFaces
    procedure :: BuildBoundaryFaces
    procedure :: BuildGhosts
    procedure :: BuildLinks
    procedure :: BuildCuboids
    procedure :: BuildMapToChild
    procedure :: BuildMapToParent
    procedure :: IdentifyClusters
    procedure :: IdentifyVertices
    procedure :: IdentifyEdges
    procedure :: IdentifyRanks

    ! import/export
    procedure :: ImportGenericMesh
    generic   :: ReadHDF5  => ReadHDF5_F, ReadHDF5_G
    generic   :: WriteHDF5 => WriteHDF5_F, WriteHDF5_G

    ! specific implementations are hided
    procedure, private :: ReadHDF5_F, ReadHDF5_G
    procedure, private :: WriteHDF5_F, WriteHDF5_G

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
      class(Mesh_3D),          intent(in)  :: mesh  !< mesh parition
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
      class(Mesh_3D), intent(inout) :: mesh !< mesh partition
    end subroutine BuildCuboids

    !---------------------------------------------------------------------------
    !> Generation of ghost elements

    module subroutine BuildGhosts(mesh)
      class(Mesh_3D), intent(inout) :: mesh !< mesh partition
    end subroutine BuildGhosts

    !---------------------------------------------------------------------------
    !> Generation of mesh faces

    module subroutine BuildFaces(mesh)
      class(Mesh_3D), intent(inout) :: mesh !< mesh partition
    end subroutine BuildFaces

    !---------------------------------------------------------------------------
    !> Generation of boundary faces

    module subroutine BuildBoundaryFaces(mesh)
      class(Mesh_3D), intent(inout) :: mesh !< mesh partition
    end subroutine BuildBoundaryFaces

    !---------------------------------------------------------------------------
    !> Generation of mesh links from global element neighbor information

    module subroutine BuildLinks(mesh)
      class(Mesh_3D), intent(inout) :: mesh !< mesh partition
    end subroutine BuildLinks

    !---------------------------------------------------------------------------
    !> Generation of map to child elements

    module subroutine BuildMapToChild(mesh)
      class(Mesh_3D), intent(inout) :: mesh !< mesh partition
    end subroutine BuildMapToChild

    !---------------------------------------------------------------------------
    !> Generation of map to child elements

    module subroutine BuildMapToParent(mesh)
      class(Mesh_3D), intent(inout) :: mesh !< mesh partition
    end subroutine BuildMapToParent

    !---------------------------------------------------------------------------
    !> Identification of sibling element clusters

    module subroutine IdentifyClusters(mesh)
      class(Mesh_3D), intent(inout) :: mesh !< mesh partition
    end subroutine IdentifyClusters

    !---------------------------------------------------------------------------
    !> Identification of mesh vertices

    module subroutine IdentifyVertices(mesh)
      class(Mesh_3D), intent(inout) :: mesh !< mesh partition
    end subroutine IdentifyVertices

    !---------------------------------------------------------------------------
    !> Generation of mesh edges using vertex IDs

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

    !---------------------------------------------------------------------------
    !> Read mesh partition from given HDF5 file

    module subroutine ReadHDF5_F(mesh, file, comm)
      class(Mesh_3D),   intent(inout) :: mesh !< mesh partition
      character(len=*), intent(in)    :: file !< name of HDF5 file
      type(MPI_Comm),   intent(in)    :: comm !< MPI "world" communicator
    end subroutine ReadHDF5_F

    !---------------------------------------------------------------------------
    !> Read mesh partition from given HDF5 group

    module subroutine ReadHDF5_G(mesh, group_id, comm)
      use HDF5_Binding
      class(Mesh_3D), target, intent(inout) :: mesh  !< mesh partition
      integer(HID_T), intent(in) :: group_id !< ID of related HDF5 group
      type(MPI_Comm), intent(in) :: comm     !< "world" communicator
    end subroutine ReadHDF5_G

    !---------------------------------------------------------------------------
    !> Write mesh partition into HDF5 file

    module subroutine WriteHDF5_F(mesh, file)
      class(Mesh_3D),   intent(in) :: mesh !< mesh partition
      character(len=*), intent(in) :: file !< name of HDF5 file
    end subroutine WriteHDF5_F

    !---------------------------------------------------------------------------
    !> Write mesh partition into given HDF5 group

    module subroutine WriteHDF5_G(mesh, group_id)
      use HDF5_Binding
      class(Mesh_3D), intent(in) :: mesh     !< mesh partition
      integer(HID_T), intent(in) :: group_id !< ID of related HDF5 group
    end subroutine WriteHDF5_G

  end interface

  !-----------------------------------------------------------------------------
  !> Type for collecting and transmitting the mesh attributes

  type MeshAttributes_3D

    integer   :: n_bound    = 0       !< number of domain boundaries
    integer   :: n_parts    = 0       !< number of non-empty partitions
    integer   :: p_geom     = 0       !< max polynomial order of element geometry
    logical   :: structured = .false. !< T if mapping to structured mesh exists
    logical   :: regular    = .false. !< T if equidistant Cartesian
    real(RNP) :: dx(3)      =  0      !< regular mesh spacing in directions 1:3
    logical   :: is_root    = .true.  !< T if root (bottom) level mesh
    logical   :: is_top     = .true.  !< T if top level mesh
    character :: refinement = ''      !< 'c' clone or 's' subdivision

    type(MeshBoundaryAttributes_3D), allocatable :: boundary(:)

  contains
    procedure :: Bcast => Bcast_MeshAttributes_3D
  end type MeshAttributes_3D

  ! constructor
  interface MeshAttributes_3D
    module procedure ExtractMeshAttributes
  end interface

  interface

    !---------------------------------------------------------------------------
    !> Get HDF5 datatype for essential static components of MeshAttributes_3D

    module subroutine Get_H5T_MeshAttributes_3D(H5T_MeshAttributes_3D)
      use HDF5_Binding
      integer(HID_T), intent(out) :: H5T_MeshAttributes_3D
    end subroutine Get_H5T_MeshAttributes_3D

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
    this % p_geom     = attrib % p_geom
    this % structured = attrib % structured
    this % regular    = attrib % regular
    this % dx         = attrib % dx
    this % is_root    = attrib % is_root
    this % is_top     = attrib % is_top
    this % refinement = attrib % refinement

    allocate(this % boundary( this%n_bound ))
    do b = 1, this % n_bound
      this % boundary(b) = MeshBoundary_3D(attrib % boundary(b))
    end do

    this % comm_world = comm

    call MPI_Comm_rank(comm, this % proc)

  end subroutine Init_Mesh_3D

  !-----------------------------------------------------------------------------
  !> Delete 3D mesh partition

  subroutine Delete_Mesh_3D(this)
    class(Mesh_3D), intent(inout) :: this

    ! reset static components ..................................................

    this % n_bound       =  0
    this % n_parts       =  0

    this % structured    = .false.
    this % regular       = .false.
    this % is_root       = .true.
    this % is_top        = .true.

    this % dx            =  0

    this % comm_world    = MPI_COMM_NULL
    this % comm_parts    = MPI_COMM_NULL

    this % proc          = -1
    this % part          = -1

    this % n_vert        =  0
    this % n_edge        =  0
    this % n_face        =  0
    this % n_elem        =  0
    this % n_elem_active =  0
    this % n_elem_frozen =  0
    this % n_cluster     =  0
    this % n_ghost       =  0
    this % n_link        =  0
    this % n_child       =  0
    this % n_parent      =  0

    this % p_geom        =  0

    this % max_vert_val  =  0
    this % max_edge_val  =  0

    this % n_elem_1      =  0
    this % n_elem_2      =  0
    this % n_elem_3      =  0

    ! release dynamic components ...............................................

    if (allocated( this % proc_part  )) deallocate( this % proc_part  )
    if (allocated( this % face       )) deallocate( this % face       )
    if (allocated( this % element    )) deallocate( this % element    )
    if (allocated( this % ghost      )) deallocate( this % ghost      )
    if (allocated( this % boundary   )) deallocate( this % boundary   )
    if (allocated( this % link       )) deallocate( this % link       )
    if (allocated( this % map_child  )) deallocate( this % map_child  )
    if (allocated( this % map_parent )) deallocate( this % map_parent )

  end subroutine Delete_Mesh_3D

  !=============================================================================
  ! MeshAttributes_3D constructor and type-bound procedures

  !-----------------------------------------------------------------------------
  !> Creates attributes by extraction from existing mesh partition

  function ExtractMeshAttributes(mesh) result(attrib)
    class(Mesh_3D), intent(in) :: mesh
    type(MeshAttributes_3D) :: attrib

    attrib % n_bound    = mesh % n_bound
    attrib % n_parts    = mesh % n_parts
    attrib % p_geom     = mesh % p_geom

    attrib % structured = mesh % structured
    attrib % regular    = mesh % regular
    attrib % dx         = mesh % dx

    attrib % is_root    = mesh % is_root
    attrib % is_top     = mesh % is_top
    attrib % refinement = mesh % refinement

    attrib % boundary = MeshBoundaryAttributes_3D(mesh % boundary)

  end function ExtractMeshAttributes

  !-----------------------------------------------------------------------------
  !> Broadcast mesh attributes

  subroutine Bcast_MeshAttributes_3D(this, root, comm)
    class(MeshAttributes_3D), intent(inout) :: this
    integer       , intent(in) :: root !< MPI root process
    type(MPI_Comm), intent(in) :: comm !< MPI communicator

    integer :: attrib_int(3), b
    logical :: attrib_log(4)

    ! mesh attributes ..........................................................

    attrib_int(1) = this % n_bound
    attrib_int(2) = this % n_parts
    attrib_int(3) = this % p_geom

    attrib_log(1) = this % structured
    attrib_log(2) = this % regular
    attrib_log(3) = this % is_root
    attrib_log(4) = this % is_top

    call XMPI_Bcast(attrib_int, root, comm)
    call XMPI_Bcast(attrib_log, root, comm)

    this % n_bound    = attrib_int(1)
    this % n_parts    = attrib_int(2)
    this % p_geom     = attrib_int(3)

    this % structured = attrib_log(1)
    this % regular    = attrib_log(2)
    this % is_root    = attrib_log(3)
    this % is_top     = attrib_log(4)

    call XMPI_Bcast(this % dx        , root, comm)
    call XMPI_Bcast(this % refinement, root, comm)

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
