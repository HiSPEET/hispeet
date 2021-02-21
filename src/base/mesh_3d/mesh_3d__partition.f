!> summary:  3D mesh partition type
!> author:   Joerg Stiller
!> date:     2020/11/13
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Mesh_3d__Partition
  use Kind_Parameters, only: IXL, IXS, RNP
  use XMPI
  use Mesh_3d__Face
  use Mesh_3d__Element
  use Mesh_3d__Boundary
  use Mesh_3d__Link
  implicit none
  private

  public :: Mesh3d_Partition

  !-----------------------------------------------------------------------------
  !> 3D mesh partition type
  !>
  !> ### Ghost elements
  !>
  !> The ghost elements stored in `ghost(1:n_ghost)` ...

  type Mesh3d_Partition

    ! dimensions
    integer :: n_vert  = 0           !< number of mesh vertices
    integer :: n_edge  = 0           !< number of mesh edges
    integer :: n_face  = 0           !< number of mesh faces
    integer :: n_elem  = 0           !< number of mesh elements
    integer :: n_ghost = 0           !< number of ghost elements
    integer :: p_geom  = 0           !< polynomial order of element geometry

    ! global attributes
    integer :: n_bound = 0           !< number of domain boundaries
    integer :: n_part  = 0           !< number of non-empty partitions

    ! local attributes
    integer :: part       = -1       !< partition ID
    logical :: structured = .false.  !< T if mapping to structured mesh exists
    logical :: regular    = .false.  !< T if equidistant Cartesian

    ! structured mesh properties
    integer :: n_elem_1 = 0          !< number of elements in direction 1
    integer :: n_elem_2 = 0          !< number of elements in direction 2
    integer :: n_elem_3 = 0          !< number of elements in direction 3

    ! regular mesh properties
    integer   :: n_face_1 =  0       !< number of faces normal to direction 1
    integer   :: n_face_2 =  0       !< number of faces normal to direction 2
    integer   :: n_face_3 =  0       !< number of faces normal to direction 3
    real(RNP) :: dx(3)    = -1       !< mesh element spacing in directions 1:3

    ! mesh components and links
    type(Mesh3d_Face)    , allocatable :: face(:)     !< mesh faces
    type(Mesh3d_Element) , allocatable :: element(:)  !< mesh elements
    type(Mesh3d_Element) , allocatable :: ghost(:)    !< ghost elements
    type(Mesh3d_Boundary), allocatable :: boundary(:) !< mesh boundaries
    type(Mesh3d_Link)    , allocatable :: link(:)     !< mesh links

    ! mesh element geometry
    real(RNP), allocatable :: x_elem(:,:,:,:,:) !< element GLL points

    ! MPI
    type(MPI_Comm) :: comm  !< communicator

  contains

    procedure :: GetPoints
    procedure :: GetPointValency

    ! automatic identification and generation of components
    procedure :: IdentifyEdges
    procedure :: IdentifyPrimaries
    procedure :: BuildFaces
    procedure :: BuildLinks
    procedure :: BuildGhosts

    ! import/export
    procedure :: ImportGenericMesh

  end type Mesh3d_Partition

  ! constructor
  interface Mesh3d_Partition
    module procedure EmptyMeshPartition
  end interface

  interface

    !---------------------------------------------------------------------------
    !> Generates Gauss-Lobatto or Gauss points to all elements of a partition

    module subroutine GetPoints(mesh, po, basis, x)
      class(Mesh3d_Partition), intent(in)  :: mesh         !< mesh parition
      integer,                 intent(in)  :: po           !< polynomial order
      character(len=*),        intent(in)  :: basis        !< 'GL' or 'GLL'
      real(RNP), allocatable,  intent(out) :: x(:,:,:,:,:) !< mesh points
    end subroutine GetPoints

    !---------------------------------------------------------------------------
    !> Computes the valency of mesh points

    module subroutine GetPointValency(mesh, v)
      class(Mesh3d_Partition), intent(in)  :: mesh          !< mesh parition
      integer,                 intent(out) :: v(0:,0:,0:,:) !< point valency
    end subroutine GetPointValency

    !---------------------------------------------------------------------------
    !> Identification of mesh edges
    !>
    !> Requires
    !>   - mesh % element % vertex % id
    !>
    !> Generates
    !>   - mesh % n_edge
    !>   - mesh % element % edge % id
    !>   - mesh % element % edge % orientation

    module subroutine IdentifyEdges(mesh)
      class(Mesh3d_Partition), intent(inout) :: mesh !< mesh partition
    end subroutine IdentifyEdges

    !---------------------------------------------------------------------------
    !> Generation of mesh faces
    !>
    !> Requires
    !>   - mesh % element % vertex % id
    !>   - mesh % element % edge % id
    !>
    !> Generates
    !>   - mesh % n_face
    !>   - mesh % face
    !>   - mesh % face % element % id
    !>   - mesh % face % element % face
    !>   - mesh % element % face % id
    !>   - mesh % element % face % orientation

    module subroutine BuildFaces(mesh)
      class(Mesh3d_Partition), intent(inout) :: mesh !< mesh partition
    end subroutine BuildFaces

    !---------------------------------------------------------------------------
    !> Identification of primary element components
    !>
    !> Requires
    !>   - mesh % element % face   % {n_neighbor, i_neighbor}
    !>   - mesh % element % edge   % {n_neighbor, i_neighbor}
    !>   - mesh % element % vertex % {n_neighbor, i_neighbor}
    !>
    !> Generates
    !>   - mesh % element % face   % primary
    !>   - mesh % element % edge   % primary
    !>   - mesh % element % vertex % primary

    module subroutine IdentifyPrimaries(mesh)
      class(Mesh3d_Partition), intent(inout) :: mesh !< mesh partition
    end subroutine IdentifyPrimaries

    !---------------------------------------------------------------------------
    !> Generation of mesh links from global element neighbor information
    !>
    !> On entry, mesh elements must be complete and `mesh%element%neighbor%id`
    !> set to the home-partition ID of the neighbors.
    !>
    !> Using this information
    !>
    !>   - `mesh % link` :
    !>      is built
    !>
    !>   – `mesh % element % neighbor % id` :
    !>      entries referring to remote neighbors are translated into ghost IDs
    !>
    !>   - `mesh % face % element` :
    !>     is completed by inserting the ghost ID of adjacent remote elements

    module subroutine BuildLinks(mesh)
      class(Mesh3d_Partition), intent(inout) :: mesh !< local partition
    end subroutine BuildLinks

    !---------------------------------------------------------------------------
    !> Generation of ghost elements
    !>
    !> On entry, the mesh elements and mesh links must be complete. Using this
    !> information, the ghosts are created in `mesh % ghost(1:n_ghost)` and
    !> initialized as follows:
    !>
    !>   - `ghost % global_id` :
    !>      is the global ID of the corresponding mesh element
    !>
    !>   - `ghost % local_id` :
    !>      is the virtual element ID in the local mesh partition. It holds
    !>      `ghost(i) % local_id = mesh % n_elem + i`
    !>

    module subroutine BuildGhosts(mesh)
      class(Mesh3d_Partition), intent(inout) :: mesh !< local partition
    end subroutine BuildGhosts

    !---------------------------------------------------------------------------
    !> Import a generic 3d mesh

    module subroutine ImportGenericMesh(mesh, generic_mesh, comm)
      use Generic_Mesh_3d
      class(Mesh3d_Partition), intent(out) :: mesh         !< mesh partition
      class(GenericMesh3d),    intent(in)  :: generic_mesh !< generic mesh
      type(MPI_Comm),          intent(in)  :: comm         !< MPI communicator
    end subroutine ImportGenericMesh

  end interface

contains

  !-----------------------------------------------------------------------------
  !> Creates an empty 3d mesh partition

  function EmptyMeshPartition( p_geom, n_bound, n_part, comm &
                             , structured, regular, dx       ) result(mesh)

    integer,             intent(in)  :: p_geom
    integer,             intent(in)  :: n_bound
    integer,             intent(in)  :: n_part
    type(MPI_Comm),      intent(in)  :: comm
    logical  , optional, intent(in)  :: structured
    logical  , optional, intent(in)  :: regular
    real(RNP), optional, intent(in)  :: dx(3)

    type(Mesh3d_Partition) :: mesh

    mesh % p_geom  = p_geom
    mesh % n_bound = n_bound
    mesh % n_part  = n_part
    mesh % comm    = comm

    if (present(structured)) mesh % structured = structured
    if (present(regular))    mesh % regular    = regular
    if (present(dx))         mesh % dx         = dx

    allocate( mesh % face     (0) )
    allocate( mesh % element  (0) )
    allocate( mesh % boundary (0) )
    allocate( mesh % link     (0) )

  end function EmptyMeshPartition

  !=============================================================================

end module Mesh_3d__Partition
