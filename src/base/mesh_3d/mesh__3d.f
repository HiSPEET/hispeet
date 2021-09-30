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

  !-----------------------------------------------------------------------------
  !> 3D mesh partition type
  !>
  !> ### Ghost elements
  !>
  !> The ghost elements stored in `ghost(1:n_ghost)` ...
  !>
  !> ### Element geometry
  !>
  !> All element domains are defined as Lagrange polynomials of degree `p_geom`.
  !> The array `x_elem(0:p_geom,0:p_geom,0:p_geom,l,1:3)` holds the Cartesian
  !> coordinates at the Lobatto points in element `l`. For convenience, `x_cube`
  !> provides a cuboidal approximation to that element, which is given by
  !>
  !>       x(ξ,η,ζ) = x_cube(0,l,1:3)
  !>                + x_cube(1,l,1:3) ξ
  !>                + x_cube(2,l,1:3) η
  !>                + x_cube(3,l,1:3) ζ
  !>
  !> with `-1 ≤ ξ,η,ζ ≤ 1`.

  type Mesh_3D

    ! global attributes ........................................................

    integer :: n_bound = 0           !< number of domain boundaries
    integer :: n_part  = 0           !< number of non-empty partitions

    ! local attributes .........................................................

    integer :: part       = -1       !< partition ID
    logical :: structured = .false.  !< T if mapping to structured mesh exists
    logical :: regular    = .false.  !< T if equidistant Cartesian

    ! dimensions
    integer :: n_vert  = 0           !< number of mesh vertices
    integer :: n_edge  = 0           !< number of mesh edges
    integer :: n_face  = 0           !< number of mesh faces
    integer :: n_elem  = 0           !< number of mesh elements
    integer :: n_ghost = 0           !< number of ghost elements
    integer :: p_geom  = 0           !< polynomial order of element geometry

    integer :: max_vert_val = 0      !< maximum vertex valency
    integer :: max_edge_val = 0      !< maximum edge valency

    ! structured mesh properties
    integer :: n_elem_1 = 0          !< number of elements in direction 1
    integer :: n_elem_2 = 0          !< number of elements in direction 2
    integer :: n_elem_3 = 0          !< number of elements in direction 3

    ! regular mesh properties
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

    ! mesh element geometry
    real(RNP), allocatable :: x_elem(:,:,:,:,:) !< element Lobatto points
    real(RNP), allocatable :: x_cube(:,:,:)     !< approximate cuboids

    ! MPI
    type(MPI_Comm) :: comm  !< communicator

  contains

    procedure :: GetCuboids
    procedure :: GetPoints
    procedure :: GetPointValency

    ! automatic identification and generation of components
    procedure :: IdentifyEdges
    procedure :: IdentifyRanks
    procedure :: BuildFaces
    procedure :: BuildLinks
    procedure :: BuildGhosts

    ! import/export
    procedure :: ImportGenericMesh

  end type Mesh_3D

  ! constructor
  interface Mesh_3D
    module procedure EmptyMeshPartition
  end interface

  interface

    !---------------------------------------------------------------------------
    !> Generates Cuboid points to all elements of a partition

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

    !---------------------------------------------------------------------------
    !> Identification of mesh edges

    module subroutine IdentifyEdges(mesh)
      class(Mesh_3D), intent(inout) :: mesh !< mesh partition
    end subroutine IdentifyEdges

    !---------------------------------------------------------------------------
    !> Generation of mesh faces

    module subroutine BuildFaces(mesh)
      class(Mesh_3D), intent(inout) :: mesh !< mesh partition
    end subroutine BuildFaces

    !---------------------------------------------------------------------------
    !> Identification of primary element components

    module subroutine IdentifyRanks(mesh)
      class(Mesh_3D), intent(inout) :: mesh !< mesh partition
    end subroutine IdentifyRanks

    !---------------------------------------------------------------------------
    !> Generation of mesh links from global element neighbor information

    module subroutine BuildLinks(mesh)
      class(Mesh_3D), intent(inout) :: mesh !< local partition
    end subroutine BuildLinks

    !---------------------------------------------------------------------------
    !> Generation of ghost elements

    module subroutine BuildGhosts(mesh)
      class(Mesh_3D), intent(inout) :: mesh !< local partition
    end subroutine BuildGhosts

    !---------------------------------------------------------------------------
    !> Import a generic 3d mesh

    module subroutine ImportGenericMesh(mesh, generic_mesh, comm)
      use Generic_Mesh__3D
      class(Mesh_3D),        intent(out) :: mesh         !< mesh partition
      class(GenericMesh_3D), intent(in)  :: generic_mesh !< generic mesh
      type(MPI_Comm),        intent(in)  :: comm         !< MPI communicator
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

    type(Mesh_3D) :: mesh

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

end module Mesh__3D
