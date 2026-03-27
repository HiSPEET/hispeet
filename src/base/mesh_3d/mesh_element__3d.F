!> summary:  3D mesh element type
!> author:   Joerg Stiller
!> date:     2020/11/12
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> @remark
!>    Always update `Get_MPI_MeshElement_3D` when changing  `MeshElement_3D` or
!>    any of the subordinate types!
!===============================================================================

module Mesh_Element__3D
  use Kind_Parameters, only: IXS, RNP
  use Constants,       only: THIRD
  use XMPI
  use Mesh_Element_Indexing__3D
  implicit none
  private

  public :: MeshElement_3D
  public :: MeshElementVertex_3D
  public :: MeshElementEdge_3D
  public :: MeshElementFace_3D
  public :: MeshElementNeighbor_3D
  public :: MeshElementAdaptation_3D
  public :: MeshElementGeometry_3D

  public :: TransformIndex
  public :: TransformedEdgeAlignment
  public :: Get_MPI_MeshElement_3D
  public :: Get_H5T_MeshElement_3D
  public :: Get_H5T_MeshElementNeighbor_3D

  !-----------------------------------------------------------------------------
  !> Element vertex data
  !>
  !> The `rank` component provides an ordering of all element vertices referring
  !> to the same mesh vertex. Typically the former belong to differen elements.
  !> However, due to periodicity, several vertices of one element may coincide
  !> with the same mesh vertex.
  !>
  !> The valency `val` equals the total number of element vertices sharing the
  !> corresponding mesh vertex.
  !>
  !> Ghost element vertices coinciding with local vertices are provided with the
  !> following components
  !>   - `id`
  !>   - `rank`
  !>   - `val`

  type MeshElementVertex_3D
    integer      :: id         = -1 !< local mesh vertex ID
    integer(IXS) :: n_neighbor =  0 !< number of neighbor elements
    integer(IXS) :: i_neighbor =  0 !< first entry in `neighbor` list
    integer(IXS) :: rank       =  0 !< rank among local EV ref to same mesh vert
    integer(IXS) :: val        =  0 !< vertex valency
  end type MeshElementVertex_3D

  !-----------------------------------------------------------------------------
  !> Element edge data
  !>
  !> The `orientation` describes the transformation required to align the
  !> element edge with the mesh edge.
  !>
  !>     orientation | transformation
  !>     ------------|-----------------
  !>       1         | none
  !>      -1         | flip
  !>
  !> Ghost element edges coinciding with local edges are provided with the
  !> following components
  !>   - `id`
  !>   - `orientation`
  !>   - `rank`
  !>   - `val`

  type MeshElementEdge_3D
    integer      :: id          = -1 !< local mesh edge ID
    integer(IXS) :: orientation =  1 !< orientation against mesh edge
    integer(IXS) :: n_neighbor  =  0 !< number of neighbor elements
    integer(IXS) :: i_neighbor  =  0 !< first entry in `neighbor` list
    integer(IXS) :: rank        =  0 !< rank among local EE ref to same mesh edge
    integer(IXS) :: val         =  0 !< edge valency
  contains
    generic :: AlignWithMesh => AlignEdgeData_IDK, AlignEdgeData_RNP
    generic :: AlignFromMesh => AlignEdgeData_IDK, AlignEdgeData_RNP
    procedure, private :: AlignEdgeData_IDK, AlignEdgeData_RNP
  end type MeshElementEdge_3D

  !-----------------------------------------------------------------------------
  !> Element face data
  !>
  !> With `rotation` specifies the number of quarter rotations required to align
  !> the element face with the mesh face. In case of negative normal orientation,
  !> the face must be flipped before rotating.
  !>
  !> Ghost element faces coinciding with local faces are provided with the
  !> following components
  !>   - `id`
  !>   - `normal`
  !>   - `rotation`
  !>   - `rank`
  !>   - `val`

  type MeshElementFace_3D
    integer      :: id         = -1 !< local mesh face ID
    integer      :: boundary   = -1 !< boundary ID > 0, 0 for frozen neighbor
    integer(IXS) :: normal     =  1 !< normal orientation against mesh face (±1)
    integer(IXS) :: rotation   =  0 !< 1/4-rotation for aligning with mesh face
    integer(IXS) :: n_neighbor =  0 !< number of neighbor elements
    integer(IXS) :: i_neighbor =  0 !< first entry in `neighbor` list
    integer(IXS) :: rank       =  0 !< rank among local EF ref to same mesh face
    integer(IXS) :: val        =  0 !< face valency
  contains
    procedure :: SetOrientation =>  Set_FaceOrientation
    generic   :: AlignWithMesh  =>  AlignWithMeshFace_IDK, &
                                    AlignWithMeshFace_IXS, &
                                    AlignWithMeshFace_RNP
    generic   :: AlignFromMesh  =>  AlignFromMeshFace_IDK, &
                                    AlignFromMeshFace_IXS, &
                                    AlignFromMeshFace_RNP
    procedure, private :: AlignWithMeshFace_IDK, AlignFromMeshFace_IDK
    procedure, private :: AlignWithMeshFace_IXS, AlignFromMeshFace_IXS
    procedure, private :: AlignWithMeshFace_RNP, AlignFromMeshFace_RNP
  end type MeshElementFace_3D

  !-----------------------------------------------------------------------------
  !> Neighbor element properties
  !>
  !> The neighbor element properties are usually accessed through the components
  !> of a given element. Every component is coupled to a neighbor component of
  !> identical type.
  !> If `part` coincides with the present partition, then `id` refers to a local
  !> element. Otherwise it corresponds to the ghost of a remote neighbor and can
  !> be used to access corresponding data.
  !>
  !> The `orientation` component encodes the encodes the alignment of the
  !> coordinate system of the neighboring element as follows
  !>
  !>       orientation = 10 * d1 + d2
  !>
  !> where d1 and d2 indicate which neighbor directions are aligned with the
  !> ξ and η directions of the given element:
  !>
  !>     value | direction
  !>     ------|-----------
  !>       1   |   +ξ
  !>       2   |   +η
  !>       3   |   +ζ
  !>       4   |   -ξ
  !>       5   |   -η
  !>       6   |   -ζ
  !>
  !> @note
  !> Before creating mesh links, `id` refers to the neighbor's home partition.

  type MeshElementNeighbor_3D
    integer      :: id   = -1 !< local ID of neighbor element, including ghosts
    integer      :: part = -1 !< partition owning the neighbor
    integer(IXS) :: component   = -1 !< coupled component {1:26}
    integer(IXS) :: orientation = -1 !< neighbor orientation
  end type MeshElementNeighbor_3D
  !
  ! overview of possible orientations
  !
  ! neighbor                           edge alignment     reverse
  ! orientation     d1   d2   d3       ex    ey    ez     orientation
  !
  !     12           x    y    z        +     +     +        12
  !     13           x    z   -y        +     +     -        16
  !     15           x   -y   -z        +     -     -        15
  !     16           x   -z    y        +     -     +        13
  !
  !     21           y    x   -z        +     +     -        21
  !     23           y    z    x        +     +     +        31
  !     24           y   -x    z        +     -     +        51
  !     26           y   -z   -x        +     -     -        61
  !
  !     31           z    x    y        +     +     +        23
  !     32           z    y   -x        +     +     -        62
  !     34           z   -x   -y        +     -     -        56
  !     35           z   -y    x        +     -     +        35
  !
  !     42          -x    y   -z        -     +     -        42
  !     43          -x    z    y        -     +     +        43
  !     45          -x   -y    z        -     -     +        45
  !     46          -x   -z   -y        -     -     -        46
  !
  !     51          -y    x    z        -     +     +        24
  !     53          -y    z   -x        -     +     -        64
  !     54          -y   -x   -z        -     -     -        54
  !     56          -y   -z    x        -     -     +        34
  !
  !     61          -z    x   -y        -     +     -        26
  !     62          -z    y    x        -     +     +        32
  !     64          -z   -x    y        -     -     +        53
  !     65          -z   -y   -x        -     -     -        65

  !-----------------------------------------------------------------------------
  !> Element adaptation information
  !>
  !> Two types of adaptation are supported:
  !>   - 'c':  global or zonal cloning of elements
  !>   - 's':  global of local subdivision of elements
  !>
  !> The elements selected for cloning or subdividing generate 1 or 8 "active"
  !> children in the next higher mesh level. In addition, the surrounding
  !> elements create "frozen" children enclosing the adaptation zone.
  !>
  !> The `refinement` component is used to identify the the adaptation type,
  !> number and state of the children and the associated element component:
  !>
  !>   | value | type  | component | children |  state  |
  !>   | -----:|:-----:|:--------- |:--------:| ------- |
  !>   |     0 |       | none      |     0    |         |
  !>   |   100 |  'c'  | volume    |     1    | frozen  |
  !>   |   101 |  's'  | vertex 1  |     1    | frozen  |
  !>   |    .. |  ..'  | ..        |    ..    |   ..    |
  !>   |   108 |  's'  | vertex 8  |     1    | frozen  |
  !>   |   201 |  's'  | edge 1    |     2    | frozen  |
  !>   |    .. |  ..'  | ..        |    ..    |   ..    |
  !>   |   212 |  's'  | edge 12   |     2    | frozen  |
  !>   |   401 |  's'  | face 1    |     4    | frozen  |
  !>   |    .. |  ..'  | ..        |    ..    |   ..    |
  !>   |   406 |  's'  | face 6    |     4    | frozen  |
  !>   |   800 |  's'  | volume    |     8    | frozen  |
  !>   |  1000 |  'c'  | volume    |     1    | active  |
  !>   |  8000 |  's'  | volume    |     8    | active  |
  !>
  !> `mark` is used for different purposes. Before adaptation, it is set to
  !> indicate the intended refinement of the existing elements. In course of
  !> the adaptation process it may take different roles to assist the creation
  !> of new elements or to remap sustained data.

  type MeshElementAdaptation_3D
    integer :: parent_proc   = -1 !< parent partition ID in comm_world
    integer :: parent_id     = -1 !< parent element ID
    integer :: refinement    = -1 !< refinement type
    integer :: child_proc    = -1 !< child partition ID in comm_world
    integer :: sublevels     = -1 !< number of active sublevels
    integer :: mark          = -1 !< adaptation mark or auxiliary data
  end type MeshElementAdaptation_3D

  !-----------------------------------------------------------------------------
  !> Element geometry data
  !>
  !> The element shape is defined by the Lobatto points `x_e(0:po,0:po,0:po,3)`,
  !> where `po` is the order of the corresponding Lagrange polynomials.
  !> Transforming the Lagrange form into a Legendre series and restricting the
  !> latter to the trilinear part yields the approximate cuboid, which is given
  !> by
  !>
  !>       x(ξ,η,ζ,:) = x_c(0,:) + x_c(1,:) ξ + x_c(2,:) η + x_c(3,:) ζ
  !>
  !> with `-1 ≤ ξ,η,ζ ≤ 1`.
  !>
  !> The array `dx_m` provides the harmonic mean spacing across the element
  !> faces based on the corresponding cuboids. Given the normal spacings `dx₁`
  !> and `dx₂` of the cuboids adjacent to a given face, the mean spacing is
  !> defined as
  !>
  !>       dx_m = 2 / (1/dx₁ + 1/dx₂)
  !>
  !> Note that mean spacing is a unique property of the shared mesh face, which
  !> is stored element-wise for convenience.

  type MeshElementGeometry_3D
    integer                :: po           !< polynomial order
    real(RNP), allocatable :: x_e(:,:,:,:) !< element Lobatto points
    real(RNP)              :: x_c(0:3,1:3) !< corresponding cuboid
    real(RNP)              :: dx_m(6)      !< mean element spacing across faces
  end type MeshElementGeometry_3D

  !-----------------------------------------------------------------------------
  !> 3D mesh element
  !>
  !> ### Numbering
  !>
  !>   - vertices:
  !>
  !>         index |  standard vertex coordinates
  !>         ------|-----------------------------
  !>           1   |  -1,-1,-1
  !>           2   |   1,-1,-1
  !>           3   |  -1, 1,-1
  !>           4   |   1, 1,-1
  !>           5   |  -1,-1, 1
  !>           6   |   1,-1, 1
  !>           7   |  -1, 1, 1
  !>           8   |   1, 1, 1
  !>
  !>   - edges:
  !>
  !>         index |  vertices |  direction
  !>         ------|-----------|-----------
  !>           1   |  1, 2     |  ξ
  !>           2   |  3, 4     |  ξ
  !>           3   |  5, 6     |  ξ
  !>           4   |  7, 8     |  ξ
  !>           5   |  1, 3     |  η
  !>           6   |  2, 4     |  η
  !>           7   |  5, 7     |  η
  !>           8   |  6, 8     |  η
  !>           9   |  1, 5     |  ζ
  !>          10   |  2, 6     |  ζ
  !>          11   |  3, 7     |  ζ
  !>          12   |  4, 8     |  ζ
  !>
  !>   - faces:
  !>
  !>         index |  vertices    |  edges          |  normal direction
  !>         ------|--------------|-----------------|------------------
  !>           1   |  1, 3, 5, 7  |  5,  7,  9, 11  |  ξ
  !>           2   |  2, 4, 6, 8  |  6,  8, 10, 12  |  ξ
  !>           3   |  1, 2, 5, 6  |  1,  3,  9, 10  |  η
  !>           4   |  3, 4, 7, 8  |  2,  4, 11, 12  |  η
  !>           5   |  1, 2, 3, 4  |  1,  2,  5,  6  |  ζ
  !>           6   |  5, 6, 7, 8  |  3,  4,  7,  8  |  ζ
  !>
  !> ### Ghost element mode
  !>
  !> Ghost elements are equipped with a reduced set of information for all
  !> faces, edges and vertices that coincide with local entities.
  !> A detailed description of the provided data is given in the documentation
  !> of the corresponding data types.
  !>
  !> No neighbor and geometry data is available for ghosts.
  !>
  !> ### SFC rank and path
  !>
  !> If a space filling curve (SFC) is available, the elements are ordered
  !> according to their position on the curve. This yields a global ranking
  !> of the active elements, which is stored in `sfc_rank`.
  !> Inside the element, the SCF follows a path through all children, from
  !> `i` to `j`. This path is encoded in `sfc_path = 10*i + j`.

  type MeshElement_3D

    integer :: id          = -1 !< local element ID
    integer :: cluster_id  = -1 !< sibling cluster ID
    integer :: cluster_oct = -1 !< sibling cluster octant (1:8)

    integer :: sfc_rank    = -1 !< rank according to position on SFC
    integer :: sfc_path    = -1 !< SFC path within element

    logical :: frozen = .false. !< switch indicating inactive elements

    type(MeshElementVertex_3D)  :: vertex(8) !< vertex data
    type(MeshElementEdge_3D)    :: edge(12)  !< edge data
    type(MeshElementFace_3D)    :: face(6)   !< face data

    type(MeshElementNeighbor_3D), allocatable :: neighbor(:) !< neighbor data

    type(MeshElementAdaptation_3D) :: adaptation !< adaptation info
    type(MeshElementGeometry_3D)   :: geometry   !< geometry data

  contains

    procedure :: IsActive
    procedure :: IsLeaf

    procedure :: GetCuboidDimensions

    procedure :: MeshFaceSide
    procedure :: NeigborEdgeIsAligned

    generic   :: AlignFromNeighborEdge => AlignFromNeighborEdge_IDK, &
                                          AlignFromNeighborEdge_RNP

    generic   :: AlignFromNeighborFace => AlignFromNeighborFace_IDK, &
                                          AlignFromNeighborFace_RNP

    generic   :: AlignWithNeighborEdge => AlignWithNeighborEdge_IDK, &
                                          AlignWithNeighborEdge_RNP

    generic   :: AlignWithNeighborFace => AlignWithNeighborFace_IDK, &
                                          AlignFromNeighborFace_RNP
    procedure :: MarkForRefinement
    procedure :: MarkForRemoval
    procedure :: Unmark

    procedure, private :: AlignFromNeighborEdge_IDK
    procedure, private :: AlignFromNeighborEdge_RNP
    procedure, private :: AlignFromNeighborFace_IDK
    procedure, private :: AlignFromNeighborFace_RNP

    procedure, private :: AlignWithNeighborEdge_IDK
    procedure, private :: AlignWithNeighborEdge_RNP
    procedure, private :: AlignWithNeighborFace_IDK
    procedure, private :: AlignWithNeighborFace_RNP

  end type MeshElement_3D

  !=============================================================================
  ! Interfaces to submodule procedures

  interface

    !---------------------------------------------------------------------------
    !> Get MPI datatype for mesh elements
    !>
    !> The datatype allows to transfer the static components of elements, e.g.
    !>
    !>     type(MeshElement_3D) :: element(n)
    !>     type(MPI_Datatype)   :: MPI_MeshElement_3D
    !>     ...
    !>     call Get_MPI_MeshElement_3D(MPI_MeshElement_3D)
    !>     call MPI_Send(element(1)%id, n, MPI_MeshElement_3D, dest, tag, comm)

    module subroutine Get_MPI_MeshElement_3D(MPI_MeshElement_3D)
      type(MPI_Datatype), intent(out) :: MPI_MeshElement_3D
    end subroutine Get_MPI_MeshElement_3D

    !---------------------------------------------------------------------------
    !> Get HDF5 datatype for essential static components of MeshElement_3D

    module subroutine Get_H5T_MeshElement_3D(H5T_MeshElement_3D)
      use HDF5_Binding
      integer(hid_t), intent(out) :: H5T_MeshElement_3D
    end subroutine Get_H5T_MeshElement_3D

    !---------------------------------------------------------------------------
    !> Get HDF5 datatype of MeshElementNeighbor_3D

    module subroutine Get_H5T_MeshElementNeighbor_3D(H5T_MeshElementNeighbor_3D)
      use HDF5_Binding
      integer(hid_t), intent(out) :: H5T_MeshElementNeighbor_3D
    end subroutine Get_H5T_MeshElementNeighbor_3D

  end interface

contains

  !=============================================================================
  ! Edge alignment procedures
  !
  ! For convenience only. For better perfomance inline into application, which
  ! saves the procedure overhead.

  !-----------------------------------------------------------------------------
  !> Switch edge data between element and mesh orientations -- integer scalar

  pure subroutine AlignEdgeData_IDK(edge, v, va)
    class(MeshElementEdge_3D), intent(in) :: edge  !< mesh element edge
    integer, intent(in)  :: v(:)  !< given edge data
    integer, intent(out) :: va(:) !< aligned edge data

    integer :: i, m

    if (edge%orientation == 1_IXS) then ! edges aligned
      va = v
    else ! element and mesh edges have reverse orientation
      m = size(v,1)
      forall(i=1:m) va(m+1-i) = v(i)
    end if

  end subroutine AlignEdgeData_IDK

  !-----------------------------------------------------------------------------
  !> Switch edge data between element and mesh orientations -- real(RNP) scalar

  pure subroutine AlignEdgeData_RNP(edge, v, va)
    class(MeshElementEdge_3D), intent(in) :: edge  !< mesh element edge
    real(RNP), intent(in)  :: v(:)  !< given edge data
    real(RNP), intent(out) :: va(:) !< aligned edge data

    integer :: i, m

    if (edge%orientation == 1_IXS) then ! edges aligned
      va = v
    else ! element and mesh edges have reverse orientation
      m = size(v,1)
      do i = 1, m
        va(m+1-i) = v(i)
      end do
    end if

  end subroutine AlignEdgeData_RNP

  !=============================================================================
  ! MeshElementFace_3D procedures

  !-----------------------------------------------------------------------------
  !> Set element face normal and rotation with respect to mesh face

  pure subroutine Set_FaceOrientation(this, efv, mfv)
    class(MeshElementFace_3D), intent(inout) :: this
    integer, intent(in) :: efv(4) !< element-face vertex IDs
    integer, intent(in) :: mfv(4) !< mesh-face vertex IDs

    associate(normal => this % normal, rotation => this % rotation)

      if (mfv(1) == efv(1)) then
        if (mfv(2) == efv(2)) then
          normal   =  1
          rotation =  0
        else
          normal   = -1
          rotation =  3
        end if

      else if (mfv(1) == efv(2)) then
        if (mfv(2) == efv(4)) then
          normal   =  1
          rotation =  1
        else
          normal   = -1
          rotation =  2
        end if

      else if (mfv(1) == efv(4)) then
        if (mfv(2) == efv(3)) then
          normal   =  1
          rotation =  2
        else
          normal   = -1
          rotation =  1
        end if

      else ! mfv(1) == efv(3)
        if (mfv(2) == efv(1)) then
          normal   =  1
          rotation =  3
        else
          normal   = -1
          rotation =  0
        end if

      end if

    end associate

  end subroutine Set_FaceOrientation

  !-----------------------------------------------------------------------------
  !> Transforms face data from element to mesh orientation -- integer scalar

  pure subroutine AlignWithMeshFace_IDK(face, ve, vm)
    class(MeshElementFace_3D), intent(in) :: face  !< mesh element face
    integer, intent(in)  :: ve(:,:)  !< element face data
    integer, intent(out) :: vm(:,:)  !< mesh face data

#include "mesh__3d__mp_alignwithmeshface-body.f"

  end subroutine AlignWithMeshFace_IDK

  !-----------------------------------------------------------------------------
  !> Transforms face data from element to mesh orientation -- integer(IXS) scalar

  pure subroutine AlignWithMeshFace_IXS(face, ve, vm)
    class(MeshElementFace_3D), intent(in) :: face  !< mesh element face
    integer(IXS), intent(in)  :: ve(:,:)  !< element face data
    integer(IXS), intent(out) :: vm(:,:)  !< mesh face data

#include "mesh__3d__mp_alignwithmeshface-body.f"

  end subroutine AlignWithMeshFace_IXS

  !-----------------------------------------------------------------------------
  !> Transforms face data from element to mesh orientation -- real(RNP) scalar

  pure subroutine AlignWithMeshFace_RNP(face, ve, vm)
    class(MeshElementFace_3D), intent(in) :: face  !< mesh element face
    real(RNP), intent(in)  :: ve(:,:)  !< element face data
    real(RNP), intent(out) :: vm(:,:)  !< mesh face data

#include "mesh__3d__mp_alignwithmeshface-body.f"

  end subroutine AlignWithMeshFace_RNP

  !-----------------------------------------------------------------------------
  !> Transforms face data from mesh to element orientation -- integer scalar

  pure subroutine AlignFromMeshFace_IDK(face, vm, ve)
    class(MeshElementFace_3D), intent(in) :: face  !< mesh element face
    integer, intent(in)  :: vm(:,:)  !< mesh face data
    integer, intent(out) :: ve(:,:)  !< element face data

#include "mesh__3d__mp_alignfrommeshface-body.f"

  end subroutine AlignFromMeshFace_IDK

  !-----------------------------------------------------------------------------
  !> Transforms face data from mesh to element orientation -- integer(IXS) scalar

  pure subroutine AlignFromMeshFace_IXS(face, vm, ve)
    class(MeshElementFace_3D), intent(in) :: face  !< mesh element face
    integer(IXS), intent(in)  :: vm(:,:)  !< mesh face data
    integer(IXS), intent(out) :: ve(:,:)  !< element face data

#include "mesh__3d__mp_alignfrommeshface-body.f"

  end subroutine AlignFromMeshFace_IXS

  !-----------------------------------------------------------------------------
  !> Transforms face data from mesh to element orientation -- real(RNP) scalar

  pure subroutine AlignFromMeshFace_RNP(face, vm, ve)
    class(MeshElementFace_3D), intent(in) :: face  !< mesh element face
    real(RNP), intent(in)  :: vm(:,:)  !< mesh face data
    real(RNP), intent(out) :: ve(:,:)  !< element face data

#include "mesh__3d__mp_alignfrommeshface-body.f"

  end subroutine AlignFromMeshFace_RNP

  !=============================================================================
  ! MeshElement_3D procedures

  !-----------------------------------------------------------------------------
  !> Checks whether the element is active

  elemental logical function IsActive(element)
    class(MeshElement_3D), intent(in) :: element !< mesh element

    IsActive = .not. element % frozen

  end function IsActive

  !-----------------------------------------------------------------------------
  !> Checks whether the element is a leaf

  elemental logical function IsLeaf(element)
    class(MeshElement_3D), intent(in) :: element !< mesh element

    IsLeaf = .not. element % frozen &
             .and. element % adaptation % refinement < 1000

  end function IsLeaf

  !-----------------------------------------------------------------------------
  !> Compute the dimensions of the approximate cuboid

  pure subroutine GetCuboidDimensions(element, dx)
    class(MeshElement_3D), intent(in) :: element !< mesh element
    real(RNP), intent(out) :: dx(3) !< cuboid length in directions ξ₁,ξ₂,ξ₃

    real(RNP) :: vol

    associate(x_c => element % geometry % x_c)

      ! side lengths of the parallelepiped .....................................

      dx(1) = 2 * sqrt( x_c(1,1) ** 2 + x_c(1,2) ** 2 + x_c(1,3) ** 2 )
      dx(2) = 2 * sqrt( x_c(2,1) ** 2 + x_c(2,2) ** 2 + x_c(2,3) ** 2 )
      dx(3) = 2 * sqrt( x_c(3,1) ** 2 + x_c(3,2) ** 2 + x_c(3,3) ** 2 )

      ! volume of the parallelepiped ............................................

      vol = 8 * ( x_c(1,1) * (x_c(2,2) * x_c(3,3) - x_c(2,3)*x_c(3,2)) &
                + x_c(1,2) * (x_c(2,3) * x_c(3,1) - x_c(2,1)*x_c(3,3)) &
                + x_c(1,3) * (x_c(2,1) * x_c(3,2) - x_c(2,2)*x_c(3,1)) )

      ! side lengths of approximate cuboid with identical volume ...............

      dx = dx * (vol / product(dx)) ** THIRD

    end associate

  end subroutine GetCuboidDimensions

  !-----------------------------------------------------------------------------
  !> Side of the adjacent mesh face on which the element is located

  pure integer function MeshFaceSide(element, f) result(side)
   class(MeshElement_3D), intent(in) :: element !< mesh element
   integer, intent(in) :: f  !< element face ID

   select case(f * element%face(f)%normal)
   case(1, 3, 5, -2, -4, -6) ! mesh face normal points into element
     side = 2
   case default
     side = 1
   end select

  end function MeshFaceSide

  !-----------------------------------------------------------------------------
  !> Checks if the neighbor is aligned to given element edge

  pure logical function NeigborEdgeIsAligned(element, e, n) result(aligned)
    class(MeshElement_3D), intent(in) :: element  !< mesh element
    integer, intent(in)  :: e !< element edge ID
    integer, intent(in)  :: n !< element neighbor ID

    integer :: d1, d2

    d1 = element % neighbor(n) % orientation / 10
    d2 = element % neighbor(n) % orientation - 10 * d1

    if (e <= 4) then ! ξ-edge

      aligned = d1 <= 3

    else if (e <= 8) then  ! η-edge

      aligned = d2 <= 3

    else  ! ζ-edge

      select case(d1)
      case(1)
        aligned = d2 == 2 .or. d2 == 6
      case(2)
        aligned = d2 == 3 .or. d2 == 4
      case(3)
        aligned = d2 == 1 .or. d2 == 5
      case(4)
        aligned = d2 == 5 .or. d2 == 3
      case(5)
        aligned = d2 == 6 .or. d2 == 1
      case(6)
        aligned = d2 == 4 .or. d2 == 2
      end select

    end if

  end function NeigborEdgeIsAligned

  !-----------------------------------------------------------------------------
  !> Transforms neighbor edge data to element orientation -- integer scalar

  pure subroutine AlignFromNeighborEdge_IDK(element, e, n, vn, ve)
    class(MeshElement_3D), intent(in) :: element  !< mesh element
    integer, intent(in)  :: e     !< element edge ID
    integer, intent(in)  :: n     !< element neighbor ID
    integer, intent(in)  :: vn(:) !< neighbor edge data
    integer, intent(out) :: ve(:) !< element edge data

    integer :: i, m

    if (element % NeigborEdgeIsAligned(e, n)) then
      ve = vn
    else
      m = size(vn)
      do i = 1, m
        ve(m+1-i) = vn(i)
      end do
    end if

  end subroutine AlignFromNeighborEdge_IDK

  !-----------------------------------------------------------------------------
  !> Transforms neighbor edge data to element orientation -- real(RNP) scalar

  pure subroutine AlignFromNeighborEdge_RNP(element, e, n, vn, ve)
    class(MeshElement_3D), intent(in) :: element  !< mesh element
    integer  , intent(in)  :: e     !< element edge ID
    integer  , intent(in)  :: n     !< element neighbor ID
    real(RNP), intent(in)  :: vn(:) !< neighbor edge data
    real(RNP), intent(out) :: ve(:) !< element edge data

    integer :: i, m

    if (element % NeigborEdgeIsAligned(e, n)) then
      ve = vn
    else
      m = size(vn)
      do i = 1, m
        ve(m+1-i) = vn(i)
      end do
    end if

  end subroutine AlignFromNeighborEdge_RNP

  !-----------------------------------------------------------------------------
  !> Transforms neighbor face data to element orientation -- integer scalar

  pure subroutine AlignFromNeighborFace_IDK(element, f, n, vn, ve)
    class(MeshElement_3D), intent(in) :: element  !< mesh element
    integer, intent(in)  :: f        !< element face ID
    integer, intent(in)  :: n        !< element neighbor ID
    integer, intent(in)  :: vn(:,:)  !< neighbor face data
    integer, intent(out) :: ve(:,:)  !< element face data

#include "mesh__3d__mp_alignfromneighborface-body.f"

  end subroutine AlignFromNeighborFace_IDK

  !-----------------------------------------------------------------------------
  !> Transforms neighbor face data to element orientation -- real(RNP) scalar

  pure subroutine AlignFromNeighborFace_RNP(element, f, n, vn, ve)
    class(MeshElement_3D), intent(in) :: element  !< mesh element
    integer,   intent(in)  :: f        !< element face ID
    integer,   intent(in)  :: n        !< element neighbor ID
    real(RNP), intent(in)  :: vn(:,:)  !< neighbor face data
    real(RNP), intent(out) :: ve(:,:)  !< element face data

#include "mesh__3d__mp_alignfromneighborface-body.f"

  end subroutine AlignFromNeighborFace_RNP

  !-----------------------------------------------------------------------------
  !> Transforms element edge data to neighbor orientation -- integer scalar

  pure subroutine AlignWithNeighborEdge_IDK(element, e, n, ve, vn)
    class(MeshElement_3D), intent(in) :: element  !< mesh element
    integer, intent(in)  :: e     !< element edge ID
    integer, intent(in)  :: n     !< element neighbor ID
    integer, intent(in)  :: ve(:) !< element edge data
    integer, intent(out) :: vn(:) !< neighbor edge data

    integer :: i, m

    if (element % NeigborEdgeIsAligned(e, n)) then
      vn = ve
    else
      m = size(vn)
      do i = 1, m
        vn(m+1-i) = ve(i)
      end do
    end if

  end subroutine AlignWithNeighborEdge_IDK

  !-----------------------------------------------------------------------------
  !> Transforms element edge data to neighbor orientation -- real(RNP) scalar

  pure subroutine AlignWithNeighborEdge_RNP(element, e, n, ve, vn)
    class(MeshElement_3D), intent(in) :: element  !< mesh element
    integer  , intent(in)  :: e     !< element edge ID
    integer  , intent(in)  :: n     !< element neighbor ID
    real(RNP), intent(in)  :: ve(:) !< element edge data
    real(RNP), intent(out) :: vn(:) !< neighbor edge data

    integer :: i, m

    if (element % NeigborEdgeIsAligned(e, n)) then
      vn = ve
    else
      m = size(vn)
      do i = 1, m
        vn(m+1-i) = ve(i)
      end do
    end if

  end subroutine AlignWithNeighborEdge_RNP

  !-----------------------------------------------------------------------------
  !> Transforms neighbor face data to element orientation -- integer scalar

  pure subroutine AlignWithNeighborFace_IDK(element, f, n, ve, vn)
    class(MeshElement_3D), intent(in) :: element  !< mesh element
    integer, intent(in)  :: f        !< element face ID
    integer, intent(in)  :: n        !< element neighbor ID
    integer, intent(in)  :: ve(:,:)  !< element face data
    integer, intent(out) :: vn(:,:)  !< neighbor face data

#include "mesh__3d__mp_alignwithneighborface-body.f"

  end subroutine AlignWithNeighborFace_IDK

  !-----------------------------------------------------------------------------
  !> Transforms neighbor face data to element orientation -- real(RNP) scalar

  pure subroutine AlignWithNeighborFace_RNP(element, f, n, ve, vn)
    class(MeshElement_3D), intent(in) :: element  !< mesh element
    integer,   intent(in)  :: f        !< element face ID
    integer,   intent(in)  :: n        !< element neighbor ID
    real(RNP), intent(in)  :: ve(:,:)  !< element face data
    real(RNP), intent(out) :: vn(:,:)  !< neighbor face data

#include "mesh__3d__mp_alignwithneighborface-body.f"

  end subroutine AlignWithNeighborFace_RNP

  !-----------------------------------------------------------------------------
  !> Mark element for refinement

  elemental subroutine MarkForRefinement(element)
    class(MeshElement_3D), intent(inout) :: element  !< mesh element

    element % adaptation % mark = 1000

  end subroutine MarkForRefinement

  !-----------------------------------------------------------------------------
  !> Mark element for removal

  elemental subroutine MarkForRemoval(element)
    class(MeshElement_3D), intent(inout) :: element  !< mesh element

    element % adaptation % mark = -1

  end subroutine MarkForRemoval

  !-----------------------------------------------------------------------------
  !> Unmark element

  elemental subroutine Unmark(element)
    class(MeshElement_3D), intent(inout) :: element  !< mesh element

    element % adaptation % mark = 0

  end subroutine Unmark

  !-----------------------------------------------------------------------------
  !> Transforms the given triple index range into a rotated coordinate system

  pure subroutine TransformIndex(orientation, n, o0, o, t0, t, s)

    integer(IXS), intent(in) :: orientation  !< orientation of transformed coords
    integer, intent(in)  :: n                !< number of items per direction
    integer, intent(in)  :: o0               !< start of original indices
    integer, intent(in)  :: o(3)             !< original indices
    integer, intent(in)  :: t0               !< start of transformed indices
    integer, intent(out) :: t(3)             !< transformed indices
    integer, intent(out), optional :: s(3,3) !< strides of transformed indices

    integer :: m, s_(3,3)

    m = n - 1

    ! transformed indices
    select case(orientation)

    case(12_IXS) ! 1 →  1, 2 →  2:  (  1,  2,  3) ← (  1,  2,  3)
      t(1) = t0 + (o(1) - o0)     ;  s_(:,1) = [  1,  0,  0 ]
      t(2) = t0 + (o(2) - o0)     ;  s_(:,2) = [  0,  1,  0]
      t(3) = t0 + (o(3) - o0)     ;  s_(:,3) = [  0,  0,  1]
    case(13_IXS) ! 1 →  1, 2 →  3:  (  1,  2,  3) ← (  1, -3,  2)
      t(1) = t0 + (o(1) - o0)     ;  s_(:,1) = [  1,  0,  0 ]
      t(2) = t0 - (o(3) - o0) + m ;  s_(:,2) = [  0,  0, -1 ]
      t(3) = t0 + (o(2) - o0)     ;  s_(:,3) = [  0,  1,  0]
    case(15_IXS) ! 1 →  1, 2 → -2:  (  1,  2,  3) ← (  1, -2, -3)
      t(1) = t0 + (o(1) - o0)     ;  s_(:,1) = [  1,  0,  0 ]
      t(2) = t0 - (o(2) - o0) + m ;  s_(:,2) = [  0, -1,  0 ]
      t(3) = t0 - (o(3) - o0) + m ;  s_(:,3) = [  0,  0, -1]
    case(16_IXS) ! 1 →  1, 2 → -3:  (  1,  2,  3) ← (  1,  3, -2)
      t(1) = t0 + (o(1) - o0)     ;  s_(:,1) = [  1,  0,  0 ]
      t(2) = t0 + (o(3) - o0)     ;  s_(:,2) = [  0,  0,  1 ]
      t(3) = t0 - (o(2) - o0) + m ;  s_(:,3) = [  0, -1,  0 ]

    case(21_IXS) ! 1 →  2, 2 →  1:  (  1,  2,  3) ← (  2,  1, -3)
      t(1) = t0 + (o(2) - o0)     ;  s_(:,1) = [  0,  1,  0 ]
      t(2) = t0 + (o(1) - o0)     ;  s_(:,2) = [  1,  0,  0 ]
      t(3) = t0 - (o(3) - o0) + m ;  s_(:,3) = [  0,  0, -1 ]
    case(23_IXS) ! 1 →  2, 2 →  3:  (  1,  2,  3) ← (  3,  1,  2)
      t(1) = t0 + (o(3) - o0)     ;  s_(:,1) = [  0,  0,  1 ]
      t(2) = t0 + (o(1) - o0)     ;  s_(:,2) = [  1,  0,  0 ]
      t(3) = t0 + (o(2) - o0)     ;  s_(:,3) = [  0,  1,  0 ]
    case(24_IXS) ! 1 →  2, 2 → -1:  (  1,  2,  3) ← ( -2,  1,  3)
      t(1) = t0 - (o(2) - o0) + m ;  s_(:,1) = [  0, -1,  0 ]
      t(2) = t0 + (o(1) - o0)     ;  s_(:,2) = [  1,  0,  0 ]
      t(3) = t0 + (o(3) - o0)     ;  s_(:,3) = [  0,  0,  1 ]
    case(26_IXS) ! 1 →  2, 2 → -3:  (  1,  2,  3) ← ( -3,  1, -2)
      t(1) = t0 - (o(3) - o0) + m ;  s_(:,1) = [  0,  0, -1 ]
      t(2) = t0 + (o(1) - o0)     ;  s_(:,2) = [  1,  0,  0 ]
      t(3) = t0 - (o(2) - o0) + m ;  s_(:,3) = [  0, -1,  0 ]

    case(31_IXS) ! 1 →  3, 2 →  1:  (  1,  2,  3) ← (  2,  3,  1)
      t(1) = t0 + (o(2) - o0)     ;  s_(:,1) = [  0,  1,  0 ]
      t(2) = t0 + (o(3) - o0)     ;  s_(:,2) = [  0,  0,  1 ]
      t(3) = t0 + (o(1) - o0)     ;  s_(:,3) = [  1,  0,  0 ]
    case(32_IXS) ! 1 →  3, 2 →  2:  (  1,  2,  3) ← ( -3,  2,  1)
      t(1) = t0 - (o(3) - o0) + m ;  s_(:,1) = [  0,  0, -1 ]
      t(2) = t0 + (o(2) - o0)     ;  s_(:,2) = [  0,  1,  0 ]
      t(3) = t0 + (o(1) - o0)     ;  s_(:,3) = [  1,  0,  0 ]
    case(34_IXS) ! 1 →  3, 2 → -1:  (  1,  2,  3) ← ( -2, -3,  1)
      t(1) = t0 - (o(2) - o0) + m ;  s_(:,1) = [  0, -1,  0 ]
      t(2) = t0 - (o(3) - o0) + m ;  s_(:,2) = [  0,  0, -1 ]
      t(3) = t0 + (o(1) - o0)     ;  s_(:,3) = [  1,  0,  0 ]
    case(35_IXS) ! 1 →  3, 2 → -2:  (  1,  2,  3) ← (  3, -2,  1)
      t(1) = t0 + (o(3) - o0)     ;  s_(:,1) = [  0,  0,  1 ]
      t(2) = t0 - (o(2) - o0) + m ;  s_(:,2) = [  0, -1,  0 ]
      t(3) = t0 + (o(1) - o0)     ;  s_(:,3) = [  1,  0,  0 ]

    case(42_IXS) ! 1 → -1, 2 →  2:  (  1,  2,  3) ← ( -1,  2, -3)
      t(1) = t0 - (o(1) - o0) + m ;  s_(:,1) = [ -1,  0,  0 ]
      t(2) = t0 + (o(2) - o0)     ;  s_(:,2) = [  0,  1,  0 ]
      t(3) = t0 - (o(3) - o0) + m ;  s_(:,3) = [  0,  0, -1 ]
    case(43_IXS) ! 1 → -1, 2 →  3:  (  1,  2,  3) ← ( -1,  3,  2)
      t(1) = t0 - (o(1) - o0) + m ;  s_(:,1) = [ -1,  0,  0 ]
      t(2) = t0 + (o(3) - o0)     ;  s_(:,2) = [  0,  0,  1 ]
      t(3) = t0 + (o(2) - o0)     ;  s_(:,3) = [  0,  1,  0 ]
    case(45_IXS) ! 1 → -1, 2 → -2:  (  1,  2,  3) ← ( -1, -2,  3)
      t(1) = t0 - (o(1) - o0) + m ;  s_(:,1) = [ -1,  0,  0 ]
      t(2) = t0 - (o(2) - o0) + m ;  s_(:,2) = [  0, -1,  0 ]
      t(3) = t0 + (o(3) - o0)     ;  s_(:,3) = [  0,  0,  1 ]
    case(46_IXS) ! 1 → -1, 2 → -3:  (  1,  2,  3) ← ( -1, -3, -2)
      t(1) = t0 - (o(1) - o0) + m ;  s_(:,1) = [ -1,  0,  0 ]
      t(2) = t0 - (o(3) - o0) + m ;  s_(:,2) = [  0,  0, -1 ]
      t(3) = t0 - (o(2) - o0) + m ;  s_(:,3) = [  0, -1,  0 ]

    case(51_IXS) ! 1 → -2, 2 →  1:  (  1,  2,  3) ← (  2, -1,  3)
      t(1) = t0 + (o(2) - o0)     ;  s_(:,1) = [  0,  1,  0 ]
      t(2) = t0 - (o(1) - o0) + m ;  s_(:,2) = [ -1,  0,  0 ]
      t(3) = t0 + (o(3) - o0)     ;  s_(:,3) = [  0,  0,  1 ]
    case(53_IXS) ! 1 → -2, 2 →  3:  (  1,  2,  3) ← ( -3, -1,  2)
      t(1) = t0 - (o(3) - o0) + m ;  s_(:,1) = [  0,  0, -1 ]
      t(2) = t0 - (o(1) - o0) + m ;  s_(:,2) = [ -1,  0,  0 ]
      t(3) = t0 + (o(2) - o0)     ;  s_(:,3) = [  0,  1,  0 ]
    case(54_IXS) ! 1 → -2, 2 → -1:  (  1,  2,  3) ← ( -2, -1, -3)
      t(1) = t0 - (o(2) - o0) + m ;  s_(:,1) = [  0, -1,  0 ]
      t(2) = t0 - (o(1) - o0) + m ;  s_(:,2) = [ -1,  0,  0 ]
      t(3) = t0 - (o(3) - o0) + m ;  s_(:,3) = [  0,  0, -1 ]
    case(56_IXS) ! 1 → -2, 2 → -3:  (  1,  2,  3) ← (  3, -1, -2)
      t(1) = t0 + (o(3) - o0)     ;  s_(:,1) = [  0,  0,  1 ]
      t(2) = t0 - (o(1) - o0) + m ;  s_(:,2) = [ -1,  0,  0 ]
      t(3) = t0 - (o(2) - o0) + m ;  s_(:,3) = [  0, -1,  0 ]

    case(61_IXS) ! 1 → -3, 2 →  1:  (  1,  2,  3) ← (  2, -3, -1)
      t(1) = t0 + (o(2) - o0)     ;  s_(:,1) = [  0,  1,  0 ]
      t(2) = t0 - (o(3) - o0) + m ;  s_(:,2) = [  0,  0, -1 ]
      t(3) = t0 - (o(1) - o0) + m ;  s_(:,3) = [ -1,  0,  0 ]
    case(62_IXS) ! 1 → -3, 2 →  2:  (  1,  2,  3) ← (  3,  2, -1)
      t(1) = t0 + (o(3) - o0)     ;  s_(:,1) = [  0,  0,  1 ]
      t(2) = t0 + (o(2) - o0)     ;  s_(:,2) = [  0,  1,  0 ]
      t(3) = t0 - (o(1) - o0) + m ;  s_(:,3) = [ -1,  0,  0 ]
    case(64_IXS) ! 1 → -3, 2 → -1:  (  1,  2,  3) ← ( -2,  3, -1)
      t(1) = t0 - (o(2) - o0) + m ;  s_(:,1) = [  0, -1,  0 ]
      t(2) = t0 + (o(3) - o0)     ;  s_(:,2) = [  0,  0,  1 ]
      t(3) = t0 - (o(1) - o0) + m ;  s_(:,3) = [ -1,  0,  0 ]
    case(65_IXS) ! 1 → -3, 2 → -2:  (  1,  2,  3) ← ( -3, -2, -1)
      t(1) = t0 - (o(3) - o0) + m ;  s_(:,1) = [  0,  0, -1 ]
      t(2) = t0 - (o(2) - o0) + m ;  s_(:,2) = [  0, -1,  0 ]
      t(3) = t0 - (o(1) - o0) + m ;  s_(:,3) = [ -1,  0,  0 ]

    end select

    if (present(s)) s = s_

  end subroutine TransformIndex

  !-----------------------------------------------------------------------------
  !> Alignment of an element edge with respect to rotated coordinate system
  !>
  !> Returns `+1` if the edge is aligned with one of the transformed directions
  !> and `-1` if the corresponding direction has the opposite sense.

  pure integer function TransformedEdgeAlignment(orientation, e) result(a)

    integer(IXS), intent(in) :: orientation  !< orientation of transformed coords
    integer,      intent(in) :: e            !< element edge ID

    select case(e)

    case(1:4) ! ξ edges
      select case(int(orientation))
      case(12,13,15,16,21,23,24,26,31,32,34,35)
        a = +1
      case default
        a = -1
      end select

    case(5:8) ! η edges
      select case(int(orientation))
      case(12,13,21,23,31,32,42,43,51,53,61,62)
        a = +1
      case default
        a = -1
      end select

    case default ! ζ edges
      select case(int(orientation))
      case(12,16,23,24,31,35,43,45,51,56,62,64)
        a = +1
      case default
        a = -1
      end select

    end select

  end function TransformedEdgeAlignment

  !=============================================================================

end module Mesh_Element__3D
