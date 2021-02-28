!> summary:  3D mesh element type
!> author:   Joerg Stiller
!> date:     2020/11/12
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Mesh_3d__Element
  use Kind_Parameters, only: IXL, IXS, RNP
  use Mesh_3d__Element_Indexing
  implicit none
  private

  public :: Mesh3d_Element
  public :: Mesh3d_ElementNeighbor

  !-----------------------------------------------------------------------------
  !> Element vertex data

  type Mesh3d_ElementVertex
    integer      :: id         = -1 !< local mesh vertex ID
    integer(IXS) :: n_neighbor =  0 !< number of neighbor elements
    integer(IXS) :: i_neighbor =  0 !< first entry in `neighbor` list
    integer(IXS) :: rank       =  0 !< rank among local EV ref to same mesh vert
    integer(IXS) :: val        =  0 !< vertex valency
  end type Mesh3d_ElementVertex

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

  type Mesh3d_ElementEdge
    integer      :: id          = -1 !< local mesh edge ID
    integer(IXS) :: orientation =  1 !< orientation against mesh edge
    integer(IXS) :: n_neighbor  =  0 !< number of neighbor elements
    integer(IXS) :: i_neighbor  =  0 !< first entry in `neighbor` list
    integer(IXS) :: rank        =  0 !< rank among local EE ref to same mesh edge
    integer(IXS) :: val         =  0 !< edge valency
  contains
    generic :: AlignWithMeshEdge    => AlignEdgeData_IDK, AlignEdgeData_RNP
    generic :: AlignWithElementEdge => AlignEdgeData_IDK, AlignEdgeData_RNP
    procedure, private :: AlignEdgeData_IDK, AlignEdgeData_RNP
  end type Mesh3d_ElementEdge

  !-----------------------------------------------------------------------------
  !> Element face data
  !>
  !> With `rotation` specifies the number of quarter rotations required to align
  !> the element face with the mesh face. In case of negative normal orientation,
  !> the face must be flipped before rotating.

  type Mesh3d_ElementFace
    integer      :: id         = -1 !< local mesh face ID
    integer      :: boundary   = -1 !< mesh boundary ID, if > 0
    integer(IXS) :: normal     =  1 !< normal orientation against mesh face (±1)
    integer(IXS) :: rotation   =  0 !< 1/4-rotation for aligning with mesh face
    integer(IXS) :: n_neighbor =  0 !< number of neighbor elements
    integer(IXS) :: i_neighbor =  0 !< first entry in `neighbor` list
    integer(IXS) :: rank       =  0 !< rank among local EF ref to same mesh face
    integer(IXS) :: val        =  0 !< face valency
  contains
    generic :: AlignWithMeshFace    => AlignWithMeshFace_IDK, &
                                       AlignWithMeshFace_IXS, &
                                       AlignWithMeshFace_RNP
    generic :: AlignWithElementFace => AlignWithElementFace_IDK, &
                                       AlignWithElementFace_IXS, &
                                       AlignWithElementFace_RNP
    procedure, private :: AlignWithMeshFace_IDK, AlignWithElementFace_IDK
    procedure, private :: AlignWithMeshFace_IXS, AlignWithElementFace_IXS
    procedure, private :: AlignWithMeshFace_RNP, AlignWithElementFace_RNP
  end type Mesh3d_ElementFace

  !-----------------------------------------------------------------------------
  !> Neighbor element properties
  !>
  !> The neighbor element properties are usually accessed through the components
  !> of a given element. Every component is coupled to a neighbor component of
  !> identical type. Thus, if the component is a face, then `cc` refers to the
  !> coupled element face of the neighbor, etc.
  !> If `part` coincides with the present partition, then `id` refers to a local
  !> element. Otherwise it corresponds to the ghost of a remote neighbor and can
  !> be used to access corresponding data.
  !>
  !> @note
  !> Before creating mesh links, `id` refers to the neighbor's home partition.

  type Mesh3d_ElementNeighbor
    integer :: id   = -1 !< local ID of neighbor element, including ghosts
    integer :: part = -1 !< partition owning the neighbor
    integer :: cc   = -1 !< coupled neighbor component
  end type Mesh3d_ElementNeighbor

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
  !>           1   |  1, 2     |  xi
  !>           2   |  3, 4     |  xi
  !>           3   |  5, 6     |  xi
  !>           4   |  7, 8     |  xi
  !>           5   |  1, 3     |  eta
  !>           6   |  2, 4     |  eta
  !>           7   |  5, 7     |  eta
  !>           8   |  6, 8     |  eta
  !>           9   |  1, 5     |  zeta
  !>          10   |  2, 6     |  zeta
  !>          11   |  3, 7     |  zeta
  !>          12   |  4, 8     |  zeta
  !>
  !>   - faces:
  !>
  !>         index |  vertices    |  edges          |  normal direction
  !>         ------|--------------|-----------------|------------------
  !>           1   |  1, 3, 5, 7  |  5,  7,  9, 11  |  xi
  !>           2   |  2, 4, 6, 8  |  6,  8, 10, 12  |  xi
  !>           3   |  1, 2, 5, 6  |  1,  3,  9, 10  |  eta
  !>           4   |  3, 4, 7, 8  |  2,  4, 11, 12  |  eta
  !>           5   |  1, 2, 3, 4  |  1,  2,  5,  6  |  zeta
  !>           6   |  5, 6, 7, 8  |  3,  4,  7,  8  |  zeta
  !>
  !> ### Ghost element mode
  !>
  !>

  type Mesh3d_Element

    integer(IXL) :: global_id = -1  !< global element ID
    integer      :: local_id  = -1  !< local  element ID

    type(Mesh3d_ElementVertex)  :: vertex(8)  !< vertex data
    type(Mesh3d_ElementEdge)    :: edge(12)   !< edge data
    type(Mesh3d_ElementFace)    :: face(6)    !< face data

    type(Mesh3d_ElementNeighbor), allocatable :: neighbor(:) !< neighbor data

  contains

    procedure, public :: DetermineFaceRank
    procedure, public :: DetermineEdgeRank
    procedure, public :: DetermineVertexRank

  end type Mesh3d_Element

contains

  !=============================================================================
  ! Edge alignment procedures
  !
  ! For convenience only. For better perfomance inline into application, which
  ! saves the procedure overhead.

  !-----------------------------------------------------------------------------
  !> Switch edge data between element and mesh orientations -- integer scalar

  pure subroutine AlignEdgeData_IDK(edge, v, va)
    class(Mesh3d_ElementEdge), intent(in) :: edge  !< mesh element edge
    integer, intent(in)  :: v(:)  !< given edge data
    integer, intent(out) :: va(:) !< aligned edge data

    integer :: i, np

    if (edge%orientation == 1_IXS) then ! edges aligned
      va = v
    else ! element and mesh edges have reverse orientation
      np = size(v,1)
      forall(i=1:np) va(np+1-i) = v(i)
    end if

  end subroutine AlignEdgeData_IDK

  !-----------------------------------------------------------------------------
  !> Switch edge data between element and mesh orientations -- real(RNP) scalar

  pure subroutine AlignEdgeData_RNP(edge, v, va)
    class(Mesh3d_ElementEdge), intent(in) :: edge  !< mesh element edge
    real(RNP), intent(in)  :: v(:)  !< given edge data
    real(RNP), intent(out) :: va(:) !< aligned edge data

    integer :: i, np

    if (edge%orientation == 1_IXS) then ! edges aligned
      va = v
    else ! element and mesh edges have reverse orientation
      np = size(v,1)
      do i = 1, np
        va(np+1-i) = v(i)
      end do
    end if

  end subroutine AlignEdgeData_RNP

  !=============================================================================
  ! Face alignment procedures

  !-----------------------------------------------------------------------------
  !> Transforms face data from element to mesh orientation -- integer scalar

  pure subroutine AlignWithMeshFace_IDK(face, ve, vm)
    class(Mesh3d_ElementFace), intent(in) :: face  !< mesh element face
    integer, intent(in)  :: ve(:,:)  !< element face data
    integer, intent(out) :: vm(:,:)  !< mesh face data

    integer :: i, j, l, np

    np = size(ve,1)
    l  = np + 1

    if (face % normal == 1_IXS) then
      select case(face % rotation)
      case(0_IXS)
        vm = ve
      case(1_IXS)
        forall (i=1:np, j=1:np)  vm(i,j) = ve(l-j,   i)
      case(2_IXS)
        forall (i=1:np, j=1:np)  vm(i,j) = ve(l-i, l-j)
      case default
        forall (i=1:np, j=1:np)  vm(i,j) = ve(  j, l-i)
      end select
    else
      select case(face % rotation)
      case(0_IXS)
        forall (i=1:np, j=1:np)  vm(i,j) = ve(  i, l-j)
      case(1_IXS)
        forall (i=1:np, j=1:np)  vm(i,j) = ve(l-j, l-i)
      case(2_IXS)
        forall (i=1:np, j=1:np)  vm(i,j) = ve(l-i,   j)
      case default
        forall (i=1:np, j=1:np)  vm(i,j) = ve(  j,   i)
      end select
    end if

  end subroutine AlignWithMeshFace_IDK

  !-----------------------------------------------------------------------------
  !> Transforms face data from element to mesh orientation -- integer scalar

  pure subroutine AlignWithMeshFace_IXS(face, ve, vm)
    class(Mesh3d_ElementFace), intent(in) :: face  !< mesh element face
    integer(IXS), intent(in)  :: ve(:,:)  !< element face data
    integer(IXS), intent(out) :: vm(:,:)  !< mesh face data

    integer :: i, j, l, np

    np = size(ve,1)
    l  = np + 1

    if (face % normal == 1_IXS) then
      select case(face % rotation)
      case(0_IXS)
        vm = ve
      case(1_IXS)
        forall (i=1:np, j=1:np)  vm(i,j) = ve(l-j,   i)
      case(2_IXS)
        forall (i=1:np, j=1:np)  vm(i,j) = ve(l-i, l-j)
      case default
        forall (i=1:np, j=1:np)  vm(i,j) = ve(  j, l-i)
      end select
    else
      select case(face % rotation)
      case(0_IXS)
        forall (i=1:np, j=1:np)  vm(i,j) = ve(  i, l-j)
      case(1_IXS)
        forall (i=1:np, j=1:np)  vm(i,j) = ve(l-j, l-i)
      case(2_IXS)
        forall (i=1:np, j=1:np)  vm(i,j) = ve(l-i,   j)
      case default
        forall (i=1:np, j=1:np)  vm(i,j) = ve(  j,   i)
      end select
    end if

  end subroutine AlignWithMeshFace_IXS

  !-----------------------------------------------------------------------------
  !> Transforms face data from element to mesh orientation -- real(RNP) scalar

  pure subroutine AlignWithMeshFace_RNP(face, ve, vm)
    class(Mesh3d_ElementFace), intent(in) :: face  !< mesh element face
    real(RNP), intent(in)  :: ve(:,:)  !< element face data
    real(RNP), intent(out) :: vm(:,:)  !< mesh face data

    integer :: i, j, l, np

    np = size(ve,1)
    l  = np + 1

    if (face % normal == 1_IXS) then
      select case(face % rotation)
      case(0_IXS)
        vm = ve
      case(1_IXS)
        forall (i=1:np, j=1:np)  vm(i,j) = ve(l-j,   i)
      case(2_IXS)
        forall (i=1:np, j=1:np)  vm(i,j) = ve(l-i, l-j)
      case default
        forall (i=1:np, j=1:np)  vm(i,j) = ve(  j, l-i)
      end select
    else
      select case(face % rotation)
      case(0_IXS)
        forall (i=1:np, j=1:np)  vm(i,j) = ve(  i, l-j)
      case(1_IXS)
        forall (i=1:np, j=1:np)  vm(i,j) = ve(l-j, l-i)
      case(2_IXS)
        forall (i=1:np, j=1:np)  vm(i,j) = ve(l-i,   j)
      case default
        forall (i=1:np, j=1:np)  vm(i,j) = ve(  j,   i)
      end select
    end if

  end subroutine AlignWithMeshFace_RNP

  !-----------------------------------------------------------------------------
  !> Transforms face data from mesh to element orientation -- integer scalar

  pure subroutine AlignWithElementFace_IDK(face, vm, ve)
    class(Mesh3d_ElementFace), intent(in) :: face  !< mesh element face
    integer, intent(in)  :: vm(:,:)  !< mesh face data
    integer, intent(out) :: ve(:,:)  !< element face data

    integer :: i, j, l, np

    np = size(vm,1)
    l  = np + 1

    if (face % normal == 1_IXS) then
      select case(face % rotation)
      case(0_IXS)
        ve = vm
      case(1_IXS)
        forall (i=1:np, j=1:np)  ve(l-j,   i) = vm(i,j)
      case(2_IXS)
        forall (i=1:np, j=1:np)  ve(l-i, l-j) = vm(i,j)
      case default
        forall (i=1:np, j=1:np)  ve(  j, l-i) = vm(i,j)
      end select
    else
      select case(face % rotation)
      case(0_IXS)
        forall (i=1:np, j=1:np)  ve(  i, l-j) = vm(i,j)
      case(1_IXS)
        forall (i=1:np, j=1:np)  ve(l-j, l-i) = vm(i,j)
      case(2_IXS)
        forall (i=1:np, j=1:np)  ve(l-i,   j) = vm(i,j)
      case default
        forall (i=1:np, j=1:np)  ve(  j,   i) = vm(i,j)
      end select
    end if

  end subroutine AlignWithElementFace_IDK

  !-----------------------------------------------------------------------------
  !> Transforms face data from mesh to element orientation -- integer scalar

  pure subroutine AlignWithElementFace_IXS(face, vm, ve)
    class(Mesh3d_ElementFace), intent(in) :: face  !< mesh element face
    integer(IXS), intent(in)  :: vm(:,:)  !< mesh face data
    integer(IXS), intent(out) :: ve(:,:)  !< element face data

    integer :: i, j, l, np

    np = size(vm,1)
    l  = np + 1

    if (face % normal == 1_IXS) then
      select case(face % rotation)
      case(0_IXS)
        ve = vm
      case(1_IXS)
        forall (i=1:np, j=1:np)  ve(l-j,   i) = vm(i,j)
      case(2_IXS)
        forall (i=1:np, j=1:np)  ve(l-i, l-j) = vm(i,j)
      case default
        forall (i=1:np, j=1:np)  ve(  j, l-i) = vm(i,j)
      end select
    else
      select case(face % rotation)
      case(0_IXS)
        forall (i=1:np, j=1:np)  ve(  i, l-j) = vm(i,j)
      case(1_IXS)
        forall (i=1:np, j=1:np)  ve(l-j, l-i) = vm(i,j)
      case(2_IXS)
        forall (i=1:np, j=1:np)  ve(l-i,   j) = vm(i,j)
      case default
        forall (i=1:np, j=1:np)  ve(  j,   i) = vm(i,j)
      end select
    end if

  end subroutine AlignWithElementFace_IXS

  !-----------------------------------------------------------------------------
  !> Transforms face data from mesh to element orientation -- real(RNP) scalar

  pure subroutine AlignWithElementFace_RNP(face, vm, ve)
    class(Mesh3d_ElementFace), intent(in) :: face  !< mesh element face
    real(RNP), intent(in)  :: vm(:,:)  !< mesh face data
    real(RNP), intent(out) :: ve(:,:)  !< element face data

    integer :: i, j, l, np

    np = size(vm,1)
    l  = np + 1

    if (face % normal == 1_IXS) then
      select case(face % rotation)
      case(0_IXS)
        ve = vm
      case(1_IXS)
        forall (i=1:np, j=1:np)  ve(l-j,   i) = vm(i,j)
      case(2_IXS)
        forall (i=1:np, j=1:np)  ve(l-i, l-j) = vm(i,j)
      case default
        forall (i=1:np, j=1:np)  ve(  j, l-i) = vm(i,j)
      end select
    else
      select case(face % rotation)
      case(0_IXS)
        forall (i=1:np, j=1:np)  ve(  i, l-j) = vm(i,j)
      case(1_IXS)
        forall (i=1:np, j=1:np)  ve(l-j, l-i) = vm(i,j)
      case(2_IXS)
        forall (i=1:np, j=1:np)  ve(l-i,   j) = vm(i,j)
      case default
        forall (i=1:np, j=1:np)  ve(  j,   i) = vm(i,j)
      end select
    end if

  end subroutine AlignWithElementFace_RNP

  !=============================================================================
  ! Procedures to determine the element component ranks
  !
  ! Note:
  ! These routines are required for mesh generation. Once completed the ranks
  ! are available via the corresponding components of the Mesh3d_ElementFace,
  ! Mesh3d_ElementEdge and Mesh3d_ElementVertex data structures.

  !-----------------------------------------------------------------------------
  !> Identifies the rank of an element face
  !>
  !> This routine determines the rank of face `f` among the all local element
  !> faces referring to the same mesh face. The `rank - 1` equals the number
  !> of neighbor elements with their local ID lower or equal than `l`, which
  !> defaults to `element % local_id`.
  !> Passing for `l` the local ID of an adjoining ghost element yields the rank
  !> of the corresponding face of the latter.
  !> If requested, the valency `val` is determined as well.

  module subroutine DetermineFaceRank(element, f, l, rank, val)
    class(Mesh3D_Element), intent(in)  :: element !< mesh element
    integer,               intent(in)  :: f       !< element face ID
    integer,     optional, intent(in)  :: l       !< reference element ID
    integer(IXS)         , intent(out) :: rank    !< element face rank
    integer(IXS),optional, intent(out) :: val     !< element face rank

    integer :: i, n
    integer :: l_

    if (present(l)) then
      l_ = l
    else
      l_ = element % local_id
    end if

    n = element % face(f) % n_neighbor
    i = element % face(f) % i_neighbor
    rank = 1_IXS + count(l_ >= element % neighbor(i:i+n-1) % id, kind=IXS)

print '(99(G0,1X))','f =',f,', l_ =',l_,', neighbor%id =',element%neighbor(i:i+n-1)%id

    if (present(val)) then
      val = int(1 + n, kind=IXS)
    end if

  end subroutine DetermineFaceRank

  !-----------------------------------------------------------------------------
  !> Identifies the rank of an element edge
  !>
  !> This routine determines the rank of edge `f` among the all local element
  !> edges referring to the same mesh edge. The `rank - 1` equals the number
  !> of neighbor elements with their local ID lower or equal than `l`, which
  !> defaults to `element % local_id`.
  !> Passing for `l` the local ID of an adjoining ghost element yields the rank
  !> of the corresponding edge of the latter.
  !> If requested, the valency `val` is determined as well.

  module subroutine DetermineEdgeRank(element, e, l, rank, val)
    class(Mesh3D_Element), intent(in)  :: element !< mesh element
    integer,               intent(in)  :: e       !< element edge ID
    integer,     optional, intent(in)  :: l       !< reference element ID
    integer(IXS)         , intent(out) :: rank    !< element edge rank
    integer(IXS),optional, intent(out) :: val     !< element edge valency

    integer :: i, k, f, n
    integer :: l_, val_

    if (present(l)) then
      l_ = l
    else
      l_ = element % local_id
    end if

    ! probe edge neighbors
    n = element % edge(e) % n_neighbor
    i = element % edge(e) % i_neighbor
    rank = 1_IXS + count(l_ >= element % neighbor(i:i+n-1) % id, kind=IXS)
    val_ = 1_IXS + n

    ! probe neighbors via adjoining faces
    do k = 1, 2
      f = F_EDGE(k,e)
      n = element % face(f) % n_neighbor
      i = element % face(f) % i_neighbor
      rank = rank + count(l_ >= element % neighbor(i:i+n-1) % id, kind=IXS)
      val_ = val_ + n
    end do

    if (present(val)) val = val_

  end subroutine DetermineEdgeRank

  !-----------------------------------------------------------------------------
  !> Identifies the rank of an element vertex
  !>
  !> This routine determines the rank of vertex `f` among the all local element
  !> vertexs referring to the same mesh vertex. The rank `r - 1` equals the
  !> number of neighbor elements with their local ID lower or equal than `l`,
  !> which defaults to `element % local_id`.
  !> Passing for `l` the local ID of an adjoining ghost element yields the rank
  !> of the corresponding vertex of the latter.
  !> If requested, the valency `val` is determined as well.

  module subroutine DetermineVertexRank(element, v, l, rank, val)
    class(Mesh3D_Element), intent(in)  :: element !< mesh element
    integer,               intent(in)  :: v       !< element vertex ID
    integer,     optional, intent(in)  :: l       !< reference element ID
    integer(IXS)         , intent(out) :: rank    !< element edge rank
    integer(IXS),optional, intent(out) :: val     !< element edge valency

    integer :: e, i, k, f, n
    integer :: l_, val_

    if (present(l)) then
      l_ = l
    else
      l_ = element % local_id
    end if

    ! probe vertex neighbors
    n = element % vertex(v) % n_neighbor
    i = element % vertex(v) % i_neighbor
    rank = 1_IXS + count(l_ >= element % neighbor(i:i+n-1) % id, kind=IXS)
    val_ = 1_IXS + n

    ! probe neighbors via adjoining edges
    do k = 1, 3
      e = E_VERT(k,v)
      n = element % edge(e) % n_neighbor
      i = element % edge(e) % i_neighbor
      rank = rank + count(l_ >= element % neighbor(i:i+n-1) % id, kind=IXS)
      val_ = val_ + n
    end do

    ! probe neighbors via adjoining faces
    do k = 1, 3
      f = F_VERT(k,v)
      n = element % face(f) % n_neighbor
      i = element % face(f) % i_neighbor
      rank = rank + count(l_ >= element % neighbor(i:i+n-1) % id, kind=IXS)
      val_ = val_ + n
    end do

    if (present(val)) val = val_

  end subroutine DetermineVertexRank

  !=============================================================================

end module Mesh_3d__Element
