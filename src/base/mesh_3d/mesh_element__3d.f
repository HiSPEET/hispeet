!> summary:  3D mesh element type
!> author:   Joerg Stiller
!> date:     2020/11/12
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Mesh_Element__3D
  use Kind_Parameters, only: IXL, IXS, RNP
  use Mesh_Element_Indexing__3D
  implicit none
  private

  public :: MeshElement_3D
  public :: MeshElementVertex_3D
  public :: MeshElementEdge_3D
  public :: MeshElementFace_3D
  public :: MeshElementNeighbor_3D

  !-----------------------------------------------------------------------------
  !> Element vertex data

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

  type MeshElementEdge_3D
    integer      :: id          = -1 !< local mesh edge ID
    integer(IXS) :: orientation =  1 !< orientation against mesh edge
    integer(IXS) :: n_neighbor  =  0 !< number of neighbor elements
    integer(IXS) :: i_neighbor  =  0 !< first entry in `neighbor` list
    integer(IXS) :: rank        =  0 !< rank among local EE ref to same mesh edge
    integer(IXS) :: val         =  0 !< edge valency
  contains
    generic :: AlignWithMesh    => AlignEdgeData_IDK, AlignEdgeData_RNP
    generic :: AlignWithElement => AlignEdgeData_IDK, AlignEdgeData_RNP
    procedure, private :: AlignEdgeData_IDK, AlignEdgeData_RNP
  end type MeshElementEdge_3D

  !-----------------------------------------------------------------------------
  !> Element face data
  !>
  !> With `rotation` specifies the number of quarter rotations required to align
  !> the element face with the mesh face. In case of negative normal orientation,
  !> the face must be flipped before rotating.

  type MeshElementFace_3D
    integer      :: id         = -1 !< local mesh face ID
    integer      :: boundary   = -1 !< mesh boundary ID, if > 0
    integer(IXS) :: normal     =  1 !< normal orientation against mesh face (±1)
    integer(IXS) :: rotation   =  0 !< 1/4-rotation for aligning with mesh face
    integer(IXS) :: n_neighbor =  0 !< number of neighbor elements
    integer(IXS) :: i_neighbor =  0 !< first entry in `neighbor` list
    integer(IXS) :: rank       =  0 !< rank among local EF ref to same mesh face
    integer(IXS) :: val        =  0 !< face valency
  contains
    procedure :: Side             =>  MeshFaceSide
    generic   :: AlignWithMesh    =>  AlignWithMeshFace_IDK, &
                                      AlignWithMeshFace_IXS, &
                                      AlignWithMeshFace_RNP
    generic   :: AlignWithElement =>  AlignWithElementFace_IDK, &
                                      AlignWithElementFace_IXS, &
                                      AlignWithElementFace_RNP
    procedure, private :: AlignWithMeshFace_IDK, AlignWithElementFace_IDK
    procedure, private :: AlignWithMeshFace_IXS, AlignWithElementFace_IXS
    procedure, private :: AlignWithMeshFace_RNP, AlignWithElementFace_RNP
  end type MeshElementFace_3D

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

  type MeshElementNeighbor_3D
    integer :: id   = -1 !< local ID of neighbor element, including ghosts
    integer :: part = -1 !< partition owning the neighbor
    integer :: cc   = -1 !< coupled neighbor component
  end type MeshElementNeighbor_3D

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

  type MeshElement_3D

    integer(IXL) :: global_id = -1  !< global element ID
    integer      :: local_id  = -1  !< local  element ID

    type(MeshElementVertex_3D)  :: vertex(8)  !< vertex data
    type(MeshElementEdge_3D)    :: edge(12)   !< edge data
    type(MeshElementFace_3D)    :: face(6)    !< face data

    type(MeshElementNeighbor_3D), allocatable :: neighbor(:) !< neighbor data

  end type MeshElement_3D

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
    class(MeshElementEdge_3D), intent(in) :: edge  !< mesh element edge
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
  ! MeshElementFace_3D procedures

  !-----------------------------------------------------------------------------
  !> Touched side of the adjacent mesh face

  pure integer function MeshFaceSide(face)
    class(MeshElementFace_3D), intent(in) :: face  !< mesh element face
    integer, parameter :: side(-1:1) = [2,0,1]
    MeshFaceSide = side(face % normal)
  end function MeshFaceSide

  !-----------------------------------------------------------------------------
  !> Transforms face data from element to mesh orientation -- integer scalar

  pure subroutine AlignWithMeshFace_IDK(face, ve, vm)
    class(MeshElementFace_3D), intent(in) :: face  !< mesh element face
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
    class(MeshElementFace_3D), intent(in) :: face  !< mesh element face
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
    class(MeshElementFace_3D), intent(in) :: face  !< mesh element face
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
    class(MeshElementFace_3D), intent(in) :: face  !< mesh element face
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
    class(MeshElementFace_3D), intent(in) :: face  !< mesh element face
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
    class(MeshElementFace_3D), intent(in) :: face  !< mesh element face
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

end module Mesh_Element__3D
