!> summary:  3D mesh element type
!> author:   Joerg Stiller
!> date:     2020/11/12
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Mesh_3d__Element
  use Kind_Parameters, only: IXL, IXS, RNP
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
    integer(IXS) :: primary    =  0 !< 0/1 if not/ first reference to mesh vertex
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
    integer(IXS) :: primary     =  0 !< 0/1 if not/ first reference to mesh edge
  contains
    generic :: AlignWithMeshEdge    => AlignEdgeData_IS, AlignEdgeData_RS
    generic :: AlignWithElementEdge => AlignEdgeData_IS, AlignEdgeData_RS
    procedure, private :: AlignEdgeData_IS, AlignEdgeData_RS
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
    integer(IXS) :: primary    =  0 !< 0/1 if not/ first reference to mesh face
  contains
    generic :: AlignWithMeshFace    => AlignWithMeshFace_IS, &
                                       AlignWithMeshFace_RS
    generic :: AlignWithElementFace => AlignWithElementFace_IS, &
                                       AlignWithElementFace_RS
    procedure, private :: AlignWithMeshFace_IS, AlignWithElementFace_IS
    procedure, private :: AlignWithMeshFace_RS, AlignWithElementFace_RS
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

  end type Mesh3d_Element

contains

  !=============================================================================
  ! Edge alignment procedures
  !
  ! For convenience only. For better perfomance inline into application, which
  ! saves the procedure overhead.

  !-----------------------------------------------------------------------------
  !> Switch edge data between element and mesh orientations -- integer scalar

  pure subroutine AlignEdgeData_IS(edge, v, va)
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

  end subroutine AlignEdgeData_IS

  !-----------------------------------------------------------------------------
  !> Switch edge data between element and mesh orientations -- real(RNP) scalar

  pure subroutine AlignEdgeData_RS(edge, v, va)
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

  end subroutine AlignEdgeData_RS

  !=============================================================================
  ! Face alignment procedures

  !-----------------------------------------------------------------------------
  !> Transforms face data from element to mesh orientation -- integer scalar

  pure subroutine AlignWithMeshFace_IS(face, ve, vm)
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

  end subroutine AlignWithMeshFace_IS

  !-----------------------------------------------------------------------------
  !> Transforms face data from element to mesh orientation -- real(RNP) scalar

  pure subroutine AlignWithMeshFace_RS(face, ve, vm)
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

  end subroutine AlignWithMeshFace_RS

  !-----------------------------------------------------------------------------
  !> Transforms face data from mesh to element orientation -- integer scalar

  pure subroutine AlignWithElementFace_IS(face, vm, ve)
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

  end subroutine AlignWithElementFace_IS

  !-----------------------------------------------------------------------------
  !> Transforms face data from mesh to element orientation -- real(RNP) scalar

  pure subroutine AlignWithElementFace_RS(face, vm, ve)
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

  end subroutine AlignWithElementFace_RS

  !=============================================================================

end module Mesh_3d__Element
