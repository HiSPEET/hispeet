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

  public :: TransformIndex

  !-----------------------------------------------------------------------------
  !> Element vertex data
  !>
  !> The orientation flag `o_neighbor` indicates the direction of the neighbor
  !> parallel to the ξ direction of the current element. This feature is only
  !> available for regular vertices enclosed by 8 elements. In other cases
  !> `o_neighbor` is set to 0.
  !>
  !> The `rank` component provides an ordering of all element vertices referring
  !> to the same mesh vertex. Typically the former belong to differen elements.
  !> However, due to periodicity, several vertices of one element may coincide
  !> with the same mesh vertex.
  !>
  !> The `valency` equals the total number of element vertices sharing the
  !> corresponding mesh vertex.

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
    procedure :: Side          =>  MeshFaceSide
    generic   :: AlignWithMesh =>  AlignWithMeshFace_IDK, &
                                   AlignWithMeshFace_IXS, &
                                   AlignWithMeshFace_RNP
    generic   :: AlignFromMesh =>  AlignFromMeshFace_IDK, &
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

  !-----------------------------------------------------------------------------
  !> Element adaptation information

  type MeshElementAdaptation_3D
    integer      :: parent_proc = -1 !< parent partition ID in comm_world
    integer      :: parent_id   = -1 !< parent element ID
    integer      :: child_proc  = -1 !< child partition ID in comm_world
    integer      :: child_id    = -1 !< first child element ID
    integer(IXS) :: refinement  = -1 !< refinement type: -1,0,1 = none,copy,full
    integer(IXS) :: levels      = -1 !< number of subdivision levels
    integer(IXS) :: mark        = -1 !< adaptation mark (targeted refinement)
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
  !>       dx_me = 2 / (1/dx₁ + 1/dx₂)
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
  !>

  type MeshElement_3D

    integer(IXL) :: global_id = -1  !< global element ID
    integer      :: local_id  = -1  !< local  element ID

    type(MeshElementVertex_3D)  :: vertex(8)  !< vertex data
    type(MeshElementEdge_3D)    :: edge(12)   !< edge data
    type(MeshElementFace_3D)    :: face(6)    !< face data

    type(MeshElementNeighbor_3D), allocatable :: neighbor(:) !< neighbor data

    type(MeshElementAdaptation_3D) :: adaptation !< adaptation info
    type(MeshElementGeometry_3D)   :: geometry   !< geometry data

  contains

    procedure :: AlignFromNeighborFace

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

    integer :: i, j, l, m

    m = size(ve,1)
    l = m + 1

    if (face % normal == 1_IXS) then
      select case(face % rotation)
      case(0_IXS)
        vm = ve
      case(1_IXS)
        forall (i=1:m, j=1:m)  vm(i,j) = ve(l-j,   i)
      case(2_IXS)
        forall (i=1:m, j=1:m)  vm(i,j) = ve(l-i, l-j)
      case default
        forall (i=1:m, j=1:m)  vm(i,j) = ve(  j, l-i)
      end select
    else
      select case(face % rotation)
      case(0_IXS)
        forall (i=1:m, j=1:m)  vm(i,j) = ve(  i, l-j)
      case(1_IXS)
        forall (i=1:m, j=1:m)  vm(i,j) = ve(l-j, l-i)
      case(2_IXS)
        forall (i=1:m, j=1:m)  vm(i,j) = ve(l-i,   j)
      case default
        forall (i=1:m, j=1:m)  vm(i,j) = ve(  j,   i)
      end select
    end if

  end subroutine AlignWithMeshFace_IDK

  !-----------------------------------------------------------------------------
  !> Transforms face data from element to mesh orientation -- integer scalar

  pure subroutine AlignWithMeshFace_IXS(face, ve, vm)
    class(MeshElementFace_3D), intent(in) :: face  !< mesh element face
    integer(IXS), intent(in)  :: ve(:,:)  !< element face data
    integer(IXS), intent(out) :: vm(:,:)  !< mesh face data

    integer :: i, j, l, m

    m = size(ve,1)
    l = m + 1

    if (face % normal == 1_IXS) then
      select case(face % rotation)
      case(0_IXS)
        vm = ve
      case(1_IXS)
        forall (i=1:m, j=1:m)  vm(i,j) = ve(l-j,   i)
      case(2_IXS)
        forall (i=1:m, j=1:m)  vm(i,j) = ve(l-i, l-j)
      case default
        forall (i=1:m, j=1:m)  vm(i,j) = ve(  j, l-i)
      end select
    else
      select case(face % rotation)
      case(0_IXS)
        forall (i=1:m, j=1:m)  vm(i,j) = ve(  i, l-j)
      case(1_IXS)
        forall (i=1:m, j=1:m)  vm(i,j) = ve(l-j, l-i)
      case(2_IXS)
        forall (i=1:m, j=1:m)  vm(i,j) = ve(l-i,   j)
      case default
        forall (i=1:m, j=1:m)  vm(i,j) = ve(  j,   i)
      end select
    end if

  end subroutine AlignWithMeshFace_IXS

  !-----------------------------------------------------------------------------
  !> Transforms face data from element to mesh orientation -- real(RNP) scalar

  pure subroutine AlignWithMeshFace_RNP(face, ve, vm)
    class(MeshElementFace_3D), intent(in) :: face  !< mesh element face
    real(RNP), intent(in)  :: ve(:,:)  !< element face data
    real(RNP), intent(out) :: vm(:,:)  !< mesh face data

    integer :: i, j, l, m

    m = size(ve,1)
    l = m + 1

    if (face % normal == 1_IXS) then
      select case(face % rotation)
      case(0_IXS)
        vm = ve
      case(1_IXS)
        forall (i=1:m, j=1:m)  vm(i,j) = ve(l-j,   i)
      case(2_IXS)
        forall (i=1:m, j=1:m)  vm(i,j) = ve(l-i, l-j)
      case default
        forall (i=1:m, j=1:m)  vm(i,j) = ve(  j, l-i)
      end select
    else
      select case(face % rotation)
      case(0_IXS)
        forall (i=1:m, j=1:m)  vm(i,j) = ve(  i, l-j)
      case(1_IXS)
        forall (i=1:m, j=1:m)  vm(i,j) = ve(l-j, l-i)
      case(2_IXS)
        forall (i=1:m, j=1:m)  vm(i,j) = ve(l-i,   j)
      case default
        forall (i=1:m, j=1:m)  vm(i,j) = ve(  j,   i)
      end select
    end if

  end subroutine AlignWithMeshFace_RNP

  !-----------------------------------------------------------------------------
  !> Transforms face data from mesh to element orientation -- integer scalar

  pure subroutine AlignFromMeshFace_IDK(face, vm, ve)
    class(MeshElementFace_3D), intent(in) :: face  !< mesh element face
    integer, intent(in)  :: vm(:,:)  !< mesh face data
    integer, intent(out) :: ve(:,:)  !< element face data

    integer :: i, j, l, m

    m = size(vm,1)
    l = m + 1

    if (face % normal == 1_IXS) then
      select case(face % rotation)
      case(0_IXS)
        ve = vm
      case(1_IXS)
        forall (i=1:m, j=1:m)  ve(l-j,   i) = vm(i,j)
      case(2_IXS)
        forall (i=1:m, j=1:m)  ve(l-i, l-j) = vm(i,j)
      case default
        forall (i=1:m, j=1:m)  ve(  j, l-i) = vm(i,j)
      end select
    else
      select case(face % rotation)
      case(0_IXS)
        forall (i=1:m, j=1:m)  ve(  i, l-j) = vm(i,j)
      case(1_IXS)
        forall (i=1:m, j=1:m)  ve(l-j, l-i) = vm(i,j)
      case(2_IXS)
        forall (i=1:m, j=1:m)  ve(l-i,   j) = vm(i,j)
      case default
        forall (i=1:m, j=1:m)  ve(  j,   i) = vm(i,j)
      end select
    end if

  end subroutine AlignFromMeshFace_IDK

  !-----------------------------------------------------------------------------
  !> Transforms face data from mesh to element orientation -- integer scalar

  pure subroutine AlignFromMeshFace_IXS(face, vm, ve)
    class(MeshElementFace_3D), intent(in) :: face  !< mesh element face
    integer(IXS), intent(in)  :: vm(:,:)  !< mesh face data
    integer(IXS), intent(out) :: ve(:,:)  !< element face data

    integer :: i, j, l, m

    m = size(vm,1)
    l = m + 1

    if (face % normal == 1_IXS) then
      select case(face % rotation)
      case(0_IXS)
        ve = vm
      case(1_IXS)
        forall (i=1:m, j=1:m)  ve(l-j,   i) = vm(i,j)
      case(2_IXS)
        forall (i=1:m, j=1:m)  ve(l-i, l-j) = vm(i,j)
      case default
        forall (i=1:m, j=1:m)  ve(  j, l-i) = vm(i,j)
      end select
    else
      select case(face % rotation)
      case(0_IXS)
        forall (i=1:m, j=1:m)  ve(  i, l-j) = vm(i,j)
      case(1_IXS)
        forall (i=1:m, j=1:m)  ve(l-j, l-i) = vm(i,j)
      case(2_IXS)
        forall (i=1:m, j=1:m)  ve(l-i,   j) = vm(i,j)
      case default
        forall (i=1:m, j=1:m)  ve(  j,   i) = vm(i,j)
      end select
    end if

  end subroutine AlignFromMeshFace_IXS

  !-----------------------------------------------------------------------------
  !> Transforms face data from mesh to element orientation -- real(RNP) scalar

  pure subroutine AlignFromMeshFace_RNP(face, vm, ve)
    class(MeshElementFace_3D), intent(in) :: face  !< mesh element face
    real(RNP), intent(in)  :: vm(:,:)  !< mesh face data
    real(RNP), intent(out) :: ve(:,:)  !< element face data

    integer :: i, j, l, m

    m = size(vm,1)
    l = m + 1

    if (face % normal == 1_IXS) then
      select case(face % rotation)
      case(0_IXS)
        ve = vm
      case(1_IXS)
        forall (i=1:m, j=1:m)  ve(l-j,   i) = vm(i,j)
      case(2_IXS)
        forall (i=1:m, j=1:m)  ve(l-i, l-j) = vm(i,j)
      case default
        forall (i=1:m, j=1:m)  ve(  j, l-i) = vm(i,j)
      end select
    else
      select case(face % rotation)
      case(0_IXS)
        forall (i=1:m, j=1:m)  ve(  i, l-j) = vm(i,j)
      case(1_IXS)
        forall (i=1:m, j=1:m)  ve(l-j, l-i) = vm(i,j)
      case(2_IXS)
        forall (i=1:m, j=1:m)  ve(l-i,   j) = vm(i,j)
      case default
        forall (i=1:m, j=1:m)  ve(  j,   i) = vm(i,j)
      end select
    end if

  end subroutine AlignFromMeshFace_RNP

  !=============================================================================
  ! MeshElement_3D procedures

  !-----------------------------------------------------------------------------
  !> Transforms neighbor face data to element orientation

  pure subroutine AlignFromNeighborFace(element, f, n, vn, ve)
    class(MeshElement_3D), intent(in) :: element  !< mesh element
    integer,   intent(in)  :: f        !< element face ID
    integer,   intent(in)  :: n        !< element neighbor ID
    real(RNP), intent(in)  :: vn(:,:)  !< neighbor face data
    real(RNP), intent(out) :: ve(:,:)  !< element face data

    integer :: i, j, l, m

    m = size(vn,1)
    l = m + 1

    select case(f)
    case(1,2)
      select case(element % neighbor(n) % orientation)
      case(12_IXS, 31_IXS, 51_IXS)
        ve = vn
      case(16_IXS, 35_IXS, 56_IXS)
        forall (i=1:m, j=1:m)  ve(l-j,   i) = vn(i,j)
      case(15_IXS, 34_IXS, 54_IXS)
        forall (i=1:m, j=1:m)  ve(l-i, l-j) = vn(i,j)
      case(13_IXS, 32_IXS, 53_IXS)
        forall (i=1:m, j=1:m)  ve(  j, l-i) = vn(i,j)
      case(42_IXS, 61_IXS, 21_IXS)
        forall (i=1:m, j=1:m)  ve(  i, l-j) = vn(i,j)
      case(46_IXS, 65_IXS, 26_IXS)
        forall (i=1:m, j=1:m)  ve(l-j, l-i) = vn(i,j)
      case(45_IXS, 64_IXS, 24_IXS)
        forall (i=1:m, j=1:m)  ve(l-i,   j) = vn(i,j)
      case default ! 43_IXS, 62_IXS, 23_IXS
        forall (i=1:m, j=1:m)  ve(  j,   i) = vn(i,j)
      end select
    case(3,4)
      select case(element % neighbor(n) % orientation)
      case(12_IXS, 16_IXS, 24_IXS)
        ve = vn
      case(62_IXS, 56_IXS, 64_IXS)
        forall (i=1:m, j=1:m)  ve(l-j,   i) = vn(i,j)
      case(42_IXS, 46_IXS, 54_IXS)
        forall (i=1:m, j=1:m)  ve(l-i, l-j) = vn(i,j)
      case(32_IXS, 26_IXS, 34_IXS)
        forall (i=1:m, j=1:m)  ve(  j, l-i) = vn(i,j)
      case(15_IXS, 13_IXS, 21_IXS)
        forall (i=1:m, j=1:m)  ve(  i, l-j) = vn(i,j)
      case(35_IXS, 23_IXS, 31_IXS)
        forall (i=1:m, j=1:m)  ve(l-j, l-i) = vn(i,j)
      case(45_IXS, 43_IXS, 51_IXS)
        forall (i=1:m, j=1:m)  ve(l-i,   j) = vn(i,j)
      case default ! 65_IXS, 53_IXS, 61_IXS
        forall (i=1:m, j=1:m)  ve(  j,   i) = vn(i,j)
      end select
    case default ! 5,6
      select case(element % neighbor(n) % orientation)
      case(12_IXS, 13_IXS, 23_IXS)
        ve = vn
      case(51_IXS, 61_IXS, 62_IXS)
        forall (i=1:m, j=1:m)  ve(l-j,   i) = vn(i,j)
      case(45_IXS, 46_IXS, 56_IXS)
        forall (i=1:m, j=1:m)  ve(l-i, l-j) = vn(i,j)
      case(24_IXS, 34_IXS, 35_IXS)
        forall (i=1:m, j=1:m)  ve(  j, l-i) = vn(i,j)
      case(15_IXS, 16_IXS, 26_IXS)
        forall (i=1:m, j=1:m)  ve(  i, l-j) = vn(i,j)
      case(21_IXS, 31_IXS, 32_IXS)
        forall (i=1:m, j=1:m)  ve(l-j, l-i) = vn(i,j)
      case(42_IXS, 43_IXS, 53_IXS)
        forall (i=1:m, j=1:m)  ve(l-i,   j) = vn(i,j)
      case default ! 54_IXS, 64_IXS, 65_IXS
        forall (i=1:m, j=1:m)  ve(  j,   i) = vn(i,j)
      end select
    end select

  end subroutine AlignFromNeighborFace

  !-----------------------------------------------------------------------------
  !> Transforms the given triple index range into a rotated coordinate system

  pure subroutine TransformIndex( orientation, n, o0, o, t0, t, s)

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

  !=============================================================================

end module Mesh_Element__3D
