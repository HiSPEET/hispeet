!> summary:  Generic tree-dimensional mesh
!> author:   Joerg Stiller
!> date:     2014/08/05, revised 2020/12/06
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!==============================================================================

module Generic_Mesh__3D
  use Kind_Parameters
  use Mesh_Element_Indexing__3D, only: V_FACE
  implicit none
  private

  !-----------------------------------------------------------------------------
  ! public constants

  ! numbering schemes for hexahedra
  public :: ROTATIONAL_NUMBERING
  public :: LEXICAL_NUMBERING

  ! element types
  public :: TETRAHEDRAL_ELEMENT
  public :: PRISMATIC_ELEMENT
  public :: PYRAMIDIC_ELEMENT
  public :: HEXAHEDRAL_ELEMENT

  ! basis functions for represention the element geometry
  public :: BEZIER_BERNSTEIN_BASIS
  public :: GAUSS_LOBATTO_BASIS

  !-----------------------------------------------------------------------------
  !> Structure defining a generic 3d mesh vertex

  type, public :: GenericMeshVertex_3D
    integer   :: id = 0 !< identifer
    real(RNP) :: x(3)   !< Cartesian coordinates
  end type GenericMeshVertex_3D

  !-----------------------------------------------------------------------------
  !> Structure defining a generic 3d mesh element
  !>
  !> The following element types are defined:
  !>
  !>   element type | specifier           | #vertices
  !>   -------------|---------------------|----------
  !>   tetrahedron  | TETRAHEDRAL_ELEMENT |     4
  !>   pyramid      | PYRAMIDIC_ELEMENT   |     5
  !>   wedge        | PRISMATIC_ELEMENT   |     6
  !>   hexahedron   | HEXAHEDRAL_ELEMENT  |     8
  !>
  !> The pyramid has a quadrilateral base, the wedge is a triangular prism.
  !>
  !> The leading entries of vertex comprise the IDs of the corresponding
  !> vertices defined in GenericMesh_3D, whereas the unused entries keep
  !> their initial value of 0.
  !>
  !> Wedges, pyramids and tetrahedrons are numbered as follows:
  !>
  !>               6-----5                                           4
  !>               |\   /|                    5                     /|\
  !>               | \ / |                  ./ \.                  / | \
  !>               |  4  |                .´/   \`.               /  |  \
  !>               |  |  |              .´ /     \ `.            3---|---2
  !>               3--|--2             4--/-------\--3            \  |  /
  !>                \ | /              | /         \ |             \ | /
  !>                 \|/               |/           \|              \|/
  !>                  1                1-------------2               1
  !>
  !> Note that numbering starts at 1.
  !>
  !> Hexahedra allow for two different vertex numbering schemes:
  !> lexical numbering (left) and rotational numbering (right)
  !>
  !>               7-------8             8-------7
  !>               |\      |\            |\      |\
  !>               | 5-------6           | 5-------6
  !>               | |     | |           | |     | |
  !>               3-|-----4 |           4-|-----3 |
  !>                \|      \|            \|      \|
  !>                 1-------2             1-------2
  !>
  !> The scheme is set globally and applied to all hexahedra within the mesh.
  !>
  !> With lexical numbering, the faces are associated with the element
  !> coordinates xi, eta, zeta:
  !>   * hexahedron, lexical numbering
  !>     - face 1  <-->  vertices 1, 3, 5, 7  <-->  xi   = -1
  !>     - face 2  <-->  vertices 2, 4, 6, 8  <-->  xi   =  1
  !>     - face 3  <-->  vertices 1, 2, 5, 6  <-->  eta  = -1
  !>     - face 4  <-->  vertices 3, 4, 7, 8  <-->  eta  =  1
  !>     - face 5  <-->  vertices 1, 2, 3, 4  <-->  zeta = -1
  !>     - face 6  <-->  vertices 5, 6, 7, 8  <-->  zeta =  1
  !>
  !>     (encoded in V_FACE adopted from Mesh_Element_Indexing__3D)
  !>
  !> In rotational, mode face nodes are numbered about the external normal:
  !>
  !>   * hexahedron
  !>     - face 1  <-->  vertices  1, 2, 6, 5
  !>     - face 2  <-->  vertices  2, 3, 7, 6
  !>     - face 3  <-->  vertices  3, 4, 8, 7
  !>     - face 4  <-->  vertices  4, 1, 5, 8
  !>     - face 5  <-->  vertices  1, 4, 3, 2
  !>     - face 6  <-->  vertices  5, 6, 7, 8
  !>
  !>   * wedge
  !>     - face 1  <-->  vertices  1, 2, 5, 4
  !>     - face 2  <-->  vertices  2, 3, 6, 5
  !>     - face 3  <-->  vertices  3, 1, 4, 6
  !>     - face 4  <-->  vertices  1, 3, 2
  !>     - face 5  <-->  vertices  4, 5, 6
  !>
  !>   * pyramid
  !>     - face 1  <-->  vertices  1, 2, 5
  !>     - face 2  <-->  vertices  2, 3, 5
  !>     - face 3  <-->  vertices  3, 4, 5
  !>     - face 4  <-->  vertices  4, 1, 5
  !>     - face 5  <-->  vertices  1, 4, 3, 2
  !>
  !>   * tetrahedron
  !>     - face 1  <-->  vertices  1, 2, 4
  !>     - face 2  <-->  vertices  2, 3, 4
  !>     - face 3  <-->  vertices  3, 1, 4
  !>     - face 4  <-->  vertices  1, 3, 2
  !>
  !> Optionally, the element geometry can be defined by specifying the basis
  !> functions, polynomial order and point coordinates.
  !>
  !> Supported basis functions:
  !>   * BEZIER_BERNSTEIN:  Bezier-form using Bernstein polynomials
  !>   * GAUSS_LOBATTO:     Lobatto-based Lagrange polynomials (hexahedra only)
  !>
  !> Point numbering
  !>   * hexahedron: lexical numbering based on triple index

  type, public :: GenericMeshElement_3D
    integer :: id        = 0          !< identifer
    integer :: typ       = 0          !< element type
    integer :: vertex(8) = 0          !< vertex indices (not IDs!)
    integer :: basis     = 0          !< type of basis functions
    integer :: order     = 0          !< polynomial order
    real(RNP), allocatable :: x(:,:)  !< control/collocation points, x(:,1:3)
  end type GenericMeshElement_3D

  !-----------------------------------------------------------------------------
  !> Structure for accessing a generic 3d mesh element face

  type, public :: GenericMeshElementFace_3D
    integer :: element_id    !< element ID
    integer :: element_face  !< element face
  end type GenericMeshElementFace_3D

  !-----------------------------------------------------------------------------
  !> Structure defining a generic 3d mesh boundary

  type, public :: GenericMeshBoundary_3D
    integer           :: id   = 0     !< identifer
    character(len=80) :: name = ''    !< boundary name
    integer           :: coupled = 0  !< ID of coupled boundary, 0 if none
    integer           :: polarity = 0 !< position WRT to periodic direction
    type(GenericMeshElementFace_3D), allocatable :: face(:) !< element faces
  end type GenericMeshBoundary_3D

  !-----------------------------------------------------------------------------
  !> Structure defining a generic 3d mesh
  !>
  !> For minimizing the overhead we consider only vertices, elements and mesh
  !> boundaries. Numbering starts at 1 throughout.
  !>
  !> Periodicity is supported through one-by-one coupling: Corresponding element
  !> faces must appear in the same order and numbered consistently. Coupling of
  !> boundaries b and c implies that
  !>
  !>   - boundary(b)%polarity = -boundary(c)%polarity ≠ 0
  !>
  !>   - boundary(b)%face(f) matches boundary(c)%face(f), and
  !>
  !>   - the vertices of the corresponding faces of
  !>
  !>       *  element( boundary(b)%face(f)%element ) and
  !>       *  element( boundary(c)%face(f)%element )
  !>
  !>     are coupled to each other according to element vertex numbering.
  !>     For example, with lexically numbered hexadra:
  !>
  !>       *  if face 1 of element i is coupled with face 2 of element j,
  !>       *  then vertices 1,3,5,7 of element i are coupled with vertices
  !>          2,4,5,8 of element j in exactly this order.
  !>
  !> All boundaries residing on the same extremity of a periodic direction
  !> must have an identical polarity.
  !>
  !> Compliance to these conventions is the precondition for constructing
  !> mappings between coupled mesh vertices on a purely topological basis.

  type, public :: GenericMesh_3D
    ! attributes
    integer :: numbering = 0  !< hexahedron numbering scheme
    ! components
    type(GenericMeshVertex_3D),   allocatable :: vertex(:)    !< vertices
    type(GenericMeshElement_3D),  allocatable :: element(:)   !< elements
    type(GenericMeshBoundary_3D), allocatable :: boundary(:)  !< boundaries
  contains
    ! type-bound procedures
    procedure :: SwitchToLexicalNumbering
    procedure :: SwitchToRotationalNumbering
    procedure :: GenerateConsistentVertexIDs
    procedure :: CreateCylinder
    procedure :: CreateAnnularGap
  end type GenericMesh_3D

  interface

    !---------------------------------------------------------------------------
    !> Creates a generic mesh for a cylinder

    module subroutine CreateCylinder(mesh, nr, nz, po, periodic)
      class(GenericMesh_3D), intent(out) :: mesh  !< cylindrical 3d mesh
      integer, intent(in) :: nr  !< num intervals in radial section
      integer, intent(in) :: nz  !< num intervals in axial  direction
      integer, intent(in) :: po  !< polynomial order of mesh elements
      logical, intent(in) :: periodic !< switch for axial periodicity
    end subroutine CreateCylinder

    !---------------------------------------------------------------------------
    !> Creates a generic mesh for an annular gap

    module subroutine CreateAnnularGap(mesh, r0, r1, h, nr, np, nz, po, periodic)
      class(GenericMesh_3D), intent(out) :: mesh  !< cylindrical 3d mesh
      real(RNP), intent(in) :: r0  !< inner radius
      real(RNP), intent(in) :: r1  !< outer radius
      real(RNP), intent(in) :: h   !< height
      integer,   intent(in) :: nr  !< num elements in radial    (r)   direction
      integer,   intent(in) :: np  !< num elements in azimuthal (phi) direction
      integer,   intent(in) :: nz  !< num elements in axial     (z)   direction
      integer,   intent(in) :: po  !< polynomial order of mesh elements
      logical,   intent(in) :: periodic !< switch for axial periodicity
    end subroutine CreateAnnularGap

  end interface

  !=============================================================================
  ! Parameters

  !-----------------------------------------------------------------------------
  ! Numbering schemes for hexahedra

  integer, parameter :: ROTATIONAL_NUMBERING = 1 !< rotational numbering
  integer, parameter :: LEXICAL_NUMBERING    = 2 !< lexical numbering

  !-----------------------------------------------------------------------------
  ! Element types

  integer, parameter :: TETRAHEDRAL_ELEMENT = 1 !< tetrahedron specifier
  integer, parameter :: PRISMATIC_ELEMENT   = 2 !< prism (wedge) specifier
  integer, parameter :: PYRAMIDIC_ELEMENT   = 3 !< pyramid specifier
  integer, parameter :: HEXAHEDRAL_ELEMENT  = 4 !< hexahedron specifier

  !-----------------------------------------------------------------------------
  ! Basis functions for represention the element geometry

  integer, parameter :: BEZIER_BERNSTEIN_BASIS = 1 ! Bezier-Bernstein basis
  integer, parameter :: GAUSS_LOBATTO_BASIS    = 2 ! Lagrangian (hexahedra only)

contains

  !-----------------------------------------------------------------------------
  !> Transforms hexahedra to lexical numbering scheme

  subroutine SwitchToLexicalNumbering(mesh)
    class(GenericMesh_3D), intent(inout) :: mesh

    integer :: b, i, k

    if (mesh%numbering == LEXICAL_NUMBERING) return

    ! element numbering ........................................................

    do i = 1, size(mesh%element)

      if (mesh%element(i)%typ /= HEXAHEDRAL_ELEMENT) cycle

      k = mesh%element(i)%vertex(3)
      mesh%element(i)%vertex(3) = mesh%element(i)%vertex(4)
      mesh%element(i)%vertex(4) = k

      k = mesh%element(i)%vertex(7)
      mesh%element(i)%vertex(7) = mesh%element(i)%vertex(8)
      mesh%element(i)%vertex(8) = k

    end do

    ! boundary face numbering ..................................................

    do b = 1, size(mesh%boundary)
      do i = 1, size(mesh%boundary(b)%face)

        associate(face => mesh%boundary(b)%face(i))
          if (mesh%element( face%element_id )%typ /= HEXAHEDRAL_ELEMENT) cycle
          select case(face%element_face)
          case(1)
            face%element_face = 3
          case(3)
            face%element_face = 4
          case(4)
            face%element_face = 1
          end select
        end associate

      end do
    end do

    mesh%numbering = LEXICAL_NUMBERING

  end subroutine SwitchToLexicalNumbering

  !-----------------------------------------------------------------------------
  !> Transforms hexahedra to rotational numbering scheme

  subroutine SwitchToRotationalNumbering(mesh)
    class(GenericMesh_3D), intent(inout) :: mesh

    integer :: b, i, k

    if (mesh%numbering == ROTATIONAL_NUMBERING) return

    ! element numbering ........................................................

    do i = 1, size(mesh%element)

      if (mesh%element(i)%typ /= HEXAHEDRAL_ELEMENT) cycle

      k = mesh%element(i)%vertex(3)
      mesh%element(i)%vertex(3) = mesh%element(i)%vertex(4)
      mesh%element(i)%vertex(4) = k

      k = mesh%element(i)%vertex(7)
      mesh%element(i)%vertex(7) = mesh%element(i)%vertex(8)
      mesh%element(i)%vertex(8) = k

    end do

    ! boundary face numbering ..................................................

    do b = 1, size(mesh%boundary)
      do i = 1, size(mesh%boundary(b)%face)

        associate(face => mesh%boundary(b)%face(i))
          if (mesh%element( face%element_id )%typ /= HEXAHEDRAL_ELEMENT) cycle
          select case(face%element_face)
          case(1)
            face%element_face = 4
          case(3)
            face%element_face = 1
          case(4)
            face%element_face = 3
          end select
        end associate

      end do
    end do

    mesh%numbering = ROTATIONAL_NUMBERING

  end subroutine SwitchToRotationalNumbering

  !-----------------------------------------------------------------------------
  !> Generates consistent vertex IDs including unique identifiers for periodic
  !> vertices.

  subroutine GenerateConsistentVertexIDs(mesh)
    class(GenericMesh_3D), intent(inout) :: mesh

    integer :: i, k
    integer :: b, lb, fb
    integer :: c, lc, fc

    associate(vertex => mesh % vertex, element => mesh % element)

      ! initialization .........................................................

      call mesh % SwitchToLexicalNumbering()

      do i = 1, size(mesh % vertex)
        mesh % vertex(i) % id = i
      end do

      ! periodic boundaries ....................................................

      do b = 1, size(mesh%boundary)
        c = mesh % boundary(b) % coupled
        if (c > b) then
          do i = 1, size(mesh % boundary(b) % face)
            lb = mesh % boundary(b) % face(i) % element_id
            fb = mesh % boundary(b) % face(i) % element_face
            lc = mesh % boundary(c) % face(i) % element_id
            fc = mesh % boundary(c) % face(i) % element_face
            do k = 1, 4
              vertex( element(lc) % vertex(V_FACE(k,fc)) ) % id = &
              vertex( element(lb) % vertex(V_FACE(k,fb)) ) % id
            end do
          end do
        end if
      end do

    end associate

  end subroutine GenerateConsistentVertexIDs

  !=============================================================================

end module Generic_Mesh__3D
