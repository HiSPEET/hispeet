!> summary:  Creation of a 3x3x3 Cartesian mesh with the center element rotated
!> author:   Joerg Stiller
!> date:     2021/08/29
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!==============================================================================

submodule(Generic_Mesh__3D) MP_CreateOneRotated
  use Constants
  use Gauss_Jacobi
  use Execution_Control
  use Structured_Mesh_Indexing__3D
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Creates a 3x3x3 Cartesian mesh with the center element rotated
  !>
  !> After creating a structured mesh, the center element is rotated
  !> `rotation(i)` times about axes `i=1,2,3`, one after the other.

  module subroutine CreateOneRotated(mesh, xo, lx, po, rotation, periodic)
    class(GenericMesh_3D), intent(out) :: mesh  !< unstructured regular mesh
    real(RNP), intent(in) :: xo(3)       !< corner closest to -infinity
    real(RNP), intent(in) :: lx(3)       !< domain extensions
    integer,   intent(in) :: po          !< polynomial order of mesh elements
    integer,   intent(in) :: rotation(3) !< rotation applied to center element
    logical,   intent(in) :: periodic(3) !< F/T for non/periodic directions

    ! local variables
    integer   :: nv, ne, nb       ! number of vertices, elements and boundaries
    integer   :: i, j, k, p, q, r ! loop indices
    integer   :: l                ! linear element index
    integer   :: v                ! linear vertex index
    integer   :: f                ! face counter
    integer   :: b                ! boundary ID
    real(RNP) :: dx(3)            ! element extensions
    real(RNP) :: xc(0:po)         ! Lobatto points in [0,1]
    real(RNP) :: x1(3)            ! position vector

    ! prerequisites ...........................................................

    nv = 4 ** 3
    ne = 3 ** 3
    nb = 6

    dx = THIRD * lx

    ! Lobatto points in [0,1]
    xc = HALF * (ONE + LobattoPoints(po))

    ! allocate mesh components
    allocate(mesh%vertex(nv), mesh%element(ne), mesh%boundary(nb))

    ! set lexical numbering
    mesh%numbering = LEXICAL_NUMBERING

    ! vertices .................................................................

    do k = 0, 3
    do j = 0, 3
    do i = 0, 3
      v = LexicalVertexIndex(i, j, k, 3, 3)
      mesh % vertex(v) % id = LexicalVertexIndex(i, j, k, 3, 3, 3, periodic)
      mesh % vertex(v) % x  = xo + dx * [ i, j, k ]
    end do
    end do
    end do

    ! elements .................................................................

    mesh%element%typ   = HEXAHEDRAL_ELEMENT
    mesh%element%basis = GAUSS_LOBATTO_BASIS
    mesh%element%order = po

    do k = 1, 3
    do j = 1, 3
    do i = 1, 3
      l = LexicalElementIndex(i, j, k, 3, 3)

      mesh%element(l)%id = l

      ! vertices
      mesh%element(l)%vertex = [ LexicalVertexIndex(i-1, j-1, k-1, 3, 3), &
                                 LexicalVertexIndex(i  , j-1, k-1, 3, 3), &
                                 LexicalVertexIndex(i-1, j  , k-1, 3, 3), &
                                 LexicalVertexIndex(i  , j  , k-1, 3, 3), &
                                 LexicalVertexIndex(i-1, j-1, k  , 3, 3), &
                                 LexicalVertexIndex(i  , j-1, k  , 3, 3), &
                                 LexicalVertexIndex(i-1, j  , k  , 3, 3), &
                                 LexicalVertexIndex(i  , j  , k  , 3, 3)  ]

      ! element collocation points
      allocate(mesh%element(l)%x((po+1)**3, 3))
      x1 = mesh%vertex(mesh%element(l)%vertex(1)) % x ! coordinates of vertex 1
      do r = 0, po
      do q = 0, po
      do p = 0, po
        v = 1 + p + (po+1) * (q + (po+1) * r) ! linear point index
        mesh%element(l) % x(v,1) = x1(1) + dx(1) * xc(p)
        mesh%element(l) % x(v,2) = x1(2) + dx(2) * xc(q)
        mesh%element(l) % x(v,3) = x1(3) + dx(3) * xc(r)
      end do
      end do
      end do

      ! rotate center element
      if (i == 2 .and. j == 2 .and. k == 2) then
        call RotateElement(mesh%element(l), rotation)
      end if

    end do
    end do
    end do

    ! boundaries ...............................................................

    b = 1
    mesh%boundary(b)%id = b
    mesh%boundary(b)%name = 'west'
    allocate(mesh % boundary(b) % face(9))
    f = 0
    i = 1
    do k = 1, 3
    do j = 1, 3
      f = f + 1
      mesh%boundary(b)%face(f)%element_id   = LexicalElementIndex(i, j, k, 3, 3)
      mesh%boundary(b)%face(f)%element_face = 1
    end do
    end do

    b = b + 1
    mesh%boundary(b)%id = b
    mesh%boundary(b)%name = 'east'
    allocate(mesh % boundary(b) % face(9))
    f = 0
    i = 3
    do k = 1, 3
    do j = 1, 3
      f = f + 1
      mesh%boundary(b)%face(f)%element_id   = LexicalElementIndex(i, j, k, 3, 3)
      mesh%boundary(b)%face(f)%element_face = 2
    end do
    end do

    b = b + 1
    mesh%boundary(b)%id = b
    mesh%boundary(b)%name = 'south'
    allocate(mesh % boundary(b) % face(9))
    f = 0
    j = 1
    do k = 1, 3
    do i = 1, 3
      f = f + 1
      mesh%boundary(b)%face(f)%element_id   = LexicalElementIndex(i, j, k, 3, 3)
      mesh%boundary(b)%face(f)%element_face = 3
    end do
    end do

    b = b + 1
    mesh%boundary(b)%id = b
    mesh%boundary(b)%name = 'north'
    allocate(mesh % boundary(b) % face(9))
    f = 0
    j = 3
    do k = 1, 3
    do i = 1, 3
      f = f + 1
      mesh%boundary(b)%face(f)%element_id   = LexicalElementIndex(i, j, k, 3, 3)
      mesh%boundary(b)%face(f)%element_face = 4
    end do
    end do

    b = b + 1
    mesh%boundary(b)%id = b
    mesh%boundary(b)%name = 'bottom'
    allocate(mesh % boundary(b) % face(9))
    f = 0
    k = 1
    do j = 1, 3
    do i = 1, 3
      f = f + 1
      mesh%boundary(b)%face(f)%element_id   = LexicalElementIndex(i, j, k, 3, 3)
      mesh%boundary(b)%face(f)%element_face = 5
    end do
    end do

    b = b + 1
    mesh%boundary(b)%id = b
    mesh%boundary(b)%name = 'top'
    allocate(mesh % boundary(b) % face(9))
    f = 0
    k = 3
    do j = 1, 3
    do i = 1, 3
      f = f + 1
      mesh%boundary(b)%face(f)%element_id   = LexicalElementIndex(i, j, k, 3, 3)
      mesh%boundary(b)%face(f)%element_face = 6
    end do
    end do

    ! identify coupled boundaries and their polarity
    do i = 1, 3
      if (periodic(i)) then
        j = 2 * i - 1
        k = 2 * i
        mesh%boundary(j) % coupled  =  k
        mesh%boundary(j) % polarity = -i
        mesh%boundary(j) % map(i,4) =  lx(i)
        mesh%boundary(k) % coupled  =  j
        mesh%boundary(k) % polarity =  i
        mesh%boundary(k) % map(i,4) = -lx(i)
      end if
    end do

  end subroutine CreateOneRotated

  !-----------------------------------------------------------------------------
  !> Rotates the given element about the ξ, η and ζ-axes

  subroutine RotateElement(element, rotation)
    type(GenericMeshElement_3D), intent(inout) :: element
    integer, intent(in) :: rotation(3) !< num rotations about ξ/η/ζ-axes

    real(RNP) :: x(0:element%order,0:element%order,0:element%order,3)
    integer   :: np, po
    integer   :: i, j, k, r

    po = element % order
    np = po + 1

    ! rotation about ξ-axis
    do r = 1, rotation(1)
      ! vertices
      element % vertex = element % vertex([5,6,1,2,7,8,3,4])
      ! coordinates
      x = reshape(element % x, shape(x))
      do k = 0, po
      do j = 0, po
      do i = 0, po
        element % x(1 + i + np*(j + np*k),1:3) = x(i,k,po-j,1:3)
      end do
      end do
      end do
    end do

    ! rotation about η-axis
    do r = 1, rotation(2)
      ! vertices
      element % vertex = element % vertex([2,6,4,8,1,5,3,7])
      ! coordinates
      x = reshape(element % x, shape(x))
      do k = 0, po
      do j = 0, po
      do i = 0, po
        element % x(1 + i + np*(j + np*k),1:3) = x(po-k,j,i,1:3)
      end do
      end do
      end do
    end do

    ! rotation about ζ-axis
    do r = 1, rotation(3)
      ! vertices
      element % vertex = element % vertex([3,1,4,2,7,5,8,6])
      ! coordinates
      x = reshape(element % x, shape(x))
      do k = 0, po
      do j = 0, po
      do i = 0, po
        element % x(1 + i + np*(j + np*k),1:3) = x(j,po-i,k,1:3)
      end do
      end do
      end do
    end do

  end subroutine RotateElement

  !=============================================================================

end submodule MP_CreateOneRotated
