!> summary:  Creation of a 3d generic mesh build from non-Cartesian unit cells
!> author:   Joerg Stiller
!> date:     2021/08/21
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!==============================================================================

submodule(Generic_Mesh__3D) MP_CreateDiamonds
  use Constants
  use Gauss_Jacobi
  use Execution_Control
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Creates a generic mesh of cells, each containing a rhombic "diamond"
  !> element surrounded by four trapezoidal elements.
  !>
  !> Top view of a unit cell comprising 5 elements.
  !>          ____________
  !>         |     /\     |
  !>         | 2  /  \  4 |
  !>         |   /    \   |
  !>         |--⟨   5  ⟩--|
  !>         |   \    /   |
  !>         | 1  \  /  3 |
  !>         |_____\/_____|
  !>  x ↑
  !>     → y

  module subroutine CreateDiamonds(mesh, lx, ly, lz, nx, ny, nz, po, periodic)
    class(GenericMesh_3D), intent(out) :: mesh  !< "diamond" mesh
    real(RNP), intent(in) :: lx  !< domain length in x-direction
    real(RNP), intent(in) :: ly  !< domain length in y-direction
    real(RNP), intent(in) :: lz  !< domain length in z-direction
    integer,   intent(in) :: nx  !< number of cells in x-direction
    integer,   intent(in) :: ny  !< number of cells in y-direction
    integer,   intent(in) :: nz  !< number of elements in z-direction
    integer,   intent(in) :: po  !< polynomial order of mesh elements
    logical,   intent(in) :: periodic !< F/T for non/periodic z-direction

    logical,   allocatable :: vm(:,:,:,:) ! mask indicating orinal vertices
    integer,   allocatable :: vc(:,:,:,:) ! cell-vertex indices

    integer   :: nv        ! number of vertices
    integer   :: ne        ! number of elements
    integer   :: nb        ! number of boundaries
    real(RNP) :: dx(3)     ! cell dimensions

    integer :: b, e, f, i, j, k, l, v

    ! prerequisites ...........................................................

    ! number of elements and boundaries
    ne = 5 * nx * ny * nz
    nb = 6

    ! cell dimensions
    dx = [ lx/nx, ly/ny, lz/nz ]

    ! set lexical numbering
    mesh%numbering = LEXICAL_NUMBERING

    ! vertices .................................................................

    allocate(vc(1:10, 1:nx, 1:ny, 0:nz))
    allocate(vm(1:10, 1:nx, 1:ny, 0:nz))

    call IdentifyVertices(periodic, vm, vc, nv)

    allocate(mesh % vertex(nv))
    do k = 0, nz
    do j = 1, ny
    do i = 1, nx
    do l = 1, 10
      if (vm(l,i,j,k)) then
        v = vc(l,i,j,k)
        mesh % vertex(v) % id = v
        mesh % vertex(v) % x  = VertexCoordinates(l, i, j, k, dx)
      end if
    end do
    end do
    end do
    end do

    ! elements ...................................................................

    allocate(mesh % element(ne))

    mesh%element%typ   = HEXAHEDRAL_ELEMENT
    mesh%element%basis = GAUSS_LOBATTO_BASIS
    mesh%element%order = po

    do k = 1, nz
    do j = 1, ny
    do i = 1, nx
    do l = 1, 5

      e = ElementIndex(l, i, j, k, nx, ny)

      mesh % element(e) % id     = e
      mesh % element(e) % vertex = ElementVertexIDs(l, i, j, k, vc)
      mesh % element(e) % x      = ElementPoints(l, i, j, k, po, dx)

    end do
    end do
    end do
    end do

    ! boundaries ...............................................................

    allocate(mesh % boundary(nb))

    b = 1
    mesh%boundary(b)%id = b
    mesh%boundary(b)%name = 'west'
    allocate(mesh % boundary(b) % face(2*ny*nz))
    f = 0
    i = 1
    do k = 1, nz
    do j = 1, ny
      f = f + 1
      mesh%boundary(b)%face(f)%element_id   = ElementIndex(1, i, j, k, nx, ny)
      mesh%boundary(b)%face(f)%element_face = 1
      f = f + 1
      mesh%boundary(b)%face(f)%element_id   = ElementIndex(3, i, j, k, nx, ny)
      mesh%boundary(b)%face(f)%element_face = 1
    end do
    end do

    b = b + 1
    mesh%boundary(b)%id = b
    mesh%boundary(b)%name = 'east'
    allocate(mesh % boundary(b) % face(2*ny*nz))
    f = 0
    i = nx
    do k = 1, nz
    do j = 1, ny
      f = f + 1
      mesh%boundary(b)%face(f)%element_id   = ElementIndex(2, i, j, k, nx, ny)
      mesh%boundary(b)%face(f)%element_face = 2
      f = f + 1
      mesh%boundary(b)%face(f)%element_id   = ElementIndex(4, i, j, k, nx, ny)
      mesh%boundary(b)%face(f)%element_face = 2
    end do
    end do

    b = b + 1
    mesh%boundary(b)%id = b
    mesh%boundary(b)%name = 'south'
    allocate(mesh % boundary(b) % face(2*nx*nz))
    f = 0
    j = 1
    do k = 1, nz
    do i = 1, nx
      f = f + 1
      mesh%boundary(b)%face(f)%element_id   = ElementIndex(1, i, j, k, nx, ny)
      mesh%boundary(b)%face(f)%element_face = 3
      f = f + 1
      mesh%boundary(b)%face(f)%element_id   = ElementIndex(2, i, j, k, nx, ny)
      mesh%boundary(b)%face(f)%element_face = 3
    end do
    end do

    b = b + 1
    mesh%boundary(b)%id = b
    mesh%boundary(b)%name = 'north'
    allocate(mesh % boundary(b) % face(2*nx*nz))
    f = 0
    j = ny
    do k = 1, nz
    do i = 1, nx
      f = f + 1
      mesh%boundary(b)%face(f)%element_id   = ElementIndex(3, i, j, k, nx, ny)
      mesh%boundary(b)%face(f)%element_face = 4
      f = f + 1
      mesh%boundary(b)%face(f)%element_id   = ElementIndex(4, i, j, k, nx, ny)
      mesh%boundary(b)%face(f)%element_face = 4
    end do
    end do

    b = b + 1
    mesh%boundary(b)%id = b
    mesh%boundary(b)%name = 'bottom'
    allocate(mesh % boundary(b) % face(5*nx*ny))
    f = 0
    k = 1
    do j = 1, ny
    do i = 1, nx
      e = ElementIndex(1, i, j, k, nx, ny) - 1
      do l = 1, 5
        mesh%boundary(b)%face(f+l)%element_id   = e + l
        mesh%boundary(b)%face(f+l)%element_face = 5
      end do
      f = f + 5
    end do
    end do

    b = b + 1
    mesh%boundary(b)%id = b
    mesh%boundary(b)%name = 'top'
    allocate(mesh % boundary(b) % face(5*nx*ny))
    f = 0
    k = nz
    do j = 1, ny
    do i = 1, nx
      e = ElementIndex(1, i, j, k, nx, ny) - 1
      do l = 1, 5
        mesh%boundary(b)%face(f+l)%element_id   = e + l
        mesh%boundary(b)%face(f+l)%element_face = 6
      end do
      f = f + 5
    end do
    end do

    ! identify coupled boundaries and their polarity
    if (periodic) then
      mesh%boundary(5) % coupled  =  6
      mesh%boundary(5) % polarity = -3
      mesh%boundary(6) % coupled  =  5
      mesh%boundary(6) % polarity =  3
    end if

  end subroutine CreateDiamonds

  !-----------------------------------------------------------------------------
  !> Identification and numbering of cell vertices
  !>
  !> The vertex indices are stored in the array `vc(1:10,1:nx,1:ny,0:nz)`,
  !> where `nx`, `ny` are the number of cells in x- and y-direction and `nz`
  !> the number layers in z-direction. On output `vc(l,i,j,k)` contains the
  !> index of vertex `l` in cell `i,j` and z-plane `k`. The vertex mask `vm`
  !> indicates whether the vertex is an original (T) or a copy (F). Copies
  !> adopt their index from the corresponding original vertex.

  subroutine IdentifyVertices(periodic, vm, vc, nv)
    logical, intent(in)  :: periodic        !< F/T for non/periodic z-direction
    logical, intent(out) :: vm(1:,1:,1:,0:) !< T/F for original/copy vertices
    integer, intent(out) :: vc(1:,1:,1:,0:) !< global cell vertex IDs (l,i,j,k)
    integer, intent(out) :: nv              !< number of original vertices

    integer :: nx, ny, nz
    integer :: i, j, k, l

    ! bounds
    nx = ubound(vc,2)
    ny = ubound(vc,3)
    nz = ubound(vc,4)

    vm = .true.

    ! unmask copies of shared vertices on x-faces of cells
    vm(1,2:,:,:) = .false.
    vm(4,2:,:,:) = .false.
    vm(8,2:,:,:) = .false.

    ! unmask copies of shared vertices on y-faces of cells
    vm(1,:,2:,:) = .false.
    vm(2,:,2:,:) = .false.
    vm(3,:,2:,:) = .false.

    ! unmask z-periodic vertices
    if (periodic) then
      vm(:,:,:,nz) = .false.
    end if

    ! initialize vertex IDs
    vc = -1

    ! number original vertices
    nv = 0
    do k = 0, nz
    do j = 1, ny
    do i = 1, nx
    do l = 1, 10
      if (vm(l,i,j,k)) then
        nv = nv + 1
        vc(l,i,j,k) = nv
      end if
    end do
    end do
    end do
    end do

    ! identify copies of shared vertices on x-faces of cells
    vc(1,2:,:,:) = vc( 3,:nx-1,:,:)
    vc(4,2:,:,:) = vc( 7,:nx-1,:,:)
    vc(8,2:,:,:) = vc(10,:nx-1,:,:)

    ! identify copies of shared vertices on y-faces of cells
    vc(1,:,2:,:) = vc( 8,:,:ny-1,:)
    vc(2,:,2:,:) = vc( 9,:,:ny-1,:)
    vc(3,:,2:,:) = vc(10,:,:ny-1,:)

    ! identify z-periodic vertices
    if (periodic) then
      vc(:,:,:,nz) = vc(:,:,:,0)
    end if

    ! check
    if (any(vc < 0)) then
      call Error('IdentifyVertices', 'failed', &
                 'Generic_Mesh__3D:MP_CreateDiamonds')
    end if

  end subroutine IdentifyVertices

  !-----------------------------------------------------------------------------
  !> Returns the Cartesian coordinates of a vertex

  pure function VertexCoordinates(l, i, j, k, dx)  result(x)
    integer,   intent(in) :: l          !< vertex ID within xy-cell
    integer,   intent(in) :: i, j       !< xy-cell index
    integer,   intent(in) :: k          !< z-layer index
    real(RNP), intent(in) :: dx(3)      !< cell extensions in x,y,z-directions
    real(RNP) :: x(3)

    real(RNP) :: xo(3)

    xo = dx * [ i-1, j-1, k ]

    select case(l)
    case(1)
      x = xo + dx/4 * [ 0, 0, 0 ]
    case(2)
      x = xo + dx/4 * [ 2, 0, 0 ]
    case(3)
      x = xo + dx/4 * [ 4, 0, 0 ]
    case(4)
      x = xo + dx/4 * [ 0, 2, 0 ]
    case(5)
      x = xo + dx/4 * [ 2, 1, 0 ]
    case(6)
      x = xo + dx/4 * [ 2, 3, 0 ]
    case(7)
      x = xo + dx/4 * [ 4, 2, 0 ]
    case(8)
      x = xo + dx/4 * [ 0, 4, 0 ]
    case(9)
      x = xo + dx/4 * [ 2, 4, 0 ]
    case default ! 10
      x = xo + dx/4 * [ 4, 4, 0 ]
    end select

  end function VertexCoordinates

  !-----------------------------------------------------------------------------
  !> Returns the index of element `l` in cell `i,j,k`

  pure function ElementIndex(l, i, j, k, nx, ny) result(e)
    integer, intent(in) :: l       !< element ID within cell
    integer, intent(in) :: i, j, k !< triple cell ID
    integer, intent(in) :: nx      !< number of cells in x-direction
    integer, intent(in) :: ny      !< number of cells in y-direction
    integer :: e

    e = l + 5 * (i-1 + nx * (j-1 + ny * (k-1)))

  end function ElementIndex

  !-----------------------------------------------------------------------------
  !> Returns the vertex IDs of element `l` in cell `i,j,k`

  pure function ElementVertexIDs(l, i, j, k, vc) result(ve)
    integer, intent(in) :: l               !< element ID within cell
    integer, intent(in) :: i, j, k         !< triple cell ID
    integer, intent(in) :: vc(1:,1:,1:,0:) !< global cell vertex IDs (l,i,j,k)
    integer :: ve(8)

    select case(l)
    case(1)
      ve = [  vc( 1,i,j,k-1 ), vc( 2,i,j,k-1 ), &
              vc( 4,i,j,k-1 ), vc( 5,i,j,k-1 ), &
              vc( 1,i,j,k   ), vc( 2,i,j,k   ), &
              vc( 4,i,j,k   ), vc( 5,i,j,k   )  ]
    case(2)
      ve = [  vc( 2,i,j,k-1 ), vc( 3,i,j,k-1 ), &
              vc( 5,i,j,k-1 ), vc( 7,i,j,k-1 ), &
              vc( 2,i,j,k   ), vc( 3,i,j,k   ), &
              vc( 5,i,j,k   ), vc( 7,i,j,k   )  ]
    case(3)
      ve = [  vc( 4,i,j,k-1 ), vc( 6,i,j,k-1 ), &
              vc( 8,i,j,k-1 ), vc( 9,i,j,k-1 ), &
              vc( 4,i,j,k   ), vc( 6,i,j,k   ), &
              vc( 8,i,j,k   ), vc( 9,i,j,k   )  ]
    case(4)
      ve = [  vc( 6,i,j,k-1 ), vc( 7,i,j,k-1 ), &
              vc( 9,i,j,k-1 ), vc(10,i,j,k-1 ), &
              vc( 6,i,j,k   ), vc( 7,i,j,k   ), &
              vc( 9,i,j,k   ), vc(10,i,j,k   )  ]
    case default ! 5
      ve = [  vc( 4,i,j,k-1 ), vc( 5,i,j,k-1 ), &
              vc( 6,i,j,k-1 ), vc( 7,i,j,k-1 ), &
              vc( 4,i,j,k   ), vc( 5,i,j,k   ), &
              vc( 6,i,j,k   ), vc( 7,i,j,k   )  ]
    end select

  end function ElementVertexIDs

  !-----------------------------------------------------------------------------
  !> Returns the Lobatto points of element `l` in cell `i,j,k`

  pure function ElementPoints(l, i, j, k, po, dx) result(x)
    integer,   intent(in) :: l        !< element ID within cell
    integer,   intent(in) :: i, j, k  !< triple cell ID
    integer,   intent(in) :: po       !< polynomial order of mesh element
    real(RNP), intent(in) :: dx(3)    !< cell extensions in x,y,z-directions
    real(RNP) :: x((po+1)**3, 3)

    real(RNP), dimension(3) :: x000, x100, x010, x110, x001, x101, x011, x111
    real(RNP), dimension(3) :: x00, x10, x01, x11, x0, x1
    real(RNP) :: xs(0:po)
    integer   :: m(4)
    integer   :: r, p, q, v

    ! Lobatto points in [0,1]
    xs = HALF * (ONE + LobattoPoints(po))

    select case(l)
    case(1)
      m = [1,2,4,5]
    case(2)
      m = [2,3,5,7]
    case(3)
      m = [4,6,8,9]
    case(4)
      m = [6,7,9,10]
    case default ! 5
      m = [4,5,6,7]
    end select

    ! element vertex coordinates
    x000 = VertexCoordinates( m(1), i, j, k-1, dx )
    x100 = VertexCoordinates( m(2), i, j, k-1, dx )
    x010 = VertexCoordinates( m(3), i, j, k-1, dx )
    x110 = VertexCoordinates( m(4), i, j, k-1, dx )
    x001 = VertexCoordinates( m(1), i, j, k  , dx )
    x101 = VertexCoordinates( m(2), i, j, k  , dx )
    x011 = VertexCoordinates( m(3), i, j, k  , dx )
    x111 = VertexCoordinates( m(4), i, j, k  , dx )

    ! interpolation to collocation points
    do r = 0, po
    do q = 0, po
    do p = 0, po

      ! linear point index
      v = 1 + p + (po+1) * (q + (po+1) * r)

      ! trilinear interpolation: direction 1
      x00  =  x000 * (1 - xs(p))  +  x100 * xs(p)
      x01  =  x001 * (1 - xs(p))  +  x101 * xs(p)
      x10  =  x010 * (1 - xs(p))  +  x110 * xs(p)
      x11  =  x011 * (1 - xs(p))  +  x111 * xs(p)

      ! trilinear interpolation: direction 2
      x0  =  x00 * (1 - xs(q))  +  x10 * xs(q)
      x1  =  x01 * (1 - xs(q))  +  x11 * xs(q)

      ! trilinear interpolation: direction 3
      x(v,:)  =  x0 * (1 - xs(r))  +  x1 * xs(r)

    end do
    end do
    end do

  end function ElementPoints

  !============================================================================

end submodule MP_CreateDiamonds
