!> summary:  Import of generic 3d meshes
!> author:   Joerg Stiller
!> date:     2020/12/25
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(Mesh_Partition__3D) MP_ImportGenericMesh
  use Execution_Control
  use Standard_Operators__1D
  use Embedded_Interpolation__1D
  use Generic_Mesh__3D
  use Mesh_Element_Indexing__3D, only: V_FACE, E_FACE
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Import a generic 3d mesh

  module subroutine ImportGenericMesh(mesh, generic_mesh, comm)
    class(MeshPartition_3D), intent(out) :: mesh         !< mesh partition
    class(GenericMesh_3D),   intent(in)  :: generic_mesh !< generic mesh
    type(MPI_Comm),          intent(in)  :: comm         !< MPI communicator

    integer :: i, rank

    call MPI_Comm_rank(comm, rank)
    if (rank > 0) return

    ! check prerequisites ......................................................

    if (generic_mesh % numbering /= LEXICAL_NUMBERING) then
      call Error('ImportGenericMesh', 'Lexical numbering is required')
    else if (any(generic_mesh % element%typ /= HEXAHEDRAL_ELEMENT)) then
      call Error('ImportGenericMesh', 'Only hexahedral elements supported')
    end if

    mesh % n_bound     =  size(generic_mesh % boundary)
    mesh % n_elem      =  size(generic_mesh % element)
    mesh % n_part      =  1
    mesh % part        =  0
    mesh % structured  = .false.
    mesh % regular     = .false.
    mesh % comm        =  comm

    allocate(mesh % boundary( mesh%n_bound ))
    do i = 1, mesh % n_bound
      mesh % boundary(i) = MeshBoundary_3D( i                                   &
                             , name     = generic_mesh % boundary(i) % name     &
                             , coupled  = generic_mesh % boundary(i) % coupled  &
                             , polarity = generic_mesh % boundary(i) % polarity )
    end do


    call ImportElements(mesh, generic_mesh)

    call mesh % IdentifyEdges()
    call mesh % BuildFaces()

    call ImportBoundaryFaces  (mesh, generic_mesh)
    call ImportElementDomains (mesh, generic_mesh)

    call mesh % BuildLinks()
    call mesh % BuildGhosts()
    call mesh % IdentifyRanks()

  end subroutine ImportGenericMesh

  !-----------------------------------------------------------------------------
  !> Import of mesh elements
  !>
  !> Allocates the elements, sets global and local IDs, and adopts the vertex
  !> IDs from then generic mesh.

  subroutine ImportElements(mesh, generic_mesh)
    class(MeshPartition_3D), intent(inout) :: mesh         !< mesh partition
    class(GenericMesh_3D),    intent(in)    :: generic_mesh !< generic mesh

    integer :: b, i, j, k, n

    allocate(mesh % element( mesh % n_elem ))

    ! element and vertex IDs ...................................................

    n = 0
    do i = 1, mesh % n_elem
      mesh % element(i) % global_id = i
      mesh % element(i) % local_id  = i
      do j = 1, 8
        k = generic_mesh % vertex( generic_mesh % element(i) % vertex(j) ) % id
        mesh % element(i) % vertex(j) % id = k
        n = max(k, n)
      end do
    end do

    mesh % n_vert = n

    ! boundary IDs .............................................................

    do b = 1, size(generic_mesh%boundary)
      associate(boundary => generic_mesh % boundary(b))
        do i = 1, size(boundary % face)

          j = boundary % face(i) % element_id    ! corresponding element ID
          k = boundary % face(i) % element_face  ! corresponding element face

          mesh % element(j) % face(k) % boundary = b

        end do
      end associate
    end do

  end subroutine ImportElements

  !-----------------------------------------------------------------------------
  !> Import of boundary faces

  subroutine ImportBoundaryFaces(mesh, generic_mesh)
    class(MeshPartition_3D), intent(inout) :: mesh         !< mesh partition
    class(GenericMesh_3D),    intent(in)    :: generic_mesh !< generic mesh

    integer :: b, e, f, i, j, s

    do b = 1, size(generic_mesh%boundary)

      associate( mb => mesh % boundary(b)         &
               , gb => generic_mesh % boundary(b) )

        allocate(mb % face(mb % n_face))

        do i = 1, mb % n_face

          e = gb % face(i) % element_id         ! corresponding element ID
          j = gb % face(i) % element_face       ! corresponding element face
          f = mesh % element(e) % face(j) % id  ! corresponding mesh face ID

          ! side of the mesh face on which the boundary is located
          if (mesh % face(f) % element(1) % id == e) then
            s = 2
          else
            s = 1
          end if

          mb % face(i) % mesh_face % id   = f
          mb % face(i) % mesh_face % side = s

          mb % face(i) % mesh_element % id   = e
          mb % face(i) % mesh_element % face = j

        end do

      end associate
    end do

  end subroutine ImportBoundaryFaces

  !-----------------------------------------------------------------------------
  !> Import of mesh boundaries

  subroutine ImportBoundaries(mesh, generic_mesh)
    class(MeshPartition_3D), intent(inout) :: mesh         !< mesh partition
    class(GenericMesh_3D),    intent(in)    :: generic_mesh !< generic mesh

    integer :: b, e, f, i, j, s

    ! prerequisites ..............................................................

    mesh % n_bound = size(generic_mesh % boundary)

    allocate(mesh % boundary(mesh%n_bound))

    ! create mesh boundaries ...................................................

    do b = 1, size(generic_mesh%boundary)

      associate( mb => mesh % boundary(b)         &
               , gb => generic_mesh % boundary(b) )

        mb % id      = b
        mb % name    = gb % name
        mb % n_face  = size(gb % face)
        mb % coupled = gb % coupled

        allocate(mb % face(mb % n_face))

        do i = 1, mb % n_face

          e = gb % face(i) % element_id         ! corresponding element ID
          j = gb % face(i) % element_face       ! corresponding element face
          f = mesh % element(e) % face(j) % id  ! corresponding mesh face ID

          ! side of the mesh face on which the boundary is located
          if (mesh % face(f) % element(1) % id == e) then
            s = 2
          else
            s = 1
          end if

          mb % face(i) % mesh_face % id   = f
          mb % face(i) % mesh_face % side = s

          mb % face(i) % mesh_element % id   = e
          mb % face(i) % mesh_element % face = j

        end do

      end associate
    end do

  end subroutine ImportBoundaries

  !-----------------------------------------------------------------------------
  !> Import of element domains

  subroutine ImportElementDomains(mesh, generic_mesh)
    class(MeshPartition_3D), intent(inout) :: mesh         !< mesh partition
    class(GenericMesh_3D),    intent(in)    :: generic_mesh !< generic mesh

    type(StandardOperators_1D) :: sop
    type(EmbeddedInterpolation_1D), allocatable :: iop(:)
    real(RNP), allocatable :: VI(:,:)
    integer :: d, e, pg, pg_max, pg_min, po

    ! preliminaries ............................................................

    ! min/max order of element geometry representation
    pg_min = huge(pg_min)
    pg_max = 0
    do e = 1, size(generic_mesh % element)
      pg_min = min(pg_min, generic_mesh % element(e) % order)
      pg_max = max(pg_max, generic_mesh % element(e) % order)
      if (generic_mesh % element(e) % basis /= GAUSS_LOBATTO_BASIS) then
        call Error('ImportElementDomains', 'basis not supported')
      end if
    end do

    ! standard and interpolation operators
    mesh%p_geom = max(mesh%p_geom, pg_max)
    po  = mesh%p_geom
    sop = StandardOperators_1D(po)
    allocate(VI(0:po,0:po), iop(pg_min:pg_max))
    call sop % Get_Inverse_Legendre_VDM(VI)
    do pg = pg_min, pg_max
      if (pg == mesh%p_geom) cycle
      iop(pg) = EmbeddedInterpolation_1D(StandardOperators_1D(pg), sop%x)
    end do

    ! element points and approximate cuboids ...................................

    allocate(mesh % x_elem(0:po, 0:po, 0:po, mesh%n_elem, 3))
    allocate(mesh % x_cube(0:3, mesh%n_elem, 3))
    do e = 1, mesh % n_elem
      associate(xg => generic_mesh % element(e) % x)
        pg = generic_mesh % element(e) % order
        do d = 1, 3

          if (pg == po) then
            mesh % x_elem(:,:,:,e,d) = reshape(xg(:,d), [po+1,po+1,po+1])
          else
            call Interpolate(iop(pg), xg(:,d), mesh%x_elem(:,:,:,e,d))
          end if

          call LinearFit(VI, mesh%x_elem(:,:,:,e,d), mesh%x_cube(:,e,d))

        end do
      end associate
    end do

  contains

    !---------------------------------------------------------------------------
    !> Interpolates element points from generic mesh element

    subroutine Interpolate(iop, xo, xi)
      class(EmbeddedInterpolation_1D), intent(in)  :: iop
      real(RNP),                    intent(in)  :: xo(iop%no, iop%no, iop%no)
      real(RNP),                    intent(out) :: xi(iop%ni, iop%ni, iop%ni)

      real(RNP) :: z1(iop%ni, iop%no, iop%no)
      real(RNP) :: z2(iop%ni, iop%ni, iop%no)

      integer :: i, j, k, p

      ! interpolation in direction 1
      do k = 1, iop%no
      do j = 1, iop%no
      do i = 1, iop%ni
        z1(i,j,k) = 0
        do p = 1, iop%no
          z1(i,j,k) = z1(i,j,k) + iop%A(i,p) * xo(p,j,k)
        end do
      end do
      end do
      end do

      ! interpolation in direction 2
      do k = 1, iop%no
      do j = 1, iop%ni
      do i = 1, iop%ni
        z2(i,j,k) = 0
        do p = 1, iop%no
          z2(i,j,k) = z2(i,j,k) + iop%A(j,p) * z1(i,p,k)
        end do
      end do
      end do
      end do

      ! interpolation in direction 3
      do k = 1, iop%ni
      do j = 1, iop%ni
      do i = 1, iop%ni
        xi(i,j,k) = 0
        do p = 1, iop%no
          xi(i,j,k) = xi(i,j,k) + iop%A(k,p) * z1(i,j,p)
        end do
      end do
      end do
      end do

    end subroutine Interpolate

    !---------------------------------------------------------------------------
    !> Computes the constant and linear Legendre coefficients to given 3D
    !> Lobatto coefficients

    subroutine LinearFit(VI, x, y)
      real(RNP), intent(in)  :: VI(0:,0:)   !< Inverse Vandermonde matrix
      real(RNP), intent(in)  :: x(0:,0:,0:) !< Lobatto coefficients
      real(RNP), intent(out) :: y(0:3)      !< Legendre coefficients

      real(RNP) :: a(size(VI,1), size(VI,1)), b(size(VI,1)), c(0:1, 0:1, 0:1)
      integer   :: i, j, k, n

      n = size(VI,1)

      do i = 0, 1
        a = reshape(matmul(VI(i,:), reshape(x, [n,n*n])), [n,n])

        do j = 0, 1
          b = matmul(VI(j,:), a)
          do k = 0, 1
            c(i,j,k) = dot_product(VI(k,:), b)
          end do
        end do
      end do

      y = [ c(0,0,0), c(1,0,0), c(0,1,0), c(0,0,1) ]

    end subroutine LinearFit

  end subroutine ImportElementDomains

  !-----------------------------------------------------------------------------
  !> Identification of neighbor elements

  subroutine IdentifyNeighbors(mesh)
    class(MeshPartition_3D), intent(inout) :: mesh !< mesh partition

    integer, allocatable :: ne_face(:)     ! num elements per face
    integer, allocatable :: ne_edge(:)     ! num elements per edge
    integer, allocatable :: ne_vert(:)     ! num elements per vertex

    integer, allocatable :: nl_face(:)     ! num local elements per face
    integer, allocatable :: nl_edge(:)     ! num local elements per edge
    integer, allocatable :: nl_vert(:)     ! num local elements per vertex

    type NeighborElement
      integer :: id = 0  !< ID
      integer :: cc = 0  !< coupled neighbor component (face, edge or vertex)
    end type NeighborElement

    ! adjoining elements
    type(NeighborElement), allocatable :: elem_face(:,:) ! per face
    type(NeighborElement), allocatable :: elem_edge(:,:) ! per edge
    type(NeighborElement), allocatable :: elem_vert(:,:) ! per vertex

    integer :: f(6), e(12), v(8)
    integer :: i, j, k, l, m

    !---------------------------------------------------------------------------
    ! count number of adjoining elements for each vertex, edge and face

    ! local elements ...........................................................

    allocate(nl_face(mesh % n_face), source = 0)
    allocate(nl_edge(mesh % n_edge), source = 0)
    allocate(nl_vert(mesh % n_vert), source = 0)

    do l = 1, mesh % n_elem

      f = mesh % element(l) % face   % id
      e = mesh % element(l) % edge   % id
      v = mesh % element(l) % vertex % id

      nl_face(f) = nl_face(f) + 1
      nl_edge(e) = nl_edge(e) + 1
      nl_vert(v) = nl_vert(v) + 1

    end do

    !---------------------------------------------------------------------------
    ! establish lists of adjoining elements

    allocate(elem_face(maxval(nl_face), mesh%n_face))
    allocate(elem_edge(maxval(nl_edge), mesh%n_edge))
    allocate(elem_vert(maxval(nl_vert), mesh%n_vert))

    allocate(ne_face(mesh % n_face), source = 0)
    allocate(ne_edge(mesh % n_edge), source = 0)
    allocate(ne_vert(mesh % n_vert), source = 0)

    ! local neighbors ..........................................................

    do l = 1, mesh % n_elem

      f = mesh % element(l) % face   % id
      e = mesh % element(l) % edge   % id
      v = mesh % element(l) % vertex % id

      ne_face(f) = ne_face(f) + 1
      ne_edge(e) = ne_edge(e) + 1
      ne_vert(v) = ne_vert(v) + 1

      do j = 1, 6
        elem_face(ne_face(f(j)), f(j)) = NeighborElement(l, j)
      end do
      do j = 1, 12
        elem_edge(ne_edge(e(j)), e(j)) = NeighborElement(l, j)
      end do
      do j = 1, 8
        elem_vert(ne_vert(v(j)), v(j)) = NeighborElement(l, j)
      end do

    end do

    !---------------------------------------------------------------------------
    ! identify element neighbors

    do l = 1, mesh % n_elem
      associate(element => mesh % element(l))

        f = element % face   % id
        e = element % edge   % id
        v = element % vertex % id

        ! count and allocate neighbors
        i = -26  ! offset, accounting for self-references by own components
        do j = 1, 6
          i = i + ne_face(f(j))
        end do
        do j = 1, 12
          i = i + ne_edge(e(j))
        end do
        do j = 1, 8
          i = i + ne_vert(v(j))
        end do
        allocate(element % neighbor(i))

        i = 0

        ! faces
        do j = 1, 6
          element % face(j) % i_neighbor = i + 1
          m = f(j)
          ! local neighbors
          do k = 1, ne_face(m)

            if (elem_face(k,m) % id == l .and. k <= nl_face(m)) then
              cycle ! skip self-reference
            end if

            i = i + 1
            element % neighbor(i) % id   = elem_face(k,m) % id
            element % neighbor(i) % cc   = elem_face(k,m) % cc
            element % neighbor(i) % part = 0
          end do
        end do

        ! edges
        do j = 1, 12
          element % edge(j) % i_neighbor = i + 1
          m = e(j)
          ! local neighbors
          do k = 1, nl_edge(m)

            if (elem_edge(k,m) % id == l .and. k <= nl_edge(m)) then
              cycle ! skip self-reference
            end if

            i = i + 1
            element % neighbor(i) % id   = elem_edge(k,m) % id
            element % neighbor(i) % cc   = elem_edge(k,m) % cc
            element % neighbor(i) % part = 0
          end do
        end do

        ! vertices
        do j = 1, 8
          element % vertex(j) % i_neighbor = i + 1
          m = v(j)
          ! local neighbors
          do k = 1, nl_vert(m)

            if (elem_vert(k,m) % id == l .and. k <= nl_vert(m)) then
              cycle ! skip self-reference
            end if

            i = i + 1
            element % neighbor(i) % id   = elem_vert(k,m) % id
            element % neighbor(i) % cc   = elem_vert(k,m) % cc
            element % neighbor(i) % part = 0
          end do
        end do

      end associate
    end do

  end subroutine IdentifyNeighbors

  !=============================================================================

end submodule MP_ImportGenericMesh
