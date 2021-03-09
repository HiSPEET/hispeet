!> summary:  Import of generic 3d meshes
!> author:   Joerg Stiller
!> date:     2020/12/25
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(Mesh_3d__Partition) MP_ImportGenericMesh
  use Execution_Control
  use Gauss_Jacobi
  use Standard_Operators_1D
  use Embedded_Interpolation
  use Generic_Mesh_3d
  use Mesh_3d__Element_Indexing, only: V_FACE, E_FACE
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Import a generic 3d mesh

  module subroutine ImportGenericMesh(mesh, generic_mesh, comm)
    class(Mesh3d_Partition), intent(out) :: mesh         !< mesh partition
    class(GenericMesh3d),    intent(in)  :: generic_mesh !< generic mesh
    type(MPI_Comm),          intent(in)  :: comm         !< MPI communicator

    integer :: rank

    call MPI_Comm_rank(comm, rank)
    if (rank > 0) return

    ! check prerequisites ......................................................

    if (generic_mesh % numbering /= LEXICAL_NUMBERING) then
      call Error('ImportGenericMesh', 'Lexical numbering is required')
    else if (any(generic_mesh % element%typ /= HEXAHEDRAL_ELEMENT)) then
      call Error('ImportGenericMesh', 'Only hexahedral elements supported')
    end if

    mesh % n_vert      =  size(generic_mesh % vertex)
    mesh % n_elem      =  size(generic_mesh % element)
    mesh % n_part      =  1
    mesh % part        =  0
    mesh % structured  = .false.
    mesh % regular     = .false.
    mesh % comm        =  comm

    call ImportElements(mesh, generic_mesh)

    call mesh % IdentifyEdges()
    call mesh % BuildFaces()

    call ImportElementDomains (mesh, generic_mesh)
    call ImportBoundaries     (mesh, generic_mesh)

    call mesh % BuildLinks()
    call mesh % IdentifyRanks()
    call mesh % BuildGhosts()

  end subroutine ImportGenericMesh

  !-----------------------------------------------------------------------------
  !> Import of mesh elements
  !>
  !> Allocates the elements, sets global and local IDs, and adopts the vertex
  !> IDs from then generic mesh.

  subroutine ImportElements(mesh, generic_mesh)
    class(Mesh3d_Partition), intent(inout) :: mesh         !< mesh partition
    class(GenericMesh3d),    intent(in)    :: generic_mesh !< generic mesh

    integer :: i

    allocate(mesh % element( mesh % n_elem ))

    do i = 1, mesh % n_elem
      mesh % element(i) % global_id = i
      mesh % element(i) % local_id  = i
      mesh % element(i) % vertex % id = generic_mesh % element(i) % vertex
    end do

  end subroutine ImportElements

  !-----------------------------------------------------------------------------
  !> Import of mesh boundaries

  subroutine ImportBoundaries(mesh, generic_mesh)
    class(Mesh3d_Partition), intent(inout) :: mesh         !< mesh partition
    class(GenericMesh3d),    intent(in)    :: generic_mesh !< generic mesh

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
          if (mesh % face(f) % element(1) %id == e) then
            s = 2
          else
            s = 1
          end if

          mb % face(i) % mesh_face % id   = f
          mb % face(i) % mesh_face % side = s

          mb % face(i) % mesh_element % id   = e
          mb % face(i) % mesh_element % face = j

          mesh % element(e) % face(j) % boundary = i

        end do

      end associate
    end do

  end subroutine ImportBoundaries

  !-----------------------------------------------------------------------------
  !> Identification of neighbor elements

  subroutine IdentifyNeighbors(mesh)
    class(Mesh3d_Partition), intent(inout) :: mesh !< mesh partition

    integer, allocatable :: ne_face(:)     ! num elements per face
    integer, allocatable :: ne_edge(:)     ! num elements per edge
    integer, allocatable :: ne_vert(:)     ! num elements per vertex

    integer, allocatable :: nl_face(:)     ! num local elements per face
    integer, allocatable :: nl_edge(:)     ! num local elements per edge
    integer, allocatable :: nl_vert(:)     ! num local elements per vertex

    integer, allocatable :: nc_face(:)     ! num coupled elements per face
    integer, allocatable :: nc_edge(:)     ! num coupled elements per edge
    integer, allocatable :: nc_vert(:)     ! num coupled elements per vertex

    logical, allocatable :: marked_edge(:) ! marker for edges
    logical, allocatable :: marked_vert(:) ! marker for vertices

    type NeighborElement
      integer :: id = 0  !< ID
      integer :: cc = 0  !< coupled neighbor component (face, edge or vertex)
    end type NeighborElement

    ! adjoining elements
    type(NeighborElement), allocatable :: elem_face(:,:) ! per face
    type(NeighborElement), allocatable :: elem_edge(:,:) ! per edge
    type(NeighborElement), allocatable :: elem_vert(:,:) ! per vertex

    integer :: f(6), e(12), v(8)
    integer :: eb(4), fb, kb, lb, mb, nb, vb(4)
    integer :: ec(4), fc, kc, lc, mc, nc, vc(4)
    integer :: b, c, i, j, k, l, m

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

    ! coupled elements .........................................................

    allocate(nc_face(mesh % n_face), source = 0)
    allocate(nc_edge(mesh % n_edge), source = 0)
    allocate(nc_vert(mesh % n_vert), source = 0)

    allocate(marked_edge(mesh % n_edge))
    allocate(marked_vert(mesh % n_vert))

    do b = 1, mesh % n_bound
      c = mesh % boundary(b) % coupled
      if (c > b) then

        marked_edge = .false.
        marked_vert = .false.

        do i = 1, mesh % boundary(b) % n_face

          ! corresponding element ID and face on boundary b
          lb = mesh % boundary(b) % face(i) % mesh_element % id
          mb = mesh % boundary(b) % face(i) % mesh_element % face

          ! corresponding mesh face, edges and vertices on boundary b
          fb = mesh % element(lb) % face   (            mb ) % id
          eb = mesh % element(lb) % edge   (E_FACE(1:4, mb)) % id
          vb = mesh % element(lb) % vertex (V_FACE(1:4, mb)) % id

          ! corresponding element ID and face on boundary c
          lc = mesh % boundary(b) % face(i) % mesh_element % id
          mc = mesh % boundary(b) % face(i) % mesh_element % face

          ! corresponding mesh face, edges and vertices on boundary c
          fc = mesh % element(lc) % face   (            mc ) % id
          ec = mesh % element(lc) % edge   (E_FACE(1:4, mc)) % id
          vc = mesh % element(lc) % vertex (V_FACE(1:4, mc)) % id

          ! increase coupled face counts
          nc_face(fb) = nc_face(fb) + nl_face(fc)
          nc_face(fc) = nc_face(fc) + nl_face(fb)

          ! increase coupled edge counts
          where(.not. marked_edge(ec))
            nc_edge(eb) = nc_edge(eb) + nl_edge(ec)
            nc_edge(ec) = nc_edge(ec) + nl_edge(eb)
            marked_edge(ec) = .true.
          end where

          ! increase coupled vertex counts
          where(.not. marked_vert(vc))
            nc_vert(vb) = nc_vert(vb) + nl_vert(vc)
            nc_vert(vc) = nc_vert(vc) + nl_vert(vb)
            marked_vert(vc) = .true.
          end where

        end do

      end if
    end do

    !---------------------------------------------------------------------------
    ! establish lists of adjoining elements

    allocate(elem_face(maxval(nl_face) + maxval(nc_face), mesh%n_face))
    allocate(elem_edge(maxval(nl_edge) + maxval(nc_edge), mesh%n_edge))
    allocate(elem_vert(maxval(nl_vert) + maxval(nc_vert), mesh%n_vert))

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

    ! neighbors coupled via periodic boundaries ................................

    do b = 1, mesh % n_bound
      c = mesh % boundary(b) % coupled
      if (c > b) then

        marked_edge = .false.
        marked_vert = .false.

        do i = 1, mesh % boundary(b) % n_face

          ! corresponding element ID and face on boundary b
          lb = mesh % boundary(b) % face(i) % mesh_element % id
          mb = mesh % boundary(b) % face(i) % mesh_element % face

          ! corresponding mesh face, edges and vertices on boundary b
          fb = mesh % element(lb) % face   (            mb ) % id
          eb = mesh % element(lb) % edge   (E_FACE(1:4, mb)) % id
          vb = mesh % element(lb) % vertex (V_FACE(1:4, mb)) % id

          ! corresponding element ID and face on boundary c
          lc = mesh % boundary(b) % face(i) % mesh_element % id
          mc = mesh % boundary(b) % face(i) % mesh_element % face

          ! corresponding mesh face, edges and vertices on boundary c
          fc = mesh % element(lc) % face   (            mc ) % id
          ec = mesh % element(lc) % edge   (E_FACE(1:4, mc)) % id
          vc = mesh % element(lc) % vertex (V_FACE(1:4, mc)) % id

          ! append coupled face neighbors
          nb = nl_face(fb) ! number of local elements adjoining to face fb
          nc = nl_face(fc) ! number of local elements adjoining to face fc
          ! append adjoining elements of coupled face
          elem_face(ne_face(fb)+1:ne_face(fb) + nc, fb) = elem_face(1:nc, fc)
          elem_face(ne_face(fc)+1:ne_face(fc) + nb, fc) = elem_face(1:nb, fb)
          ! increase counters
          ne_face(fb) = ne_face(fb) + nc
          ne_face(fc) = ne_face(fc) + nb

          ! append coupled edge neighbors
          do j = 1, 4
            if (marked_edge(ec(j))) cycle
            kb = eb(j)
            kc = ec(j)
            nb = nl_edge(kb) ! num local elements adjoining to edge kb
            nc = nl_edge(kc) ! num local elements adjoining to edge kc
            ! append adjoining elements of coupled edge
            elem_edge(ne_edge(kb)+1:ne_edge(kb) + nc, kb) = elem_edge(1:nc, kc)
            elem_edge(ne_edge(kc)+1:ne_edge(kc) + nb, kc) = elem_edge(1:nb, kb)
            ! increase counters
            ne_edge(kb) = ne_edge(kb) + nc
            ne_edge(kc) = ne_edge(kc) + nb
            ! mark coupled edges as processed
            marked_edge(kc) = .true.
          end do

          ! append coupled vertex neighbors
          do j = 1, 4
            if (marked_vert(vc(j))) cycle
            kb = vb(j)
            kc = vc(j)
            nb = nl_vert(kb) ! num local elements adjoining to vertex kb
            nc = nl_vert(kc) ! num local elements adjoining to vertex kc
            ! append adjoining elements of coupled vertices
            elem_vert(ne_vert(kb)+1:ne_vert(kb) + nc, kb) = elem_vert(1:nc, kc)
            elem_vert(ne_vert(kc)+1:ne_vert(kc) + nb, kc) = elem_vert(1:nb, kb)
            ! increase counters
            ne_vert(kb) = ne_vert(kb) + nc
            ne_vert(kc) = ne_vert(kc) + nb
            ! mark coupled vertices as processed
            marked_vert(kc) = .true.
          end do

        end do

      end if
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

  !-----------------------------------------------------------------------------
  !> Import of element domains

  subroutine ImportElementDomains(mesh, generic_mesh)
    class(Mesh3d_Partition), intent(inout) :: mesh         !< mesh partition
    class(GenericMesh3d),    intent(in)    :: generic_mesh !< generic mesh

    type(InterpolationOperator), allocatable :: iop(:)
    real(RNP), allocatable :: xi(:)
    integer :: c, e, pg, pg_max, pg_min, po

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

    ! interpolation operators
    mesh%p_geom = max(mesh%p_geom, pg_max)
    po = mesh%p_geom
    allocate(xi(0:po), source = GLL_Points(po))
    allocate(iop(pg_min:pg_max))
    do pg = pg_min, pg_max
      if (pg == mesh%p_geom) cycle
      iop(pg) = InterpolationOperator(StandardOperators1D(pg), xi)
    end do

    ! import/interpolate element points ........................................

    allocate(mesh % x_elem(0:po, 0:po, 0:po, mesh%n_elem, 3))

    do e = 1, mesh % n_elem
      pg = generic_mesh % element(e) % order
      associate(xg => generic_mesh%element(e)%x, xe => mesh % x_elem)
        if (pg == po) then
          do c = 1, 3
            xe(:,:,:,e,c) = reshape(xg(:,c), [po+1,po+1,po+1])
          end do
        else
          do c = 1, 3
            call Interpolate(iop(pg), xg(:,c), xe(:,:,:,e,c))
          end do
        end if
      end associate
    end do

  contains

    subroutine Interpolate(iop, xo, xi)
      class(InterpolationOperator), intent(in)  :: iop
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

  end subroutine ImportElementDomains

  !=============================================================================

end submodule MP_ImportGenericMesh
