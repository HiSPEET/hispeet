module Mesh_3d__Generate_Regular_Mesh

  use Kind_Parameters  , only: IXL, IXS, RNP
  use Execution_Control, only: Error
  use Gauss_Jacobi
  use XMPI

  use Mesh_3d__Boundary
  use Mesh_3d__Element
  use Mesh_3d__Element_Indexing
  use Mesh_3d__Partition
  use Mesh_3d__Structured_Indexing
  implicit none
  private

  public :: GenerateRegularMesh

contains

  subroutine GenerateRegularMesh(mesh, np, ep, xo, dx, periodic, comm, pg)

    type(Mesh3d_Partition), intent(out) :: mesh !< local partition
    integer,   intent(in) :: np(3)       !< num partitions in directions 1:3
    integer,   intent(in) :: ep(3)       !< elements per partition and direction
    real(RNP), intent(in) :: xo(3)       !< corner closest to -infinity
    real(RNP), intent(in) :: dx(3)       !< mesh spacing per direction
    logical,   intent(in) :: periodic(3) !< set true for periodic directions
    type(MPI_Comm), intent(in) :: comm   !< MPI communicator
    integer, optional, intent(in) :: pg  !< polynomial order of geometry [1]

    ! local variables ..........................................................

    integer :: part       ! partition ID
    integer :: n_part     ! number of non-empty partitions
    integer :: n_bound    ! number of external boundaries
    integer :: p_geom     ! polynomial order of element geometry
    integer :: n1, n2, n3 ! local mesh dimensions

    ! MPI
    integer :: comm_size
    integer :: rank

    ! initialization ...........................................................

    call MPI_Comm_size(comm, comm_size)
    call MPI_Comm_rank(comm, rank)

    ! number of partitions
    n_part = product(np)
    if (n_part > comm_size .and. rank == 0) then
      call Error('GenerateRegularMesh', 'n_part > MPI communicator size')
    end if

    if (present(pg)) then
      p_geom = pg
    else
      p_geom = 1
    end if

    n_bound = 6

    ! partition ID
    if (rank < n_part) then
      part = rank
    else
      ! create empty partition
      mesh = Mesh3d_Partition( p_geom, n_bound, n_part &
                             , comm       =  comm      &
                             , structured = .true.     &
                             , regular    = .true.     &
                             , dx         =  dx        )
      return
    end if

    ! local mesh dimensions
    n1 = ep(1)
    n2 = ep(2)
    n3 = ep(3)

    ! mesh attributes and properties ...........................................

    mesh % n_vert     =  NumberOfVertices (n1, n2, n3)
    mesh % n_edge     =  NumberOfEdges    (n1, n2, n3)
    mesh % n_face     =  NumberOfFaces    (n1, n2, n3)
    mesh % n_elem     =  NumberOfElements (n1, n2, n3)
    mesh % p_geom     =  p_geom
    mesh % n_bound    =  n_bound
    mesh % n_part     =  n_part
    mesh % part       =  part

    mesh % structured = .true.
    mesh % n_elem_1   =  n1
    mesh % n_elem_2   =  n2
    mesh % n_elem_3   =  n3

    mesh % regular    = .true.
    mesh % n_face_1   =  (n1+1) * n2 * n3
    mesh % n_face_2   =  (n2+1) * n3 * n1
    mesh % n_face_3   =  (n3+1) * n1 * n2
    mesh % dx         =  dx

    mesh % comm       =  comm

    call GenerateRegularElements(mesh, np, ep, periodic)
    call GenerateRegularElementDomains(mesh, xo)
    call GenerateRegularFaces(mesh)
    call GenerateRegularMeshBoundaries(mesh, periodic)

    call mesh % BuildLinks()
    call mesh % IdentifyPrimaries()
    call mesh % BuildGhosts()

  end subroutine GenerateRegularMesh

  !-----------------------------------------------------------------------------
  !> Builds the local mesh elements

  subroutine GenerateRegularElements(mesh, np, ep, periodic)

    class(Mesh3d_Partition), intent(inout) :: mesh  !< local partition
    integer, intent(in) :: np(3) !< num partitions in directions 1:3
    integer, intent(in) :: ep(3) !< num elements per partition and direction
    logical, intent(in) :: periodic(3) !< set true for periodic directions

    integer(IXL) :: ig, jg, kg, ng(3)  ! global elements indices and counts
    integer      :: ie, je, ke         ! local element indices
    integer      :: ip, jp, kp         ! partion triple index
    integer      :: i0, j0, k0         ! offsets WRT global numbering
    integer      :: ies, jes, kes      ! shifted element indices
    integer      :: ips, jps, kps      ! shifted partition indices
    integer      :: i, e, l, r, s, t

    type(Mesh3d_ElementNeighbor) :: neighbor(26)

    ! prerequisites ............................................................

    ! global mesh dimensions
    ng = np * ep

    call TripleIndex(ip, jp, kp, np(1), np(2), 1, 1, 1, l=mesh%part+1)

    ! offset of local element indices WRT global numbering
    i0 = ep(1) * (ip - 1)
    j0 = ep(2) * (jp - 1)
    k0 = ep(3) * (kp - 1)

    ! generate elements ........................................................

    allocate(mesh%element( mesh % n_elem ))

    associate( element => mesh % element   &
             , n1      => mesh % n_elem_1  &
             , n2      => mesh % n_elem_2  &
             , n3      => mesh % n_elem_3  )

      do ke = 1, n3
      do je = 1, n2
      do ie = 1, n1

        ! local and global element IDs .........................................

        e = LexicalElementIndex(ie, je, ke, n1, n2)

        ig = i0 + ie
        jg = j0 + je
        kg = k0 + ke

        element(e) % global_id = LexicalElementIndex(ig,jg,kg,ng(1),ng(2))
        element(e) % local_id  = e

        ! element vertices .....................................................

        element(e) % vertex(1) % id = LexicalVertexIndex(ie-1,je-1,ke-1, n1,n2)
        element(e) % vertex(2) % id = LexicalVertexIndex(ie  ,je-1,ke-1, n1,n2)
        element(e) % vertex(3) % id = LexicalVertexIndex(ie-1,je  ,ke-1, n1,n2)
        element(e) % vertex(4) % id = LexicalVertexIndex(ie  ,je  ,ke-1, n1,n2)
        element(e) % vertex(5) % id = LexicalVertexIndex(ie-1,je-1,ke  , n1,n2)
        element(e) % vertex(6) % id = LexicalVertexIndex(ie  ,je-1,ke  , n1,n2)
        element(e) % vertex(7) % id = LexicalVertexIndex(ie-1,je  ,ke  , n1,n2)
        element(e) % vertex(8) % id = LexicalVertexIndex(ie  ,je  ,ke  , n1,n2)

        ! element edges ........................................................

        element(e) % edge( 1) % id = LexicalEdgeIndex(ie,je-1,ke-1, 1, n1,n2,n3)
        element(e) % edge( 2) % id = LexicalEdgeIndex(ie,je  ,ke-1, 1, n1,n2,n3)
        element(e) % edge( 3) % id = LexicalEdgeIndex(ie,je-1,ke  , 1, n1,n2,n3)
        element(e) % edge( 4) % id = LexicalEdgeIndex(ie,je  ,ke  , 1, n1,n2,n3)

        element(e) % edge( 5) % id = LexicalEdgeIndex(ie-1,je,ke-1, 2, n1,n2,n3)
        element(e) % edge( 6) % id = LexicalEdgeIndex(ie  ,je,ke-1, 2, n1,n2,n3)
        element(e) % edge( 7) % id = LexicalEdgeIndex(ie-1,je,ke  , 2, n1,n2,n3)
        element(e) % edge( 8) % id = LexicalEdgeIndex(ie  ,je,ke  , 2, n1,n2,n3)

        element(e) % edge( 9) % id = LexicalEdgeIndex(ie-1,je-1,ke, 3, n1,n2,n3)
        element(e) % edge(10) % id = LexicalEdgeIndex(ie  ,je-1,ke, 3, n1,n2,n3)
        element(e) % edge(11) % id = LexicalEdgeIndex(ie-1,je  ,ke, 3, n1,n2,n3)
        element(e) % edge(12) % id = LexicalEdgeIndex(ie  ,je  ,ke, 3, n1,n2,n3)

        ! element faces ........................................................

        element(e) % face(1) % id = LexicalFaceIndex(ie-1,je,ke, 1, n1,n2,n3)
        element(e) % face(2) % id = LexicalFaceIndex(ie  ,je,ke, 1, n1,n2,n3)

        element(e) % face(3) % id = LexicalFaceIndex(ie,je-1,ke, 2, n1,n2,n3)
        element(e) % face(4) % id = LexicalFaceIndex(ie,je  ,ke, 2, n1,n2,n3)

        element(e) % face(5) % id = LexicalFaceIndex(ie,je,ke-1, 3, n1,n2,n3)
        element(e) % face(6) % id = LexicalFaceIndex(ie,je,ke  , 3, n1,n2,n3)

        if (ie == 1  .and. ip == 1    )  element(e) % face(1) % boundary = 1
        if (ie == n1 .and. ip == np(1))  element(e) % face(2) % boundary = 2
        if (je == 1  .and. jp == 1    )  element(e) % face(3) % boundary = 3
        if (je == n2 .and. jp == np(2))  element(e) % face(4) % boundary = 4
        if (ke == 1  .and. kp == 1    )  element(e) % face(5) % boundary = 5
        if (ke == n3 .and. kp == np(3))  element(e) % face(6) % boundary = 6

        ! neighbors and boundaries .............................................

        neighbor = Mesh3d_ElementNeighbor()

        do t = -1, 1  ! ζ
        do s = -1, 1  ! η
        do r = -1, 1  ! ξ

          ! skip element center
          if (r == 0 .and. s == 0 .and. t == 0) cycle

          ! neighbor element triple index within local partition
          ies = ie + r
          jes = je + s
          kes = ke + t

          ! first index of the neighbor partition
          ips = ip
          if (ies < 1) then
            ies = n1
            ips = ip - 1
            if (ips < 1) then
              if (periodic(1)) then
                ips = np(1)
              else
                if (ElementComponentType(r,s,t) == IS_FACE) then
                  element(e) % face(ElementFaceID(r,s,t)) % boundary = 1
                end if
                ips = -1
              end if
            end if
          else if (ies > n1) then
            ies = 1
            ips = ip + 1
            if (ips > np(1)) then
              if (periodic(1)) then
                ips = 1
              else
                if (ElementComponentType(r,s,t) == IS_FACE) then
                  element(e) % face(ElementFaceID(r,s,t)) % boundary = 2
                end if
                ips = -1
              end if
            end if
          end if

          ! second index of the neighbor partition
          jps = jp
          if (jes < 1) then
            jes = n2
            jps = jp - 1
            if (jps < 1) then
              if (periodic(2)) then
                jps = np(2)
              else
                if (ElementComponentType(r,s,t) == IS_FACE) then
                  element(e) % face(ElementFaceID(r,s,t)) % boundary = 3
                end if
                jps = -1
              end if
            end if
          else if (jes > n2) then
            jes = 1
            jps = jp + 1
            if (jps > np(2)) then
              if (periodic(2)) then
                jps = 1
              else
                if (ElementComponentType(r,s,t) == IS_FACE) then
                  element(e) % face(ElementFaceID(r,s,t)) % boundary = 4
                end if
                jps = -1
              end if
            end if
          end if

          ! third index of the neighbor partition
          kps = kp
          if (kes < 1) then
            kes = n3
            kps = kp - 1
            if (kps < 1) then
              if (periodic(3)) then
                kps = np(3)
              else
                if (ElementComponentType(r,s,t) == IS_FACE) then
                  element(e) % face(ElementFaceID(r,s,t)) % boundary = 5
                end if
                kps = -1
              end if
            end if
          else if (kes > n3) then
            kes = 1
            kps = kp + 1
            if (kps > np(3)) then
              if (periodic(3)) then
                kps = 1
              else
                if (ElementComponentType(r,s,t) == IS_FACE) then
                  element(e) % face(ElementFaceID(r,s,t)) % boundary = 6
                end if
                kps = -1
              end if
            end if
          end if

          ! skip nonexisting neighbor
          if (ips < 0 .or. jps < 0 .or. kps < 0) cycle

          ! linear index of present element component
          i = ElementComponentIndex(r,s,t)

          ! neighbor element ID in its home partition
          neighbor(i) % id   = LexicalElementIndex(ies, jes, kes, n1, n2)

          ! neighbor element home partition ID
          neighbor(i) % part = LexicalIndex(ips,jps,kps,np(1),np(2),1,1,1) - 1

          ! neighbor element component
          select case(ElementComponentType(r,s,t))
          case(IS_FACE)
            neighbor(i) % cc = ElementFaceID(-r,-s,-t)
          case(IS_EDGE)
            neighbor(i) % cc = ElementEdgeID(-r,-s,-t)
          case(IS_VERTEX)
            neighbor(i) % cc = ElementVertexID(-r,-s,-t)
          end select

        end do
        end do
        end do

        allocate( element(e) % neighbor( count(neighbor % id > 0) ))

        ! store connectivity into element
        l = 0
        do i = 1, 26
          if (neighbor(i) % id > 0) then
            l = l + 1
            element(e) % neighbor(l) = neighbor(i)
            select case(ElementComponentType(i))
            case(IS_FACE)
              s = ElementFaceID(i)
              element(e) % face(s) % n_neighbor = 1
              element(e) % face(s) % i_neighbor = int(l, IXS)
            case(IS_EDGE)
              s = ElementEdgeID(i)
              element(e) % edge(s) % n_neighbor = 1
              element(e) % edge(s) % i_neighbor = int(l, IXS)
            case(IS_VERTEX)
              s = ElementVertexID(i)
              element(e) % vertex(s) % n_neighbor = 1
              element(e) % vertex(s) % i_neighbor = int(l, IXS)
            end select
          end if
        end do

      end do
      end do
      end do

    end associate

  end subroutine GenerateRegularElements

  !-----------------------------------------------------------------------------
  !> Builds the local mesh faces

  subroutine GenerateRegularFaces(mesh)

    class(Mesh3d_Partition), intent(inout) :: mesh  !< local partition

    integer :: i, j, k, l

    allocate(mesh % face( mesh % n_face ))

    associate( face => mesh % face      &
             , n1   => mesh % n_elem_1  &
             , n2   => mesh % n_elem_2  &
             , n3   => mesh % n_elem_3  )

      ! x1-faces ...............................................................

      do k = 1, n3
      do j = 1, n2
      do i = 0, n1

        l = LexicalFaceIndex(i, j, k, 1, n1, n2, n3)

        if (i > 0) then
          face(l) % element(1) % id   = LexicalElementIndex(i  , j, k, n1, n2)
          face(l) % element(1) % face = 2
        end if
        if (i < n1) then
          face(l) % element(2) % id   = LexicalElementIndex(i+1, j, k, n1, n2)
          face(l) % element(2) % face = 1
        end if

      end do
      end do
      end do

      ! x2-faces ...............................................................

      do k = 1, n3
      do j = 0, n2
      do i = 1, n1

        l = LexicalFaceIndex(i, j, k, 2, n1, n2, n3)

        if (j > 0) then
          face(l) % element(1) % id   = LexicalElementIndex(i, j  , k, n1, n2)
          face(l) % element(1) % face = 4
        end if
        if (j < n2) then
          face(l) % element(2) % id   = LexicalElementIndex(i, j+1, k, n1, n2)
          face(l) % element(2) % face = 3
        end if

      end do
      end do
      end do

      ! x3-faces ...............................................................

      do k = 0, n3
      do j = 0, n2
      do i = 0, n1

        l = LexicalFaceIndex(i, j, k, 3, n1, n2, n3)

        if (k > 0) then
          face(l) % element(1) % id   = LexicalElementIndex(i, j, k  , n1, n2)
          face(l) % element(1) % face = 6
        end if

        if (k < n3) then
          face(l) % element(2) % id   = LexicalElementIndex(i, j, k+1, n1, n2)
          face(l) % element(2) % face = 5
        end if

      end do
      end do
      end do

    end associate

  end subroutine GenerateRegularFaces

  !-----------------------------------------------------------------------------
  !> Creates the boundaries of a structured mesh partitions

  subroutine GenerateRegularMeshBoundaries(mesh, periodic)

    class(Mesh3d_Partition), intent(inout) :: mesh !< local partition
    logical, intent(in) :: periodic(3)             !< T for periodic directions

    integer :: b, e, f, i, j, k, l
    integer :: coupled(6)

    allocate(mesh % boundary(6))

    coupled = 0
    if (periodic(1)) coupled([1,2]) = [2,1]
    if (periodic(2)) coupled([3,4]) = [4,3]
    if (periodic(3)) coupled([5,6]) = [6,5]

    associate( boundary => mesh % boundary  &
             , element  => mesh%element     &
             , n1       => mesh % n_elem_1  &
             , n2       => mesh % n_elem_2  &
             , n3       => mesh % n_elem_3  )

      ! west ...................................................................

      b = 1
      i = 1
      boundary(b) = Mesh3d_Boundary( b, 'west', n2*n3, coupled(b) )
      f = 1
      do k = 1, n3
      do j = 1, n2
        e = LexicalElementIndex(i, j, k, n1, n2)
        l = LexicalFaceIndex(i-1, j, k, 1, n1, n2, n3)
        boundary(b) % face(f) % mesh_face    % id   = l
        boundary(b) % face(f) % mesh_face    % side = 2
        boundary(b) % face(f) % mesh_element % id   = e
        boundary(b) % face(f) % mesh_element % face = b
        f = f + 1
      end do
      end do

      ! east ...................................................................

      b = 2
      i = n1
      boundary(b) = Mesh3d_Boundary( b, 'east', n2*n3, coupled(b) )
      f = 1
      do k = 1, n3
      do j = 1, n2
        e = LexicalElementIndex(i, j, k, n1, n2)
        l = LexicalFaceIndex(i, j, k, 1, n1, n2, n3)
        boundary(b) % face(f) % mesh_face    % id   = l
        boundary(b) % face(f) % mesh_face    % side = 1
        boundary(b) % face(f) % mesh_element % id   = e
        boundary(b) % face(f) % mesh_element % face = b
        f = f + 1
      end do
      end do

      ! south ..................................................................

      b = 3
      j = 1
      boundary(b) = Mesh3d_Boundary( b, 'south', n1*n3, coupled(b) )
      f = 1
      do k = 1, n3
      do i = 1, n1
        e = LexicalElementIndex(i, j, k, n1, n2)
        l = LexicalFaceIndex(i, j-1, k, 2, n1, n2, n3)
        boundary(b) % face(f) % mesh_face    % id   = l
        boundary(b) % face(f) % mesh_face    % side = 2
        boundary(b) % face(f) % mesh_element % id   = e
        boundary(b) % face(f) % mesh_element % face = b
        f = f + 1
      end do
      end do

      ! north ..................................................................

      b = 4
      j = n2
      boundary(b) = Mesh3d_Boundary( b, 'north', n1*n3, coupled(b) )
      f = 1
      do k = 1, n3
      do i = 1, n1
        e = LexicalElementIndex(i, j, k, n1, n2)
        l = LexicalFaceIndex(i, j, k, 2, n1, n2, n3)
        boundary(b) % face(f) % mesh_face    % id   = l
        boundary(b) % face(f) % mesh_face    % side = 1
        boundary(b) % face(f) % mesh_element % id   = e
        boundary(b) % face(f) % mesh_element % face = b
        f = f + 1
      end do
      end do

      ! bottom .................................................................

      b = 5
      k = 1
      boundary(b) = Mesh3d_Boundary( b, 'bottom', n1*n2, coupled(b) )
      f = 1
      do j = 1, n2
      do i = 1, n1
        e = LexicalElementIndex(i, j, k, n1, n2)
        l = LexicalFaceIndex(i, j, k-1, 3, n1, n2, n3)
        boundary(b) % face(f) % mesh_face    % id   = l
        boundary(b) % face(f) % mesh_face    % side = 2
        boundary(b) % face(f) % mesh_element % id   = e
        boundary(b) % face(f) % mesh_element % face = b
        f = f + 1
      end do
      end do

      ! top ....................................................................

      b = 6
      k = n3
      boundary(b) = Mesh3d_Boundary( b, 'top', n1*n2, coupled(b) )
      f = 1
      do j = 1, n2
      do i = 1, n1
        e = LexicalElementIndex(i, j, k, n1, n2)
        l = LexicalFaceIndex(i, j, k, 3, n1, n2, n3)
        boundary(b) % face(f) % mesh_face    % id   = l
        boundary(b) % face(f) % mesh_face    % side = 1
        boundary(b) % face(f) % mesh_element % id   = e
        boundary(b) % face(f) % mesh_element % face = b
        f = f + 1
      end do
      end do

    end associate

  end subroutine GenerateRegularMeshBoundaries

  !-----------------------------------------------------------------------------
  !> Creates the element domains

  subroutine GenerateRegularElementDomains(mesh, xo)
    class(Mesh3d_Partition), intent(inout) :: mesh  !< local partition
    real(RNP)              , intent(in)    :: xo(3) !< corner closest to -∞

    real(RNP), dimension(0 : mesh%p_geom)   :: x1, x2, x3, ys
    integer :: e, i, j, k, l, m, n

    ! GLL points transformed to [-1,0]
    ys = (GLL_Points(mesh % p_geom) - 1) / 2

    associate( n1 => mesh % n_elem_1  &
             , n2 => mesh % n_elem_2  &
             , n3 => mesh % n_elem_3  &
             , pg => mesh % p_geom    &
             , dx => mesh % dx        )

      allocate(mesh % x_elem(0:pg, 0:pg, 0:pg, 1:mesh%n_elem, 1:3))

      do n = 1, n3
      do m = 1, n2
      do l = 1, n1

        e = LexicalElementIndex(l, m, n, n1, n2)

        ! 1D point distributions
        x1 = xo(1) + (l + ys) * dx(1)
        x2 = xo(2) + (m + ys) * dx(2)
        x3 = xo(3) + (n + ys) * dx(3)

        ! element points
        do k = 0, pg
        do j = 0, pg
        do i = 0, pg
          mesh % x_elem(i,j,k,e,1) = x1(i)
          mesh % x_elem(i,j,k,e,2) = x2(j)
          mesh % x_elem(i,j,k,e,3) = x3(k)
        end do
        end do
        end do

      end do
      end do
      end do

    end associate

  end subroutine GenerateRegularElementDomains

  !=============================================================================

end module Mesh_3d__Generate_Regular_Mesh
