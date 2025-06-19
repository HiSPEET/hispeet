module Generate_Regular_Mesh__3D

  use Kind_Parameters  , only: IXS, RNP
  use Constants        , only: ZERO, HALF
  use Execution_Control, only: Error
  use Gauss_Jacobi
  use XMPI

  use Affine_Transformation__3D
  use Mesh_Boundary__3D
  use Mesh_Element__3D
  use Mesh_Element_Indexing__3D
  use Mesh__3D
  use Structured_Mesh_Indexing__3D
  implicit none
  private

  public :: GenerateRegularMesh

contains

  subroutine GenerateRegularMesh(mesh, np, ep, xo, dx, periodic, comm, pg)

    type(Mesh_3D), intent(out) :: mesh   !< local partition
    integer,   intent(in) :: np(3)       !< num partitions in directions 1:3
    integer,   intent(in) :: ep(3)       !< elements per partition and direction
    real(RNP), intent(in) :: xo(3)       !< corner closest to -infinity
    real(RNP), intent(in) :: dx(3)       !< mesh spacing per direction
    logical,   intent(in) :: periodic(3) !< set true for periodic directions
    type(MPI_Comm), intent(in) :: comm   !< MPI "world" communicator
    integer, optional, intent(in) :: pg  !< polynomial order of geometry [1]

    ! local variables ..........................................................

    integer :: n_parts       ! number of non-empty partitions
    integer :: n1, n2, n3    ! local mesh dimensions
    integer :: i0, j0, k0    ! element offsets WRT global numbering
    logical :: self(3)       ! indicator wether linked to itself

    type(MeshAttributes_3D) :: attrib

    ! MPI
    integer :: comm_size
    integer :: rank

    !---------------------------------------------------------------------------
    ! initialization

    call MPI_Comm_size(comm, comm_size)
    call MPI_Comm_rank(comm, rank)

    ! number of partitions
    n_parts = product(np)
    if (n_parts > comm_size .and. rank == 0) then
      call Error('GenerateRegularMesh', 'n_parts > MPI communicator size')
    end if

    !---------------------------------------------------------------------------
    ! active partitions: essentials

    if (rank < n_parts) then

      ! essential attributes ...................................................

      mesh % n_bound     =  6
      mesh % n_parts     =  n_parts
      mesh % structured  =  .true.
      mesh % regular     =  .true.
      mesh % comm_world  =  comm
      mesh % proc        =  rank
      mesh % part        =  rank
      mesh % dx          =  dx

      if (present(pg)) then
        mesh % p_geom = pg
      else
        mesh % p_geom = 1
      end if

      ! dimensions .............................................................

      ! local mesh dimensions
      n1 = ep(1)
      n2 = ep(2)
      n3 = ep(3)

      ! indicator wether linked to itself
      self = periodic .and. np == 1

      mesh % n_vert = NumberOfVertices (n1, n2, n3, self)
      mesh % n_edge = NumberOfEdges    (n1, n2, n3, self)
      mesh % n_face = NumberOfFaces    (n1, n2, n3, self)
      mesh % n_elem = NumberOfElements (n1, n2, n3)

      ! structured mesh dimensions
      mesh % n_elem_1 = n1
      mesh % n_elem_2 = n2
      mesh % n_elem_3 = n3

      mesh % n_elem_active = mesh % n_elem
      mesh % n_elem_frozen = 0

      ! elements, faces and boundaries .........................................

      call GenerateRegularElements(mesh, np, ep, periodic, self, i0, j0, k0)
      call GenerateRegularElementGeometry(mesh, xo, i0, j0, k0)
      call GenerateRegularFaces(mesh, self)
      call GenerateRegularMeshBoundaries(mesh, periodic)

    end if

    !---------------------------------------------------------------------------
    ! comms, maps and empty partitions

    if (rank == 0) then
      attrib = MeshAttributes_3D(mesh)
    end if

    call attrib % Bcast(0, comm)

    if (rank >= n_parts) then
      mesh = Mesh_3D(attrib, comm)
    end if

    call mesh % BuildCommunicator()

    !---------------------------------------------------------------------------
    ! active partitions: remaining components

    if (rank < n_parts) then

      call mesh % BuildLinks()
      call mesh % BuildGhosts()
      call mesh % IdentifyRanks()

      allocate(mesh % map_child(0))
      allocate(mesh % map_parent(0))

    end if

  end subroutine GenerateRegularMesh

  !-----------------------------------------------------------------------------
  !> Builds the local mesh elements

  subroutine GenerateRegularElements(mesh, np, ep, periodic, self, i0, j0, k0)

    class(Mesh_3D), intent(inout) :: mesh  !< local partition
    integer, intent(in)  :: np(3)       !< num partitions in directions 1:3
    integer, intent(in)  :: ep(3)       !< num elements per partition and dir.
    logical, intent(in)  :: periodic(3) !< set true for periodic directions
    logical, intent(in)  :: self(3)     !< indicator wether linked to itself
    integer, intent(out) :: i0, j0, k0  !< element offsets WRT global numbering

    integer :: ie, je, ke         ! local element indices
    integer :: ip, jp, kp         ! partion triple index
    integer :: ies, jes, kes      ! shifted element indices
    integer :: ips, jps, kps      ! shifted partition indices
    integer :: i, e, l, r, s, t

    type(MeshElementNeighbor_3D) :: neighbor(26)

    ! prerequisites ............................................................

    call TripleIndex(ip, jp, kp, 1, 1, 1, np(1), np(2), l=mesh%part+1)

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

        element(e) % id = e

        ! element vertices .....................................................

        element(e) % vertex(1) % id = LexicalVertexIndex(ie-1,je-1,ke-1, n1,n2,n3, self)
        element(e) % vertex(2) % id = LexicalVertexIndex(ie  ,je-1,ke-1, n1,n2,n3, self)
        element(e) % vertex(3) % id = LexicalVertexIndex(ie-1,je  ,ke-1, n1,n2,n3, self)
        element(e) % vertex(4) % id = LexicalVertexIndex(ie  ,je  ,ke-1, n1,n2,n3, self)
        element(e) % vertex(5) % id = LexicalVertexIndex(ie-1,je-1,ke  , n1,n2,n3, self)
        element(e) % vertex(6) % id = LexicalVertexIndex(ie  ,je-1,ke  , n1,n2,n3, self)
        element(e) % vertex(7) % id = LexicalVertexIndex(ie-1,je  ,ke  , n1,n2,n3, self)
        element(e) % vertex(8) % id = LexicalVertexIndex(ie  ,je  ,ke  , n1,n2,n3, self)

        ! element edges ........................................................

        element(e) % edge( 1) % id = LexicalEdgeIndex(ie,je-1,ke-1, 1, n1,n2,n3, self)
        element(e) % edge( 2) % id = LexicalEdgeIndex(ie,je  ,ke-1, 1, n1,n2,n3, self)
        element(e) % edge( 3) % id = LexicalEdgeIndex(ie,je-1,ke  , 1, n1,n2,n3, self)
        element(e) % edge( 4) % id = LexicalEdgeIndex(ie,je  ,ke  , 1, n1,n2,n3, self)

        element(e) % edge( 5) % id = LexicalEdgeIndex(ie-1,je,ke-1, 2, n1,n2,n3, self)
        element(e) % edge( 6) % id = LexicalEdgeIndex(ie  ,je,ke-1, 2, n1,n2,n3, self)
        element(e) % edge( 7) % id = LexicalEdgeIndex(ie-1,je,ke  , 2, n1,n2,n3, self)
        element(e) % edge( 8) % id = LexicalEdgeIndex(ie  ,je,ke  , 2, n1,n2,n3, self)

        element(e) % edge( 9) % id = LexicalEdgeIndex(ie-1,je-1,ke, 3, n1,n2,n3, self)
        element(e) % edge(10) % id = LexicalEdgeIndex(ie  ,je-1,ke, 3, n1,n2,n3, self)
        element(e) % edge(11) % id = LexicalEdgeIndex(ie-1,je  ,ke, 3, n1,n2,n3, self)
        element(e) % edge(12) % id = LexicalEdgeIndex(ie  ,je  ,ke, 3, n1,n2,n3, self)

        ! element faces ........................................................

        element(e) % face(1) % id = LexicalFaceIndex(ie-1,je,ke, 1, n1,n2,n3, self)
        element(e) % face(2) % id = LexicalFaceIndex(ie  ,je,ke, 1, n1,n2,n3, self)

        element(e) % face(3) % id = LexicalFaceIndex(ie,je-1,ke, 2, n1,n2,n3, self)
        element(e) % face(4) % id = LexicalFaceIndex(ie,je  ,ke, 2, n1,n2,n3, self)

        element(e) % face(5) % id = LexicalFaceIndex(ie,je,ke-1, 3, n1,n2,n3, self)
        element(e) % face(6) % id = LexicalFaceIndex(ie,je,ke  , 3, n1,n2,n3, self)

        if (ie == 1  .and. ip == 1    )  element(e) % face(1) % boundary = 1
        if (ie == n1 .and. ip == np(1))  element(e) % face(2) % boundary = 2
        if (je == 1  .and. jp == 1    )  element(e) % face(3) % boundary = 3
        if (je == n2 .and. jp == np(2))  element(e) % face(4) % boundary = 4
        if (ke == 1  .and. kp == 1    )  element(e) % face(5) % boundary = 5
        if (ke == n3 .and. kp == np(3))  element(e) % face(6) % boundary = 6

        ! neighbors and boundaries .............................................

        neighbor = MeshElementNeighbor_3D()

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
          neighbor(i) % id = LexicalElementIndex(ies, jes, kes, n1, n2)

          ! neighbor element home partition ID
          neighbor(i) % part = LexicalIndex(ips,jps,kps,1,1,1,np(1),np(2),np(3)) - 1

          ! neighbor element coupled component and orientation
          neighbor(i) % component = int(ElementComponentIndex(-r,-s,-t), IXS)
          neighbor(i) % orientation = 12 ! always aligned

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

  subroutine GenerateRegularFaces(mesh, self)

    class(Mesh_3D), intent(inout) :: mesh  !< local partition
    logical, intent(in) :: self(3) !< indicator wether linked to itself

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

        l = LexicalFaceIndex(i, j, k, 1, n1, n2, n3, self)

        if (i > 0 .and. i < n1) then
          face(l) % element(1) % id   = LexicalElementIndex(i, j, k, n1, n2)
          face(l) % element(1) % face = 2
          face(l) % element(2) % id   = LexicalElementIndex(i+1, j, k, n1, n2)
          face(l) % element(2) % face = 1
        else if (i == 0) then
          if (self(1)) then
            face(l) % element(1) % id   = LexicalElementIndex(n1, j, k, n1, n2)
            face(l) % element(1) % face = 2
          end if
          face(l) % element(2) % id   = LexicalElementIndex(i+1, j, k, n1, n2)
          face(l) % element(2) % face = 1
        else if (i == n1) then
          face(l) % element(1) % id   = LexicalElementIndex(i, j, k, n1, n2)
          face(l) % element(1) % face = 2
          if (self(1)) then
            face(l) % element(2) % id   = LexicalElementIndex(1, j, k, n1, n2)
            face(l) % element(2) % face = 1
          end if
        end if

      end do
      end do
      end do

      ! x2-faces ...............................................................

      do k = 1, n3
      do j = 0, n2
      do i = 1, n1

        l = LexicalFaceIndex(i, j, k, 2, n1, n2, n3, self)

        if (j > 0 .and. j < n2) then
          face(l) % element(1) % id   = LexicalElementIndex(i, j, k, n1, n2)
          face(l) % element(1) % face = 4
          face(l) % element(2) % id   = LexicalElementIndex(i, j+1, k, n1, n2)
          face(l) % element(2) % face = 3
        else if (j == 0) then
          if (self(2)) then
            face(l) % element(1) % id   = LexicalElementIndex(i, n2, k, n1, n2)
            face(l) % element(1) % face = 4
          end if
          face(l) % element(2) % id   = LexicalElementIndex(i, j+1, k, n1, n2)
          face(l) % element(2) % face = 3
        else if (j == n2) then
          face(l) % element(1) % id   = LexicalElementIndex(i, j, k, n1, n2)
          face(l) % element(1) % face = 4
          if (self(2)) then
            face(l) % element(2) % id   = LexicalElementIndex(i, 1, k, n1, n2)
            face(l) % element(2) % face = 3
          end if
        end if

      end do
      end do
      end do

      ! x3-faces ...............................................................

      do k = 0, n3
      do j = 1, n2
      do i = 1, n1

        l = LexicalFaceIndex(i, j, k, 3, n1, n2, n3, self)

        if (k > 0 .and. k < n3) then
          face(l) % element(1) % id   = LexicalElementIndex(i, j, k, n1, n2)
          face(l) % element(1) % face = 6
          face(l) % element(2) % id   = LexicalElementIndex(i, j, k+1, n1, n2)
          face(l) % element(2) % face = 5
        else if (k == 0) then
          if (self(3)) then
            face(l) % element(1) % id   = LexicalElementIndex(i, j, n3, n1, n2)
            face(l) % element(1) % face = 6
          end if
          face(l) % element(2) % id   = LexicalElementIndex(i, j, k+1, n1, n2)
          face(l) % element(2) % face = 5
        else if (k == n3) then
          face(l) % element(1) % id   = LexicalElementIndex(i, j, k, n1, n2)
          face(l) % element(1) % face = 6
          if (self(3)) then
            face(l) % element(2) % id   = LexicalElementIndex(i, j, 1, n1, n2)
            face(l) % element(2) % face = 5
          end if
        end if

      end do
      end do
      end do

    end associate

  end subroutine GenerateRegularFaces

  !-----------------------------------------------------------------------------
  !> Creates the boundaries of a structured mesh partition

  subroutine GenerateRegularMeshBoundaries(mesh, periodic)

    class(Mesh_3D), intent(inout) :: mesh !< local partition
    logical, intent(in) :: periodic(3) !< indicator of periodic directions

    integer   :: b, e, f, i, j, k
    integer   :: coupled(6), polarity(6)
    real(RNP) :: map(4,4,6)
    character(len=6) :: name(6)

    ! preliminaries ............................................................

    name(1) = 'west'
    name(2) = 'east'
    name(3) = 'south'
    name(4) = 'north'
    name(5) = 'bottom'
    name(6) = 'top'

    coupled  = 0
    polarity = 0

    do b = 1, 6
      map(:,:,b) = AFFINE_IDENTITY_MAP_3D
    end do

    if (periodic(1)) then
      coupled (1) =  2
      polarity(1) = -1
      map (1,4,1) =  mesh%n_elem_1 * mesh%dx(1)
      coupled (2) =  1
      polarity(2) =  1
      map (1,4,2) = -mesh%n_elem_1 * mesh%dx(1)
    end if
    if (periodic(2)) then
      coupled (3) =  4
      polarity(3) = -2
      map (2,4,3) =  mesh%n_elem_2 * mesh%dx(2)
      coupled (4) =  3
      polarity(4) =  2
      map (2,4,4) = -mesh%n_elem_2 * mesh%dx(2)
    end if
    if (periodic(3)) then
      coupled (5) =  6
      polarity(5) = -3
      map (3,4,5) =  mesh%n_elem_3 * mesh%dx(3)
      coupled (6) =  5
      polarity(6) =  3
      map (3,4,6) = -mesh%n_elem_3 * mesh%dx(3)
    end if

    allocate(mesh % boundary(6))

    ! boundary attributes ......................................................

    do b = 1, 6
      mesh % boundary(b) = MeshBoundary_3D( name     = name(b)        &
                                          , id       = b              &
                                          , coupled  = coupled(b)     &
                                          , polarity = polarity(b)    &
                                          , map      = map(1:4,1:4,b) )
    end do

    if (mesh%part < 0) return

    ! faces ....................................................................

    associate( boundary => mesh % boundary  &
             , element  => mesh % element   &
             , n1       => mesh % n_elem_1  &
             , n2       => mesh % n_elem_2  &
             , n3       => mesh % n_elem_3  )

      ! west ...................................................................

      b = 1
      e = LexicalElementIndex(1, 1, 1, n1, n2)
      if (element(e) % face(b) % boundary == b) then
        boundary(b) % n_face = n2*n3
        allocate(boundary(b) % face(boundary(b) % n_face))
        i = 1
        f = 1
        do k = 1, n3
        do j = 1, n2
          e = LexicalElementIndex(i, j, k, n1, n2)
          boundary(b) % face(f) % element_id   = e
          boundary(b) % face(f) % element_face = b
          f = f + 1
        end do
        end do
      end if

      ! east ...................................................................

      b = 2
      e = LexicalElementIndex(n1, 1, 1, n1, n2)
      if (element(e) % face(b) % boundary == b) then
        boundary(b) % n_face = n2*n3
        allocate(boundary(b) % face(boundary(b) % n_face))
        i = n1
        f = 1
        do k = 1, n3
        do j = 1, n2
          e = LexicalElementIndex(i, j, k, n1, n2)
          boundary(b) % face(f) % element_id   = e
          boundary(b) % face(f) % element_face = b
          f = f + 1
        end do
        end do
      end if

      ! south ..................................................................

      b = 3
      e = LexicalElementIndex(1, 1, 1, n1, n2)
      if (element(e) % face(b) % boundary == b) then
        boundary(b) % n_face = n1*n3
        allocate(boundary(b) % face(boundary(b) % n_face))
        j = 1
        f = 1
        do k = 1, n3
        do i = 1, n1
          e = LexicalElementIndex(i, j, k, n1, n2)
          boundary(b) % face(f) % element_id   = e
          boundary(b) % face(f) % element_face = b
          f = f + 1
        end do
        end do
      end if

      ! north ..................................................................

      b = 4
      e = LexicalElementIndex(1, n2, 1, n1, n2)
      if (element(e) % face(b) % boundary == b) then
        boundary(b) % n_face = n1*n3
        allocate(boundary(b) % face(boundary(b) % n_face))
        j = n2
        f = 1
        do k = 1, n3
        do i = 1, n1
          e = LexicalElementIndex(i, j, k, n1, n2)
          boundary(b) % face(f) % element_id   = e
          boundary(b) % face(f) % element_face = b
          f = f + 1
        end do
        end do
      end if

      ! bottom .................................................................

      b = 5
      e = LexicalElementIndex(1, 1, 1, n1, n2)
      if (element(e) % face(b) % boundary == b) then
        boundary(b) % n_face = n1*n2
        allocate(boundary(b) % face(boundary(b) % n_face))
        k = 1
        f = 1
        do j = 1, n2
        do i = 1, n1
          e = LexicalElementIndex(i, j, k, n1, n2)
          boundary(b) % face(f) % element_id   = e
          boundary(b) % face(f) % element_face = b
          f = f + 1
        end do
        end do
      end if

      ! top ....................................................................

      b = 6
      e = LexicalElementIndex(1, 1, n3, n1, n2)
      if (element(e) % face(b) % boundary == b) then
        boundary(b) % n_face = n1*n2
        allocate(boundary(b) % face(boundary(b) % n_face))
        k = n3
        f = 1
        do j = 1, n2
        do i = 1, n1
          e = LexicalElementIndex(i, j, k, n1, n2)
          boundary(b) % face(f) % element_id   = e
          boundary(b) % face(f) % element_face = b
          f = f + 1
        end do
        end do
      end if

    end associate

  end subroutine GenerateRegularMeshBoundaries

  !-----------------------------------------------------------------------------
  !> Creates element domains, cuboid approximations and mean face-normal spacing

  subroutine GenerateRegularElementGeometry(mesh, xo, i0, j0, k0)
    class(Mesh_3D), intent(inout) :: mesh  !< local partition
    real(RNP), intent(in) :: xo(3)      !< corner closest to -∞
    integer  , intent(in) :: i0, j0, k0 !< element offsets WRT global numbering

    real(RNP), dimension(0 : mesh%p_geom) :: x1, x2, x3, ys
    integer :: e, i, j, k, r, s, t

    associate(po => mesh % p_geom, dx => mesh % dx)

      ! Lobatto points transformed to [-1,0]
      ys = (LobattoPoints(po) - 1) / 2

      do k = 1, mesh % n_elem_3
      do j = 1, mesh % n_elem_2
      do i = 1, mesh % n_elem_1

        e = LexicalElementIndex(i, j, k, mesh % n_elem_1, mesh % n_elem_2)

        associate(geometry => mesh % element(e) % geometry)

          ! degree
          geometry % po = po

          ! 1D point distributions
          x1 = xo(1) + (i0 + i + ys) * dx(1)
          x2 = xo(2) + (j0 + j + ys) * dx(2)
          x3 = xo(3) + (k0 + k + ys) * dx(3)

          ! element points
          allocate(geometry % x_e(0:po, 0:po, 0:po, 1:3))
          do t = 0, po
          do s = 0, po
          do r = 0, po
            geometry % x_e(r,s,t,1) = x1(r)
            geometry % x_e(r,s,t,2) = x2(s)
            geometry % x_e(r,s,t,3) = x3(t)
          end do
          end do
          end do

          ! cuboid
          geometry % x_c(0:3,1) = HALF * [ x1(0) + x1(po), dx(1), ZERO , ZERO  ]
          geometry % x_c(0:3,2) = HALF * [ x2(0) + x2(po), ZERO , dx(2), ZERO  ]
          geometry % x_c(0:3,3) = HALF * [ x3(0) + x3(po), ZERO , ZERO , dx(3) ]

          ! mean spacing normal to faces
          geometry % dx_m = dx([1,1,2,2,3,3])

        end associate

      end do
      end do
      end do

    end associate

  end subroutine GenerateRegularElementGeometry

  !=============================================================================

end module Generate_Regular_Mesh__3D
