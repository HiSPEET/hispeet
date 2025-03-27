!> summary:  Generate a generic 3d mesh from Gmsh meshfile
!> author:   Matthias Frey, Benedikt Wex, Joerg Stiller
!> date:     2023/06/16
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> @note
!>   - Expects input file compatible to GMSH format 4.1
!>   - Restricted to hexahedral elements
!>   - All boundaries must be provided with a physical tag
!===============================================================================

module Import_GMSH__3D
  use Kind_Parameters
  use Constants
  use Execution_Control
  use Generic_Mesh__3D
  implicit none
  private

  public :: ImportGMSH_3D

  !---------------------------------------------------------------------------
  !> Length of GMSH names

  integer, parameter :: MSH_LEN_NAME = 127

  !-----------------------------------------------------------------------------
  !> GMSH physical entity

  type MshPhysicalEntity
    integer                 :: dim           !< dimension
    integer                 :: physicalTag   !< tag
    character(MSH_LEN_NAME) :: physicalName  !< name
  end type MshPhysicalEntity

  !-----------------------------------------------------------------------------
  !> GMSH mesh node

  type MshNode
    integer   :: nodeTag      = 0 !< node tag
    real(RNP) :: nodeCoord(3) = 0 !< node coordinates
  end type MshNode

  !-----------------------------------------------------------------------------
  !> GMSH mesh element

  type MshElement
    integer :: id          = 0            !< element id
    integer :: elementTag  = 0            !< element tag
    integer :: vertices(8) = 0            !< element vertices
    integer, allocatable :: nodes(:)      !< 1D node vector
    integer, allocatable :: points(:,:,:) !< 3D point array
  end type MshElement

  !-----------------------------------------------------------------------------
  !> GMSH mesh boundary

  type MshBoundary
    character(MSH_LEN_NAME) :: physicalName !< boundary name
    integer :: physicalTag = 0 !< GMSH physical tag
    integer :: nFaces      = 0 !< num boundary faces
    integer :: coupledID   = 0 !< numeration ID of coupled surface
  end type MshBoundary

  !-----------------------------------------------------------------------------
  !> GMSH mesh surface
  !>
  !> Surfaces are entities composed of 2D mesh faces. Each boundary consists of
  !> one ore more surfaces. However, not surfaces are not required to be part of
  !> a boundary.

  type MshSurface
    integer :: boundaryID  = 0       !< numeration ID of corresponding boundary
    integer :: physicalTag = 0       !< tag of physical GMSH entity
    integer :: surfaceTag  = 0       !< tag of GMSH surface entity
    integer :: periodicTag = 0       !< tag of coupled surface
    integer :: periodicID  = 0       !< numeration ID of coupled surface
    integer :: nFaces      = 0       !< number of corresponding faces
    integer :: nNodes      = 0       !< number of corresponding nodes
    logical :: periodic    = .FALSE. !< flag for periodic furfaces
    logical :: master      = .FALSE. !< flag for master furfaces
    integer, allocatable :: nodes(:) !< tags of surface nodes
    real(RNP) :: A(4,4)              !< affinity transform to coupled surface
  end type MshSurface

  !-----------------------------------------------------------------------------
  !> GMSH mesh boundary face, stored by GMSH as 2D element

  type MshBoundaryFace
    integer :: id          = 0 !< numeration ID
    integer :: surfaceID   = 0 !< numeration ID of corresponding surface
    integer :: surfaceTag  = 0 !< tag of corresponding surface
    integer :: physicalTag = 0 !< corresponding physical tags
    integer :: nodes(4)    = 0 !< 4 corner nodes
    integer :: elem_id     = 0 !< element tag
    integer :: elem_face   = 0 !< element face
  end type MshBoundaryFace

  !=============================================================================
  ! External module procedures

  interface

    !---------------------------------------------------------------------------
    !> Maps shell vertex, edge and surface points to collocation points

    module subroutine MapShellNodes(v, e, s, shellOrder, points)
      integer, intent(in)  :: shellOrder !< shell order
      integer, intent(in)  :: v(:)       !< vertex nodes
      integer, intent(in)  :: e(:)       !< edge nodes without end nodes
      integer, intent(in)  :: s(:)       !< surface nodes without boundary
      integer, allocatable, intent(out) :: points(:,:,:) !< collocation points
    end subroutine MapShellNodes

    !---------------------------------------------------------------------------
    !> Merges the IDs of periodic vertices

    module subroutine MergePeriodicVertices(mesh, surface)
      class(GenericMesh_3D), intent(inout) :: mesh       !< generic mesh
      type(MshSurface),      intent(in)    :: surface(:) !< GMSH surfaces
    end subroutine MergePeriodicVertices

  end interface

contains

  !-----------------------------------------------------------------------------
  !> Reads the mesh from file and converts it to a GenericMesh_3D object

  subroutine ImportGMSH_3D(mshfile, mesh)
    character(len=*),      intent(in)  :: mshfile !< GMSH file (*.msh)
    class(GenericMesh_3D), intent(out) :: mesh    !< output generic mesh

    type(MshPhysicalEntity), allocatable :: entity(:)
    type(MshNode),           allocatable :: node(:)
    type(MshElement),        allocatable :: element(:)
    type(MshBoundary),       allocatable :: boundary(:)
    type(MshBoundaryFace),   allocatable :: bface(:)
    type(MshSurface),        allocatable :: surface(:)

    integer, allocatable :: node_tag(:)    ! map from node tags to node IDs
    integer, allocatable :: vertex_node(:) ! map from node IDs to vertices
    integer, allocatable :: points(:,:,:)  ! collocation nodes
    integer, allocatable :: p_shell(:)     ! list of shell nodes
    integer, allocatable :: p_edge(:)      ! list of edge nodes
    integer, allocatable :: p_face(:)      ! list of face nodes
    logical, allocatable :: mask(:)        ! array for masking nodes or vertices

    real(RNP):: affinityMatrix(4,4) ! affinity transformation matrix

    integer :: MSH, stat
    integer :: e, f, i, j, k, l, m, n, param, pT, s, sT
    integer :: numPoints, numCurves, numSurfaces, numVolumes
    integer :: numNodes, numElements, numVertices
    integer :: numBoundaries, numBoundaryFaces, numHexElements
    integer :: numEntityBlocks, entityDim, entityTag
    integer :: numPhysicalNames, numPhysicalTags
    integer :: numNodesInBlock, numElementsInBlock
    integer :: numPeriodicLinks, entityTagMaster, nodeTagMaster
    integer :: numCorrespondingNodes, numAffine, nodeTag, maxNodeTag
    integer :: elementType, meshOrder, numShells, shellOrder
    integer :: edgeCenter, np
    integer :: p_vert(8)
    logical :: match(8)
    character(MSH_LEN_NAME) :: mshblock
    real :: minX, minY, minZ, maxX, maxY, maxZ

    !---------------------------------------------------------------------------
    ! Read Gmsh mesh file

    ! open GMSH file
    write(*,'(2X,A)') 'opening GMSH file: ' // trim(mshfile)

    open(newunit = MSH, file = trim(mshfile)//'.msh')

    ! skip first three lines
    read(MSH,'(2/)')

    ! check if boundaries have been specified properly
    read(MSH,*) mshblock
    if(mshblock /= '$PhysicalNames') then
      call Error('ImportGMSH_3D', 'No boundaries specified', 'Import_GMSH__3D')
    end if

    ! read $PhysicalNames block ................................................

    read(MSH,*) numPhysicalNames

    allocate(entity(numPhysicalNames))
    do i = 1, numPhysicalNames
      read(MSH,*) entity(i)%dim, entity(i)%physicalTag, entity(i)%physicalName
    end do
    read(MSH,'(/)')

    ! NOTE: assuming that only boundary surfaces are listed here
    numBoundaries = count(entity%dim == 2)

    allocate(boundary(numBoundaries))
    j = 0
    do i = 1, numPhysicalNames
      if (entity(i)%dim == 2) then
        j = j + 1
        boundary(j) % physicalName = entity(i) % physicalName
        boundary(j) % physicalTag  = entity(i) % physicalTag
      end if
    end do

    ! read $Entities block .....................................................

    read(MSH,*) numPoints, numCurves, numSurfaces, numVolumes

    ! ignore points & curves
    do i = 1, numPoints + numCurves
      read(MSH,*)
    end do

    ! identify mesh boundary surfaces (without duplicates)
    allocate(surface(numSurfaces))
    j = 0
    do i = 1, numSurfaces
      read(MSH,*) sT, minX, minY, minZ, maxX, maxY, maxZ, numPhysicalTags, pT
      ! ignore duplicates or internal surfaces
      if (numPhysicalTags == 0) cycle
      j = j + 1
      surface(j) % physicalTag = pT
      surface(j) % surfaceTag  = sT
    end do
    numSurfaces = j
    surface = surface(1:numSurfaces)

    ! ignore volumes
    do i = 1, numVolumes
      read(MSH,*)
    end do

    ! read $Nodes block ........................................................

    write(*,'(2X,A)') 'reading node data'
    read(MSH,'(1/)')
    read(MSH,*) numEntityBlocks, numNodes
    allocate(node(numNodes))

    maxNodeTag = 0
    k = 0
    do i = 1, numEntityBlocks
      read(MSH,*) entityDim, entityTag, param, numNodesInBlock
      if (numNodesInBlock == 0) cycle
      do j = 1, numNodesInBlock
        k = k + 1
        read(MSH,*) node(k)%nodeTag
        maxNodeTag = max(node(k)%nodeTag, maxNodeTag)
      end do
      k = k - numNodesInBlock
      do j = 1, numNodesInBlock
        k = k + 1
        read(MSH,*) node(k)%nodeCoord
      end do
    end do

    ! map tags to nodes
    allocate(node_tag(maxNodeTag), source = -1)
    do k = 1, numNodes
      node_tag(node(k)%nodeTag) = k
    end do

    ! read $Elements block .....................................................

    write(*,'(2X,A)') 'reading boundary faces and elements'

    read(MSH,'(/)')
    read(MSH,*) numEntityBlocks, numElements

    allocate(mask(numNodes))
    allocate(element(numElements))
    allocate(bface(numElements))

    f = 0 ! boundary face counter
    e = 0 ! hexahedral element counter
    s = 0 ! surface counter

    do i = 1, numEntityBlocks
      read(MSH,*) entityDim, entityTag, elementType, numElementsInBlock

      select case(entityDim)

      case(2)
        ! surface
        s = s + 1
        if (entityTag /= surface(s)%surfaceTag) then
          call Error( 'ImportGMSH_3D'                        &
                    , 'surface tag does not match entityTag' &
                    , 'Import_GMSH__3D'                      )
        end if
        allocate(surface(s)%nodes(numElementsInBlock*4), source = -1)
        surface(s) % nFaces = numElementsInBlock
        mask = .false.
        l = 0
        do j = 1, numElementsInBlock
          f = f + 1
          bface(f)%id         = f
          bface(f)%surfaceTag = entityTag
          bface(f)%surfaceID  = s
          read(MSH,*) k, bface(f)%nodes
          do k = 1, 4
            n = bface(f)%nodes(k)
            if (mask(n)) cycle ! node already added
            l = l + 1
            surface(s)%nodes(l) = n
            mask(n) = .true.
          end do
        end do
        surface(s) % nNodes = l

      case(3)
        ! volume
        do j = 1, numElementsInBlock
          e = e + 1
          element(e)%id = e
          meshOrder = GetMeshOrder(elementType)
          allocate(element(e)%nodes((meshOrder+1)**3))
          read(MSH,*) element(e)%elementTag, element(e)%nodes
        end do

      end select
    end do

    numBoundaryFaces = f
    numHexElements   = e

    ! read $Periodic block .....................................................

    read(MSH,*)
    read(MSH,*, iostat=stat)

    if(stat == 0) then
      write(*,'(2X,A)') 'reading periodic boundary data'
      read(MSH,*) numPeriodicLinks

      do i = 1, numPeriodicLinks
        read(MSH,*) entityDim, entityTag, entityTagMaster
        read(MSH,*) numAffine, affinityMatrix ! affinity matrix with C ordering
        read(MSH,*) numCorrespondingNodes
        do j = 1, numCorrespondingNodes       ! corresponding
          read(MSH,*) nodeTag, nodeTagMaster  ! node tags
        end do                                ! are ignored

        if (entityDim == 2) then

          do j = 1, numSurfaces
            if (surface(j) % surfaceTag == entityTagMaster) m = j
            if (surface(j) % surfaceTag == entityTag      ) s = j
          end do

          ! master surface
          surface(m) % periodicTag    = entityTag
          surface(m) % periodicID     = s
          surface(m) % periodic       = .TRUE.
          surface(m) % master         = .TRUE.
          surface(m) % A              = transpose(affinityMatrix)

          ! slave surface
          surface(s) % periodicTag    = entityTagMaster
          surface(s) % periodicID     = m
          surface(s) % periodic       = .TRUE.
          surface(s) % master         = .FALSE.

        end if
      end do
    else
      write(*,'(4X,A)') 'no periodic boundaries found'
      numPeriodicLinks = 0
    end if

    close(MSH)

    !--------------------------------------------------------------------------
    ! Process Gmsh mesh data

    ! perform affine spatial transformations to find periodic node pairs
    !   - identification/numbering of vertices
    !   - identification of element points
    !   - identification of boundaries

    ! allocate high order points in 3D array for each element .................

    write(*,'(2X,A)') 'identifying collocation points'

    np         = meshOrder+1         ! number of nodes per direction
    numShells  = floor(meshOrder/2.) ! number of shells
    edgeCenter = meshOrder/2+1       ! position of center point on odd edge

    ! map element shells into 3D array points
    do i = 1, numHexElements
      allocate(element(i)%points(np,np,np), source = 1)

      element(i)%vertices = element(i)%nodes(1:8)

      if(meshOrder > 1) then
        ! subdivide nodes into shells
        shellOrder = meshOrder
        k = 1
        l = 0
        do j = 1, numShells
          allocate(points(np-l,np-l,np-l), source = 1)
          allocate(p_shell(NumShellNodes(shellOrder)))

          p_shell = element(i)%nodes(k:k+NumShellNodes(shellOrder)-1)

          ! subdivide nodes into vertex, edge and face nodes
          allocate(p_edge((shellOrder-1)*12), source = 1)
          allocate(p_face(((shellOrder-1)**2)*6), source = 1)

          p_vert = p_shell(1:8)
          p_edge = p_shell(9:8+size(p_edge))
          p_face = p_shell(size(p_edge)+9:size(p_shell))

          ! map the shell nodes into 3D array of collocation points
          call MapShellNodes(p_vert,p_edge,p_face,shellOrder,points)
          element(i)%points(j:np-(j-1),j:np-(j-1),j:np-(j-1)) = points

          ! prepare variables for next shell with order reduced by 2
          k = k + NumShellNodes(shellOrder)
          l = l + 2
          shellOrder = shellOrder - 2
          deallocate(p_shell, p_edge, p_face, points)
        end do
      end if

      ! map points in the element center
      associate(nodes => element(i)%nodes, points => element(i)%points)

        ! check if mesh order is odd or even
        select case(mod(meshOrder,2))
        case(0) ! even -> one central point
          points(edgeCenter,edgeCenter,edgeCenter) = nodes(size(nodes))

        case(1) ! odd  -> 8 corner points
          points(np/2+1, np/2  , np/2)   = nodes(size(nodes) - 7)
          points(np/2+1, np/2+1, np/2)   = nodes(size(nodes) - 6)
          points(np/2+1, np/2+1, np/2+1) = nodes(size(nodes) - 5)
          points(np/2+1, np/2  , np/2+1) = nodes(size(nodes) - 4)
          points(np/2  , np/2  , np/2)   = nodes(size(nodes) - 3)
          points(np/2  , np/2+1, np/2)   = nodes(size(nodes) - 2)
          points(np/2  , np/2+1, np/2+1) = nodes(size(nodes) - 1)
          points(np/2  , np/2  , np/2+1) = nodes(size(nodes)    )
        end select

      end associate

    end do

    ! identify vertices ........................................................

    write(*,'(2X,A)') 'identifying vertices ...'
    allocate(vertex_node(numNodes), source = 0)

    ! mark element vertex nodes
    do i = 1, numHexElements
      vertex_node(element(i)%vertices) = 1
    end do

    ! only choose nodes, where the vertex mask is true
    k = 0
    do i = 1, numNodes
      if (vertex_node(i) == 0) cycle
      k = k + 1
      vertex_node(i) = k
    end do
    numVertices = k

    ! transform surface node entries to vertex IDs
    do j = 1, numSurfaces
    do l = 1, surface(j)%nNodes
      n = surface(j)%nodes(l)
      if (n > 0) then
        surface(j)%nodes(l) = vertex_node( node_tag(n) )
      end if
    end do
    end do

    ! map physical tags to boundary faces .....................................

    do i = 1, numBoundaryFaces
    do j = 1, numSurfaces
      if (bface(i)%surfaceID == j) then
        bface(i)%physicalTag = surface(j)%physicalTag
        exit
      end if
    end do
    end do

    ! identify elements corresponing to each boundary face .....................

    write(*,'(2X,A)') 'identifying boundary elements & faces'

    ! mark all boundary nodes
    mask = .false.
    do i = 1, numBoundaryFaces
    if (bface(i)%physicalTag == 0) cycle
      mask(bface(i)%nodes) = .true.
    end do

    ! identify elements adjacent to boundary faces
    do i = 1, numBoundaryFaces
      if (bface(i)%physicalTag == 0) cycle
      do j = 1, numHexElements

        ! skip elements with no boundary face
        if (count(mask(element(j)%vertices)) < 4) cycle

        ! identify element vertices matching a vertex of the boundary face
        match = .false.
        do k = 1, 8
          match(k) = any(element(j)%vertices(k) == bface(i)%nodes)
        end do

        if (count(match) < 4) then
          cycle
        else if (match(1) .and. match(2) .and. match(5)) then
          ! face 1 with element vertices [1,2,5,6]
          m = 1
        else if (match(2) .and. match(3) .and. match(6)) then
          ! face 2 with element vertices [2,3,6,7]
          m = 2
        else if (match(3) .and. match(4) .and. match(7)) then
          ! face 3 with element vertices [3,4,7,8]
          m = 3
        else if (match(1) .and. match(4) .and. match(5)) then
          ! face 4 with element vertices [1,4,5,8]
          m = 4
        else if (match(1) .and. match(2) .and. match(3)) then
          ! face 5 with element vertices [1,2,3,4]
          m = 5
        else if (match(5) .and. match(6) .and. match(7)) then
          ! face 6 with element vertices [5,6,7,8]
          m = 6
        else
          m = 0
        end if

        if (m > 0) then
          bface(i)%elem_id   = element(j)%id
          bface(i)%elem_face = m
          exit
        end if

      end do
    end do

    ! map surfaces for boundaries and count boundary faces .....................

    do i = 1, numBoundaries
      boundary(i)%nFaces = 0
      do j = 1, numSurfaces
        if (boundary(i)%physicalTag == surface(j)%physicalTag) then
          boundary(i)%nFaces = boundary(i)%nFaces + surface(j)%nFaces
          surface(j)%boundaryID = i
        end if
      end do
    end do

    ! identify periodic boundaries .............................................

    if (numPeriodicLinks > 0) then
      do i = 1, numBoundaries
        do j = 1, numSurfaces
          if (surface(j)%boundaryID /= i) cycle
          if (surface(j)%periodic) then
            k = surface(surface(j)%periodicID)%boundaryID
            boundary(i)%coupledID = k
            boundary(k)%coupledID = i
            exit
          end if
        end do
      end do
    end if

    ! give some mesh information as display output .............................

    write(*,'(/,2X,A)') 'mesh info'
    write(*,'(4X,A,T25,G0)') 'num vertices'       , numVertices
    write(*,'(4X,A,T25,G0)') 'num elements'       , numHexElements
    write(*,'(4X,A,T25,G0)') 'num boundaries'     , numBoundaries
    write(*,'(4X,A,T25,G0)') 'num boundary faces:', numBoundaryFaces
    write(*,'(4X,A,T25,G0)') 'geometry order:'    , meshOrder
    write(*,*)

    !---------------------------------------------------------------------------
    ! Create HiSPEET generic mesh

    ! create generic mesh vertices .............................................

    write(*,'(2X,A)') 'creating generic mesh vertices'

    allocate (mesh%vertex(numVertices))
    do i = 1, numNodes
      k = vertex_node(i)
      if (k == 0) cycle ! no vertex
        mesh%vertex(k)%id = k
        mesh%vertex(k)%x  = node(i)%nodeCoord
    end do

    write(*,'(2X,A)') 'merging periodic vertices'
    if (numPeriodicLinks > 0) then
      call MergePeriodicVertices(mesh, surface)
    end if

    ! create generic mesh elements ............................................

    write(*,'(2X,A)') 'creating generic mesh elements'

    mesh%numbering = ROTATIONAL_NUMBERING

    allocate(mesh%element(numHexElements))

    do i = 1, numHexElements

      ! ID, type, order and basis
      mesh%element(i)%id    = element(i)%id
      mesh%element(i)%typ   = HEXAHEDRAL_ELEMENT
      mesh%element(i)%order = meshOrder
      mesh%element(i)%basis = EQUIDISTANT_NODAL_BASIS

      ! vertices
      do k = 1, 8
        mesh%element(i)%vertex(k) = vertex_node(element(i)%vertices(k))
      end do

      ! collocation points
      allocate(mesh%element(i)%x(np*np*np, 3))
      associate(x => mesh%element(i)%x, point => element(i)%points)
        m = 0
        do j = np, 1, -1
          do l = 1, np
            do k = 1, np
              m = m + 1
              n = node_tag(point(j,k,l))
              x(m,1) = node(n) % nodeCoord(1)
              x(m,2) = node(n) % nodeCoord(2)
              x(m,3) = node(n) % nodeCoord(3)
            end do
          end do
        end do
      end associate

    end do

    ! create generic mesh boundaries ..........................................

    write(*,'(2X,A)') 'creating generic mesh boundaries'

    allocate(mesh%boundary(numBoundaries))

    do i = 1, numBoundaries
      mesh%boundary(i)%id      = i
      mesh%boundary(i)%name    = boundary(i)%physicalName
      mesh%boundary(i)%coupled = boundary(i)%coupledID
      allocate(mesh%boundary(i)%face(boundary(i)%nFaces))

      k = 0
      do j = 1, numBoundaryFaces
        if (bface(j)%physicalTag == boundary(i)%physicalTag) then
          k = k + 1
          mesh%boundary(i)%face(k)%element_id   = bface(j)%elem_id
          mesh%boundary(i)%face(k)%element_face = bface(j)%elem_face
        end if
      end do

      ! print boundary info
      write(*,'(I6,2X)',advance='NO') mesh%boundary(i)%id
      write(*,'(A, 2X)',advance='NO') mesh%boundary(i)%name(1:60)
      if (mesh%boundary(i)%coupled > 0) then
        write(*,'(A4,I5)',advance='NO') 'c:', mesh%boundary(i)%coupled
      end if
      write(*,*)

    end do

    ! enforce lexical numbering and generate periodic vertex IDs ..............

    call mesh % SwitchToLexicalNumbering()

  end subroutine ImportGMSH_3D

  !-----------------------------------------------------------------------------
  !> Returns mesh order according to element type numbers defined by GMSH

  integer function GetMeshOrder(elementType) result(meshOrder)
    integer, intent (in) :: elementType  !< element type number defined by GMSH

    select case(elementType)
    case(5)
      meshOrder = 1
    case(12)
      meshOrder = 2
    case(92:98)
      meshOrder = elementType - 89
    case default
      call Error('GetMeshOrder','Element type not supported','Import_GMSH__3D')
    end select

  end function GetMeshOrder

  !-----------------------------------------------------------------------------
  !> Returns the number of nodes per shell for a given shell order

  pure integer function NumShellNodes(shellOrder)
    integer, intent (in) :: shellOrder !< shell order

    NumShellNodes = (shellOrder + 1)**3 - (shellOrder-1)**3

  end function NumShellNodes

  !=============================================================================

end module Import_GMSH__3D
