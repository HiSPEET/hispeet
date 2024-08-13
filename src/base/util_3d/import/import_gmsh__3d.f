!> summary:  Generate a generic 3d mesh from Gmsh meshfile
!> author:   Matthias Frey, Joerg Stiller
!> date:     2023/06/16
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> @note
!> This version is restricted to hexahedral elements.
!===============================================================================

module Import_GMSH__3D
  use Kind_Parameters
  use Execution_Control
  use Generic_Mesh__3D

  implicit none
  private

  public :: ImportGMSH_3D

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

  end interface

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
    integer, allocatable :: faces(:,:)    !< map of vertex to face IDs
    integer, allocatable :: nodes(:)      !< 1D node vector
    integer, allocatable :: points(:,:,:) !< 3D point array
  end type MshElement

  !-----------------------------------------------------------------------------
  !>  GMSH mesh boundary

  type MshBoundary
    character(len = 127) :: physicalName !< boundary name
    integer :: physicalTag = 0           !< corresponding physical tags
    integer :: nElem       = 0           !< number of corresponding elements
  end type MshBoundary

  !-----------------------------------------------------------------------------
  !> GMSH mesh boundary face, stored by GMSH as 2D element

  type MshBoundaryFace
    integer :: id          = 0 !< numeration id
    integer :: surfaceID   = 0 !< numeration id of corresponding surface
    integer :: surfaceTag  = 0 !< tag of corresponding surface
    integer :: physicalTag = 0 !< corresponding physical tags
    integer :: nodes(4)    = 0 !< 4 corner nodes
    integer :: elem_id     = 0 !< element tag
    integer :: elem_face   = 0 !< element face
  end type MshBoundaryFace

contains

  !-----------------------------------------------------------------------------
  !> Reads the mesh from file and converts it to a GenericMesh_3D object

  subroutine ImportGMSH_3D(mshfile, mesh)
    character(len=*),      intent(in)  :: mshfile !< GMSH file (*.msh)
    class(GenericMesh_3D), intent(out) :: mesh    !< output generic mesh

    type(MshNode),         allocatable :: node(:)
    type(MshElement),      allocatable :: element(:)
    type(MshBoundary),     allocatable :: boundary(:)
    type(MshBoundaryFace), allocatable :: bface(:)

    integer :: i, io, j, k, l, m, MSH, n, param, pT, sT
    integer :: numPoints, numCurves, numSurfaces, numVolumes
    integer :: numNodes, numElements, numVertices
    integer :: numBoundaries, numBfaces, numHexElem
    integer :: numEntityBlocks, entityDim, entityTag
    integer :: numPhysicalNames, numPhysicalTags
    integer :: numNodesInBlock, numElementsInBlock
    integer :: numPeriodicLinks, entityTagMaster, nodeTagMaster
    integer :: numCorrespondingNodes, numAffine, nodeTag, maxNodeTag
    integer :: elementType, meshOrder, numShells, shellOrder
    integer :: edgeCenter, np
    integer :: mask(4), p_vert(8)
    real    :: minX, minY, minZ, maxX, maxY, maxZ
    character(len = 100) :: mshblock

    integer, allocatable :: node_tag(:)   ! mask to identify the node tags
    integer, allocatable :: vertex_mask(:)    ! mask to extract vertices from nodes

    integer, allocatable :: boundaryList(:,:) ! list of boundaries surface
    integer, allocatable :: periodicSurf(:,:) ! tags of periodical surfaces

    integer, allocatable :: points(:,:,:)      !< collocation nodes
    integer, allocatable :: p_shell(:)        !< list of shell nodes
    integer, allocatable :: p_edge(:)         !< list of edge nodes
    integer, allocatable :: p_face(:)         !< list of face nodes

    !---------------------------------------------------------------------------
    ! Read Gmsh mesh file

    ! open GMSH file
    write(*,'(/,A,/)') 'opening GMSH file ' // trim(mshfile)

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
    numBoundaries = numPhysicalNames - 1  ! ignore volume entity
    allocate(boundary(numBoundaries))

    do i = 1, numBoundaries
      read(MSH,*) entityDim, boundary(i)%physicalTag, boundary(i)%physicalName
    end do

    read(MSH,'(2/)')

    ! read $Entities block .....................................................

    read(MSH,*) numPoints, numCurves, numSurfaces, numVolumes
    numEntityBlocks = numPoints+numCurves+numSurfaces+numVolumes
    allocate(boundaryList(numSurfaces,4), source = 0)

    ! ignore points & curves
    do i = 1, numPoints + numCurves
      read(MSH,*)
    end do

    ! tags of mesh boundary surfaces
    j = 0
    do i = 1, numSurfaces
      read(MSH,*) sT, minX, minY, minZ, &
           maxX, maxY, maxZ, numPhysicalTags, pT
      if (numPhysicalTags == 0) cycle  ! surface i is not a boundary
      j = j + 1
      boundaryList(j,1:3) = [j, sT, pT]
    end do

    ! ignore volumes
    do i = 1, numVolumes
      read(MSH,*)
    end do

    ! read $Nodes block.........................................................

    write(*,'(2X,A)') 'reading node data'
    read(MSH,'(1/)')
    read(MSH,*) numEntityBlocks, numNodes
    allocate(node(numNodes))

    maxNodeTag = 0
    k = 0
    do i = 1, numEntityBlocks
      read(MSH,*) entityDim, entityTag, param, numNodesInBlock
      if (numNodesInBlock == 0) then
        cycle
      else
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
      end if
    end do

    ! map tag to each node
    allocate(node_tag(maxNodeTag), source = -1)
    do k = 1, numNodes
      node_tag(node(k)%nodeTag) = k
    end do

    ! read $Elements block .....................................................

    write(*,'(2X,A)') 'reading boundary and element data'

    read(MSH,'(/)')
    read(MSH,*) numEntityBlocks, numElements

    allocate(element(numElements))
    allocate(bface(numElements))

    numBfaces  = 0
    numHexElem = 0
    do i = 1, numEntityBlocks
      read(MSH,*) entityDim, entityTag, elementType, numElementsInBlock

      if (numElementsInBlock == 0) then
        cycle
      else
        do j = 1, numElementsInBlock
          select case(entityDim)

          ! 2D element coresponds to boundary face
          case(2)
            boundaryList(i,4) = numElementsInBlock ! num elements on boundary
            numBfaces = numBfaces + 1
            bface(numBfaces)%id = numBfaces
            bface(numBfaces)%surfaceTag = entityTag
            bface(numBfaces)%surfaceID = i
            read(MSH,*) k , bface(numBfaces)%nodes

          ! 3D hexaedral element
          case(3)
            numHexElem = numHexElem + 1
            element(numHexElem)%id = numHexElem
            meshOrder = GetMeshOrder(elementType)
            allocate(element(numHexElem)%nodes((meshOrder+1)**3))

            read(MSH,*) element(numHexElem)%elementTag, &
                        element(numHexElem)%nodes
          end select
        end do
      end if
    end do

    ! read $Periodic block .....................................................

    read(MSH,*)
    read(MSH,*,iostat=io)

    if (io == 0) then
      write(*,'(2X,A)') 'reading periodic boundary data'

      read(MSH,*) numPeriodicLinks
      allocate(periodicSurf(numPeriodicLinks,6), source = 0)

      do i = 1, numPeriodicLinks
        read(MSH,*) entityDim, entityTag, entityTagMaster
        if(entityDim == 2) periodicSurf(i,1:2) = [entityTag, entityTagMaster]
        read(MSH,*) numAffine
        read(MSH,*) numCorrespondingNodes
        do j = 1, numCorrespondingNodes
          read(MSH,*) nodeTag, nodeTagMaster
        end do
      end do
    else
      write(*,'(2X,A)') 'no periodic boundaries found'
    end if

    close(MSH)

    !---------------------------------------------------------------------------
    ! Process Gmsh mesh data

    ! perform necessary transformations
    !   - identification/numbering of vertices
    !   - identification of element points
    !   - identification of boundaries

    ! allocate high order points in 3D array for each element .................

    write(*,'(2X,A)') 'identifying collocation points'

    np         = meshOrder+1         ! number of nodes per direction
    numShells  = floor(meshOrder/2.) ! number of shells
    edgeCenter = meshOrder/2+1       ! position of center point on odd edge

    ! map element shells into 3D array points
    do i = 1, numHexElem
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
          deallocate(p_shell ,p_edge, p_face, points)
        end do
      end if

      ! map points in the element center

      associate(n => element(i)%nodes, points => element(i)%points)

        ! check if mesh order is odd or even
        select case(mod(meshOrder,2))
        case(0) ! even -> one central point
          points(edgeCenter,edgeCenter,edgeCenter) = n(size(n))

        case(1) ! odd  -> 8 corner points
          points(np/2+1,np/2,np/2)     = n(size(n) - 7)
          points(np/2+1,np/2+1,np/2)   = n(size(n) - 6)
          points(np/2+1,np/2+1,np/2+1) = n(size(n) - 5)
          points(np/2+1,np/2,np/2+1)   = n(size(n) - 4)
          points(np/2,np/2,np/2)       = n(size(n) - 3)
          points(np/2,np/2+1,np/2)     = n(size(n) - 2)
          points(np/2,np/2+1,np/2+1)   = n(size(n) - 1)
          points(np/2,np/2,np/2+1)     = n(size(n)    )
        end select

      end associate

    end do

    ! identify vertices ........................................................

    write(*,'(2X,A)') 'identifying vertices'
    allocate(vertex_mask(numNodes), source = 0)

    ! mark element vertex nodes
    do i = 1, numHexElem
      vertex_mask(element(i)%vertices) = 1
    end do

    ! only choose nodes, where the vertex mask is true
    k = 0
    do i = 1, numNodes
      if (vertex_mask(i) == 0) cycle
      k = k + 1
      vertex_mask(i) = k
    end do
    numVertices = k

    ! identify boundary elements and faces .....................................

    write(*,'(2X,A)') 'identifying boundary elements & faces'

    ! map node tags to faces 1 to 6

    do i = 1,numHexElem
      allocate(element(i)%faces(6,4))
      associate(v => element(i)%vertices)
        element(i)%faces(4,:) = [v(1),v(4),v(5),v(8)]
        element(i)%faces(2,:) = [v(2),v(3),v(6),v(7)]
        element(i)%faces(1,:) = [v(1),v(2),v(5),v(6)]
        element(i)%faces(3,:) = [v(3),v(4),v(7),v(8)]
        element(i)%faces(5,:) = v(1:4)
        element(i)%faces(6,:) = v(5:8)
      end associate
    end do


    ! identify elements corresponing to each boundary face .....................

    do i = 1,numBfaces
    do j = 1,numHexElem

      ! check if all 4 face corners correspond to vertices of element j
      do k = 1,4
        do l = 1,8
          if (bface(i)%nodes(k) == element(j)%vertices(l)) mask(k) = 1
        end do
      end do

      if(all(mask == 1)) then
        bface(i)%elem_id = element(j)%id

        ! check which side of element is corresponding
        do m = 1,6
          mask(:) = 0
          do k = 1,4
            where (element(j)%faces(m,:) == bface(i)%nodes(k)) mask = 1
          end do
          if(all(mask == 1)) bface(i)%elem_face = m
        end do
      end if

      mask(:) = 0

    end do
    end do

    ! map physical tags to boundary faces ......................................

    do i = 1, numBfaces
    do j = 1, numSurfaces
      if (boundaryList(j,1) == bface(i)%surfaceID) then
        bface(i)%physicalTag = boundaryList(j,3)
      end if
    end do
    end do

    ! get total number of elements on each boundary ............................

    do i = 1, numBoundaries
      boundary(i)%nElem = 0
      do j = 1, numSurfaces
        if(boundary(i)%physicalTag == boundaryList(j,3)) then
          boundary(i)%nElem = boundary(i)%nElem + boundaryList(j,4)
        end if
      end do
    end do

    ! process periodic boundaries if specified .................................

    if (io == 0) then

      write(*,'(2X,A)') 'processing periodic boundaries'

      do i = 1, numPeriodicLinks

        ! physical Tags of periodic surfaces
        do j = 1, numSurfaces
          if (periodicSurf(i,1) == boundaryList(j,2)) then
            periodicSurf(i,3) = boundaryList(j,3)
          end if
          if (periodicSurf(i,2) == boundaryList(j,2)) then
            periodicSurf(i,4) = boundaryList(j,3)
          end if
        end do

        ! boundary IDs of periodic surfaces
        do j = 1, numBoundaries
          if(periodicSurf(i,3) == boundary(j)%physicalTag) then
            periodicSurf(i,5) = j
          end if
          if(periodicSurf(i,4) == boundary(j)%physicalTag) then
            periodicSurf(i,6) = j
          end if
        end do
      end do

    end if

    ! give some mesh information as display output .............................

    write(*,*)
    write(*,'(2X,A)') 'mesh info'
    write(*,'(4X,A,T25,G0)') 'num vertices'       , numVertices
    write(*,'(4X,A,T25,G0)') 'num elements'       , numHexElem
    write(*,'(4X,A,T25,G0)') 'num boundaries'     , numBoundaries
    write(*,'(4X,A,T25,G0)') 'num boundary faces:', numBfaces
    write(*,'(4X,A,T25,G0)') 'geometry order:'    , meshOrder
    write(*,*)

    !---------------------------------------------------------------------------
    ! Create HiSPEET generic mesh

    ! create generic mesh vertices .............................................

    write(*,'(2X,A)') 'creating generic mesh vertices'

    allocate (mesh%vertex(numVertices))
    do i = 1, numNodes
      k = vertex_mask(i)
      if (k == 0) cycle ! no vertex
        mesh%vertex(k)%id = k
        mesh%vertex(k)%x  = node(i)%nodeCoord
    end do

    ! create generic mesh elements .............................................

    write(*,'(2X,A)') 'creating generic mesh elements'

    ! set rotational numbering
    mesh%numbering = ROTATIONAL_NUMBERING

    allocate(mesh%element(numHexElem))

    do i = 1, numHexElem

      ! ID, type, order and basis
      mesh%element(i)%id    = element(i)%id
      mesh%element(i)%typ   = HEXAHEDRAL_ELEMENT
      mesh%element(i)%order = meshOrder
      mesh%element(i)%basis = EQUIDISTANT_NODAL_BASIS

      ! vertices
      do k = 1,8
        mesh%element(i)%vertex(k) = vertex_mask(element(i)%vertices(k))
      end do

      ! collocation points
      allocate(mesh%element(i)%x(np*np*np, 3))
      associate(x => mesh%element(i)%x, point => element(i)%points)
        m = 0
        do j = np,1,-1
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

    ! create generic mesh boundaries ...........................................

    write(*,'(2X,A)') 'creating generic mesh boundaries'

    allocate(mesh%boundary(numBoundaries))

    l = 1
    do i = 1,numBoundaries

      mesh%boundary(i)%id   = i
      mesh%boundary(i)%name = boundary(i)%physicalName

      ! boundary faces
      allocate(mesh%boundary(i)%face(boundary(i)%nElem))
      k = 0
      do j = 1,numBfaces
        if (bface(j)%physicalTag == boundary(i)%physicalTag) then
          k = k + 1
          mesh%boundary(i)%face(k)%element_id   = bface(j)%elem_id
          mesh%boundary(i)%face(k)%element_face = bface(j)%elem_face
        end if
      end do

      ! periodic boundaries
      if (io == 0) then
        do j = 1, numPeriodicLinks
          if (periodicSurf(j,3) == boundary(i)%physicalTag) then
            mesh%boundary(i)%coupled = periodicSurf(j,6)
          end if
          if (periodicSurf(j,4) == boundary(i)%physicalTag) then
            mesh%boundary(i)%coupled = periodicSurf(j,5)
          end if
        end do
      end if

      ! map polarity for periodic boundaries
      if (mesh%boundary(i)%coupled /= 0 .and. mesh%boundary(i)%polarity == 0) then
        mesh%boundary(i)%polarity = l
        mesh%boundary(mesh%boundary(i)%coupled)%polarity = -l
        l = l + 1
      end if

      ! print boundary info
      write(*,'(I6,2X)',advance='NO') mesh%boundary(i)%id
      if (len_trim(trim(mesh%boundary(i)%name)) <= 57) then
        write(*,'(A57)',advance='NO') adjustl(mesh%boundary(i)%name)
      else
        write(*,'(A57)',advance='NO') adjustl(mesh%boundary(i)%name(1:56)//'')
      end if
      if (mesh%boundary(i)%coupled > 0) then
        write(*,'(A4,I4)',advance='NO') 'c:' ,mesh%boundary(i)%coupled
        write(*,'(A4,I3)',advance='NO') 'p:' ,mesh%boundary(i)%polarity
      end if
      write(*,*)

    end do

    ! enforce lexical numbering and generate consistent vertex IDs ..............

    call mesh % SwitchToLexicalNumbering()
    call mesh % GenerateConsistentVertexIDs()

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
