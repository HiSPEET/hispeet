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
  use Execution_Control   ! to throw warnings or errors
  use Generic_Mesh__3D

  implicit none
  private

  public :: Import_GMSH_3D

  interface

    module subroutine getCESpoints(v,e,s,shellOrder,p)
      integer, intent(in)  :: shellOrder                !< shell order
      integer, intent(in)  :: v(:), e(:), s(:)          !< points of vertices, edges, surfaces
      integer, allocatable, intent(out) :: p(:,:,:)     !< 3D point array to be returned
    end subroutine getCESpoints

  end interface


  ! -----------------------------------------------------------------------
  !> Typ description
  type :: MshNode
    integer :: nodeTag     = 0            !< node tag
    real(RNP) :: nodeCoord(3) = 0         !< node coordinates
  end type MshNode

  type :: MshElement
    integer :: id          = 0            !< element id
    integer :: elementTag  = 0            !< element tag
    integer :: vertices(8) = 0            !< hexadedral element corner vertices
    integer, allocatable :: faces(:,:)    !< array to map vertex ids to face id
    integer, allocatable :: nodes(:)      !< 1D node vector
    integer, allocatable :: points(:,:,:) !< 3D point array
  end type MshElement

  type :: MshBoundary
    character(len = 127) :: physicalName  !< boundary name
    integer :: physicalTag = 0            !< corresponding physical tags
    integer :: nElem       = 0            !< number of corresponding elements
  end type MshBoundary

  type :: MshBoundaryFace
    integer :: id          = 0            !< numeration id
    integer :: surfaceID   = 0            !< numeration id of corresponding surface
    integer :: surfaceTag  = 0            !< tag of corresponding surface
    integer :: physicalTag = 0            !< corresponding physical tags
    integer :: nodes(4)    = 0            !< 4 corner nodes (= hexa element vertices)
    integer :: elem_id     = 0            !< element tag
    integer :: elem_face   = 0            !< element face
  end type MshBoundaryFace


contains


  !-----------------------------------------------------------------------------
  !> Reads the mesh from file and converts it to a GenericMesh_3D object

  subroutine Import_GMSH_3D(mshfile, generic_mesh)
    character(len=*),      intent(in)  :: mshfile       !< GMSH file (*.msh)
    class(GenericMesh_3D), intent(out) :: generic_mesh  !< output generic mesh

    ! MSH types declared above

    type(MshNode),         allocatable  :: node(:)
    type(MshElement),      allocatable  :: element(:)
    type(MshBoundary),     allocatable  :: boundary(:)
    type(MshBoundaryFace), allocatable  :: bface(:)


    ! auxiliary variables ......................................................

    integer :: i, io, j, k, l, m, MSH, n, param, pT, sT
    integer :: numPoints, numCurves, numSurfaces, numVolumes
    integer :: numnodes, numElements, numVertices
    integer :: numBoundaries, numBfaces, numHexElem
    integer :: numEntityBlocks, entityDim, entityTag
    integer :: numPhysicalNames, numPhysicalTags, physicalTag
    integer :: numNodesInBlock, numElementsInBlock
    integer :: numPeriodicLinks, entityTagMaster, nodeTagMaster
    integer :: numCorrespondingNodes, numAffine, nodeTag, maxNodeTag
    integer :: elementType, meshOrder, numShells_total, shellOrder
    integer :: edgecenter, np_dir, np_total
    integer :: iv(8), mask(4), p_vert(8)
    real    :: minX, minY, minZ, maxX, maxY, maxZ
    character(len = 100) :: mshblock

    integer, allocatable :: node_nodeTag(:)     !< mask to identify the node tags for efficient point coordinate mapping
    integer, allocatable :: vertexmask(:)       !< mask to identify the element corner nodes that correspond to a vertex

    integer, allocatable :: boundarymap(:,:)    !< map surface ids & tags to boundary physical tag and number of elements
    integer, allocatable :: periodicSurf(:,:)   !< array to store physical tags [...] of periodical surfaces

    integer, allocatable :: shellarray(:,:,:)   !< 3D array containing the node tags of a hexaedral cube shell
    integer, allocatable :: p_shell(:)          !< vector to extract current shell node tags from element%nodes vector
    integer, allocatable :: p_edge(:)           !< vector to extract edge node tags from p_shell vector
    integer, allocatable :: p_surface(:)        !< vector to extract surface node tags from p_shell vector


integer :: most_frequent, max_count, counter

    !---------------------------------------------------------------------------
    ! Read Gmsh mesh file

    ! open GMSH file
    write(*,'(/,A,/)') 'Opening GMSH file '//trim(mshfile)//'.msh'

    open(newunit = MSH, file = trim(mshfile)//'.msh')

    ! skip first three lines
    read(MSH,'(2/)')


    ! check if boundaries have been specified properly
    read(MSH,*) mshblock
    if(mshblock /= '$PhysicalNames') then
      print *, 'Error: No boundaries specified. Program will abort'
      stop
    end if


    ! read $PhysicalNames block ................................................

    read(MSH,*) numPhysicalNames
    numBoundaries = numPhysicalNames - 1    ! ignore volume entity
    allocate(boundary(numBoundaries))

    do i = 1, numBoundaries
      read(MSH,*) entityDim, boundary(i)%physicalTag, boundary(i)%physicalName
    end do

    read(MSH,'(2/)')



    ! read $Entities block .....................................................

    read(MSH,*) numPoints, numCurves, numSurfaces, numVolumes
    numEntityBlocks = numPoints+numCurves+numSurfaces+numVolumes
    allocate(boundarymap(numSurfaces,4), source = 0)

    ! ignore points & curves
    do i = 1, numPoints+numCurves
      read(MSH,*)
    end do

    ! get tags of mesh boundary surfaces
    j = 0
    do i = 1, numSurfaces
      read(MSH,*) sT, minX, minY, minZ, &
           maxX, maxY, maxZ, numPhysicalTags, pT
      if (numPhysicalTags == 0) cycle  ! surface i is not a boundary
      j = j+1
      boundarymap(j,1:3) = [j, sT, pT]
    end do

    ! ignore volumes
    do i = 1, numVolumes
      read(MSH,*)
    end do



    ! read $Nodes block.........................................................

    print *, 'reading node data ...'
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
    allocate(node_nodeTag(maxNodeTag), source = -1)
    do k = 1, numNodes
      node_nodeTag(node(k)%nodeTag) = k
    end do


    ! read $Elements block .....................................................

    print *, 'reading boundary and element data ...'

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

          ! 2D element faces give information about mesh boundaries
          case(2)
            boundarymap(i,4) = numElementsInBlock !number of elements on each boundary surface
            numBfaces = numBfaces + 1
            bface(numBfaces)%id = numBfaces
            bface(numBfaces)%surfaceTag = entityTag
            bface(numBfaces)%surfaceID = i
            read(MSH,*) k , bface(numBfaces)%nodes

          ! 3D hexaedral element
          case(3)
            numHexElem = numHexElem + 1
            element(numHexElem)%id = numHexElem
            meshOrder = getorder(elementType)
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

    if(io == 0) then
      print *, 'reading periodic boundary data ...'

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
      print *, 'no periodic boundaries have been imposed'
    end if


    close(MSH)


    !---------------------------------------------------------------------------
    !> Process Gmsh mesh data

    ! perform necessary transformations ...
    !   - identification/numbering of vertices
    !   - identification of element points
    !   - identification of boundaries

    ! allocate high order points in 3D matrix for each element .................

    print *, 'reallocating high order points ...'

    np_dir           = meshOrder+1            ! number of points per direction (edge with vertices included)
    np_total         = np_dir**3              ! total number of points per element
    numShells_total  = floor(meshOrder/2.)    ! total number of shells
    edgecenter       = meshOrder/2+1          ! center point of odd edge


    ! map element nodes into 3D array points
    do i = 1, numHexElem
      allocate(element(i)%points(np_dir,np_dir,np_dir), source = 1)

      element(i)%vertices = element(i)%nodes(1:8)

      if(meshOrder > 1) then
        ! subdivide nodes vector into shells p_shell
        shellOrder = meshOrder
        k = 1
        l = 0
        do j = 1, numShells_total
          allocate(shellarray(np_dir-l,np_dir-l,np_dir-l), source = 1)
          allocate(p_shell(np_shell(shellOrder)))

          p_shell = element(i)%nodes(k:k+np_shell(shellOrder)-1)

          ! subdivide p_shell into corners, edges and surfaces
          allocate(p_edge((shellOrder-1)*12), source = 1)
          allocate(p_surface(((shellOrder-1)**2)*6), source = 1)

          p_vert    = p_shell(1:8)
          p_edge    = p_shell(9:8+size(p_edge))
          p_surface = p_shell(size(p_edge)+9:size(p_shell))

          ! call subroutine to map the shell points on 3D shellmatrix

          call getCESpoints(p_vert,p_edge,p_surface,shellOrder,shellarray)
          element(i)%points(j:np_dir-(j-1),j:np_dir-(j-1),j:np_dir-(j-1)) = shellarray

          ! prepare variables for next shell corresponding to hexaedral with order reduced by 2
          k = k + np_shell(shellOrder)
          l = l + 2
          shellOrder = shellOrder - 2
          deallocate(p_shell,p_edge,p_surface,shellarray)
        end do
      end if

      ! map residual incomplete element points in center

      associate(n => element(i)%nodes, p => element(i)%points)

      ! check if mesh order is odd or even
      select case(mod(meshOrder,2))
      case(0)     ! even -> one central point
        p(edgecenter,edgecenter,edgecenter) = n(size(n))

      case(1)     ! odd  -> 8 corner points
        p(np_dir/2+1,np_dir/2,np_dir/2)        = n(size(n) - 7)
        p(np_dir/2+1,np_dir/2+1,np_dir/2)      = n(size(n) - 6)
        p(np_dir/2+1,np_dir/2+1,np_dir/2+1)    = n(size(n) - 5)
        p(np_dir/2+1,np_dir/2,np_dir/2+1)      = n(size(n) - 4)
        p(np_dir/2,np_dir/2,np_dir/2)          = n(size(n) - 3)
        p(np_dir/2,np_dir/2+1,np_dir/2)        = n(size(n) - 2)
        p(np_dir/2,np_dir/2+1,np_dir/2+1)      = n(size(n) - 1)
        p(np_dir/2,np_dir/2,np_dir/2+1)        = n(size(n))
      end select
      end associate

    end do


    ! identify vertices ........................................................

    print *, 'identifying vertices ...'
    allocate(vertexmask(numNodes), source = 0)

    ! mark element vertex nodes
    do i = 1, numHexElem
      vertexmask(element(i)%vertices) = 1
    end do

    ! identify mesh vertices
    k = 0
    do i = 1, numNodes
      if (vertexmask(i) == 0) cycle
      k = k + 1
      vertexmask(i) = k
    end do
    numVertices = k



    ! identify boundary elements and faces .....................................

    print *, 'identifying boundary elements & faces ...'

    ! map node tags to faces 1 to 6

    do i = 1,numHexElem
      allocate(element(i)%faces(6,4))
      associate(v => element(i)%vertices)
      element(i)%faces(4,:) = [v(1),v(4),v(5),v(8)] !4
      element(i)%faces(2,:) = [v(2),v(3),v(6),v(7)] !2
      element(i)%faces(1,:) = [v(1),v(2),v(5),v(6)] !1
      element(i)%faces(3,:) = [v(3),v(4),v(7),v(8)] !3
      element(i)%faces(5,:) = v(1:4) !5
      element(i)%faces(6,:) = v(5:8) !6
      end associate
!      print *, 'elemid:', element(i)%id, element(i)%vertices
!        do j = 1,6
!        print *, j, 'facenodes', element(i)%faces(j,:)
!        end do
    end do


    ! identify elements corresponing to each boundary face .....................

    do i = 1,numBfaces
!print *, bface(i)%nodes
      do j = 1,numHexElem
        ! check if all 4 face corners correspond to vertices of element j
        do k = 1,4
          do l = 1,8
            if (bface(i)%nodes(k) == element(j)%vertices(l)) mask(k) = 1
          end do
        end do
        if(all(mask == 1)) then     ! corresponding element found
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
        if (boundarymap(j,1) == bface(i)%surfaceID) then
          bface(i)%physicalTag = boundarymap(j,3)
        end if
      end do
    end do


    ! get total number of elements on each boundary ............................

    do i = 1, numBoundaries
      boundary(i)%nElem = 0
      do j = 1, numSurfaces
        if(boundary(i)%physicalTag == boundarymap(j,3)) then
          boundary(i)%nElem = boundary(i)%nElem + boundarymap(j,4)
        end if
      end do
    end do




    ! process periodic boundaries if specified .................................

    if(io == 0) then

      print *, 'processing periodic boundaries ...'

      do i = 1, numPeriodicLinks
        ! get physical Tags of periodic surfaces
        do j = 1, numSurfaces
          if (periodicSurf(i,1) == boundarymap(j,2)) periodicSurf(i,3) = boundarymap(j,3)
          if (periodicSurf(i,2) == boundarymap(j,2)) periodicSurf(i,4) = boundarymap(j,3)
        end do
        ! get boundary ids of periodic surfaces
        do j = 1, numBoundaries
          if(periodicSurf(i,3) == boundary(j)%physicalTag) periodicSurf(i,5) = j
          if(periodicSurf(i,4) == boundary(j)%physicalTag) periodicSurf(i,6) = j
        end do
      end do

    end if




    !---------------------------------------------------------------------------
    !> Create HiSPEET generic mesh

    ! create generic mesh vertices .............................................

    print *, 'generating generic mesh vertices ...'

    ! map only the corner nodes of the node list which correspond to a vertex
    allocate (generic_mesh%vertex(numVertices))
    do i = 1, numNodes
      k = vertexmask(i)
      if (k == 0) cycle
        generic_mesh%vertex(k)%id = k
        generic_mesh%vertex(k)%x  = node(i)%nodeCoord
    end do


    ! create generic mesh elements .............................................

    print *, 'generating generic mesh elements ...'

    ! set rotational numbering
    generic_mesh%numbering = ROTATIONAL_NUMBERING

    allocate(generic_mesh%element(numHexElem))

    do i = 1, numHexElem
      ! ID, type, order and basis
      generic_mesh%element(i)%id    = element(i)%id
      generic_mesh%element(i)%typ   = HEXAHEDRAL_ELEMENT
      generic_mesh%element(i)%order = meshOrder
      generic_mesh%element(i)%basis = EQUIDISTANT_NODAL_BASIS

      ! vertices
      do k = 1,8
      generic_mesh%element(i)%vertex(k)  = vertexmask(element(i)%vertices(k))
      end do

      ! high order points
      allocate(generic_mesh%element(i)%x(np_total,3))

      ! map 3D point array to generic mesh element following a lexical numbering
      associate(x => generic_mesh%element(i)%x, point => element(i)%points)
      m = 0
      do j = np_dir,1,-1
        do l = 1, np_dir
          do k = 1, np_dir
            m = m+1
            n = node_nodeTag(point(j,k,l))

            x(m,1) = node(n)%nodeCoord(1)
            x(m,2) = node(n)%nodeCoord(2)
            x(m,3) = node(n)%nodeCoord(3)
          end do
        end do
      end do
      end associate

    end do



    ! create generic mesh boundaries ...........................................

    print *, 'generating generic mesh boundaries ...'

    allocate(generic_mesh%boundary(numBoundaries))

    do i = 1,numBoundaries

      ! identifier
      generic_mesh%boundary(i)%id = i

      ! name
      generic_mesh%boundary(i)%name = boundary(i)%physicalName

      ! create boundary faces
      allocate(generic_mesh%boundary(i)%face(boundary(i)%nElem))

      k = 0
      do j = 1,numBfaces
        if (bface(j)%physicalTag == boundary(i)%physicalTag) then
          k = k+1
          generic_mesh%boundary(i)%face(k)%element_id    = bface(j)%elem_id
          generic_mesh%boundary(i)%face(k)%element_face  = bface(j)%elem_face
        end if
      end do


      ! coupled boundaries
      if (io == 0) then
        do j = 1, numPeriodicLinks
          if(periodicSurf(j,3) == boundary(i)%physicalTag) generic_mesh%boundary(i)%coupled = periodicSurf(j,6)
          if(periodicSurf(j,4) == boundary(i)%physicalTag) generic_mesh%boundary(i)%coupled = periodicSurf(j,5)
        end do
      end if

    end do

    ! enforce lexical numbering ................................................

    call generic_mesh % SwitchToLexicalNumbering()

    ! give some mesh information as display output .............................

    write(*,*)
    write(*,'(2X,A,I2)') 'Mesh order:', meshOrder
    write(*,'(2X,A)') 'Mesh contains'
    write(*,'(2X,I8,A)') numVertices,   ' vertices'
    write(*,'(2X,I8,A)') numHexElem,    ' hexaedral elements'
    write(*,'(2X,I8,A)') numBoundaries, ' boundaries'
    write(*,'(2X,I8,A)') numBfaces,     ' faces on mesh boundaries'



  end subroutine Import_GMSH_3D


  !> Get mesh order: Returns mesh order according to element type numbers defined by GMSH

  function getorder(elementType) result(meshOrder)
    integer, intent (in) :: elementType         !< element type number defined by GMSH
    integer              :: meshOrder           !< order of hexaedral mesh

    select case(elementType)
    case(5)
      meshOrder = 1
    case(12)
      meshOrder = 2
    case(92:98)
      meshOrder = elementType - 89
    case default
      print *, 'Element type is not supported. Program will abort ...'
      stop
    end select
  end function



  !> Function np_shell returns number of points per shell for a given shell order

  function np_shell(shellOrder) result(npoints)
    integer, intent (in) :: shellOrder          !< shell order
    integer              :: npoints             !< number of points per shell

    npoints = (shellOrder+1)**3 - (shellOrder-1)**3
  end function

  !=============================================================================

end module Import_GMSH__3D
