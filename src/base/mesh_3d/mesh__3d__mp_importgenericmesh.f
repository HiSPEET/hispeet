!> summary:  Import of generic 3d meshes
!> author:   Joerg Stiller
!> date:     2020/12/25
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(Mesh__3D) MP_ImportGenericMesh
  use Constants
  use Execution_Control
  use Generic_Mesh__3D
  implicit none

  interface

    !---------------------------------------------------------------------------
    !> Identification of neighbor elements

    module subroutine IdentifyElementNeighbors(mesh)
      class(Mesh_3D), intent(inout) :: mesh !< mesh partition
    end subroutine IdentifyElementNeighbors

  end interface

contains

  !-----------------------------------------------------------------------------
  !> Import a generic 3d mesh

  module subroutine ImportGenericMesh(mesh, generic_mesh, comm)
    class(Mesh_3D),        intent(out) :: mesh         !< mesh partition
    class(GenericMesh_3D), intent(in)  :: generic_mesh !< generic mesh
    type(MPI_Comm),        intent(in)  :: comm         !< MPI communicator

    type(MeshAttributes_3D) :: attrib
    integer :: i, rank

    call MPI_Comm_rank(comm, rank)

    if (rank == 0) then

      ! check prerequisites ....................................................

      if (generic_mesh % numbering /= LEXICAL_NUMBERING) then
        call Error('ImportGenericMesh', 'Lexical numbering is required')
      else if (any(generic_mesh % element%typ /= HEXAHEDRAL_ELEMENT)) then
        call Error('ImportGenericMesh', 'Only hexahedral elements supported')
      else
        write(*,'(/,A)') 'importing generic mesh'
      end if

      ! basic initialization ...................................................

      mesh % n_bound     =  size(generic_mesh % boundary)
      mesh % n_parts     =  1

      mesh % structured  = .false.
      mesh % regular     = .false.

      mesh % comm_world  =  comm
      mesh % proc        =  0
      mesh % part        =  0

      allocate(mesh % boundary( mesh%n_bound ))
      do i = 1, mesh % n_bound
        mesh % boundary(i) % name     = generic_mesh % boundary(i) % name
        mesh % boundary(i) % id       = i
        mesh % boundary(i) % coupled  = generic_mesh % boundary(i) % coupled
        mesh % boundary(i) % polarity = generic_mesh % boundary(i) % polarity
      end do

      ! mesh components ........................................................

      call ImportElements(mesh, generic_mesh)

      call mesh % IdentifyEdges()
      call mesh % BuildFaces()

      call ImportBoundaryFaces(mesh, generic_mesh)
      call ImportElementDomains(mesh, generic_mesh)
      call IdentifyElementNeighbors(mesh)

      call mesh % BuildLinks()
      call mesh % BuildGhosts()
      call mesh % BuildCuboids()
      call mesh % IdentifyRanks()

    end if

    ! empty partitions and intracommunicator ...................................

    if (rank == 0) then
      attrib = MeshAttributes_3D(mesh)
    end if

    call attrib % Bcast(0, comm)

    if (rank > 0) then
      call mesh % Init_Mesh_3D(attrib, comm)
    end if

    ! intracommunicator between active partitions
    call mesh % BuildCommunicator()

    ! print info ...............................................................

    if (rank == 0) then
      write(*,'(T3,A,T30,9(G0,1X))') 'number of elements:   ', mesh % n_elem
      write(*,'(T3,A,T30,9(G0,1X))') 'number of faces:      ', mesh % n_face
      write(*,'(T3,A,T30,9(G0,1X))') 'number of edges:      ', mesh % n_edge
      write(*,'(T3,A,T30,9(G0,1X))') 'number of vertices:   ', mesh % n_vert
      write(*,'(T3,A,T30,9(G0,1X))') 'number of boundaries: ', mesh % n_bound
    end if

  end subroutine ImportGenericMesh

  !-----------------------------------------------------------------------------
  !> Import of mesh elements
  !>
  !> Allocates the elements, sets local IDs, and adopts the vertex IDs from the
  !> generic mesh.

  subroutine ImportElements(mesh, generic_mesh)
    class(Mesh_3D),        intent(inout) :: mesh         !< mesh partition
    class(GenericMesh_3D), intent(in)    :: generic_mesh !< generic mesh

    integer :: b, i, j, k, n

    mesh % n_elem = size(generic_mesh % element)
    allocate(mesh % element( mesh % n_elem ))

    ! element and vertex IDs ...................................................

    n = 0
    do i = 1, mesh % n_elem
      mesh % element(i) % id = i
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
    class(Mesh_3D),        intent(inout) :: mesh         !< mesh partition
    class(GenericMesh_3D), intent(in)    :: generic_mesh !< generic mesh

    integer :: b, e, f, i, j, s

    do b = 1, size(generic_mesh%boundary)

      associate( generic_boundary => generic_mesh % boundary(b) &
               , mesh_boundary    => mesh % boundary(b)         )

        mesh_boundary % n_face = size(generic_boundary % face)
        allocate(mesh_boundary % face(mesh_boundary % n_face))

        do i = 1, mesh_boundary % n_face

          e = generic_boundary % face(i) % element_id   ! corresp element ID
          j = generic_boundary % face(i) % element_face ! corresp element face
          f = mesh % element(e) % face(j) % id          ! corresp mesh face ID

          ! side of the mesh face on which the boundary is located
          if (mesh % face(f) % element(1) % id == e) then
            s = 2
          else
            s = 1
          end if

          mesh_boundary % face(i) % mesh_face % id   = f
          mesh_boundary % face(i) % mesh_face % side = s

          mesh_boundary % face(i) % mesh_element % id   = e
          mesh_boundary % face(i) % mesh_element % face = j

        end do

      end associate
    end do

  end subroutine ImportBoundaryFaces

  !-----------------------------------------------------------------------------
  !> Import of element domains

  subroutine ImportElementDomains(mesh, generic_mesh)
    class(Mesh_3D),        intent(inout) :: mesh         !< mesh partition
    class(GenericMesh_3D), intent(in)    :: generic_mesh !< generic mesh

    integer :: e

    ! preliminaries ............................................................

    mesh % p_geom = -1
    do e = 1, size(generic_mesh % element)
      if (generic_mesh % element(e) % basis /= GAUSS_LOBATTO_BASIS) then
        call Error('ImportElementDomains', 'basis not supported')
      end if
      mesh % p_geom = max(mesh % p_geom, generic_mesh % element(e) % order)
    end do

    ! element degree and points ................................................

    do e = 1, mesh % n_elem
      associate( geometry => mesh % element(e) % geometry &
               , po => generic_mesh % element(e) % order  &
               , xg => generic_mesh % element(e) % x      )

        allocate(geometry % x_e(0:po, 0:po, 0:po, 1:3))

        geometry % po = po
        geometry % x_e(0:po,0:po,0:po,1:3) = reshape(xg, [po+1,po+1,po+1,3])

      end associate
    end do

  end subroutine ImportElementDomains

  !=============================================================================

end submodule MP_ImportGenericMesh
