!> summary:  Import of generic 3d meshes
!> author:   Joerg Stiller
!> date:     2020/12/25
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(Mesh__3D) MP_ImportGenericMesh
  use Constants
  use Execution_Control
  use Gauss_Jacobi
  use Embedded_Interpolation_Operator__1D
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
        mesh % boundary(i) % map      = generic_mesh % boundary(i) % map
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
      call mesh % BuildSFC()

      allocate(mesh % map_child(0))
      allocate(mesh % map_parent(0))

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
    mesh % n_elem_active = mesh % n_elem
    mesh % n_elem_frozen = 0

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

    integer :: b, e, i, j

    do b = 1, size(generic_mesh%boundary)

      associate( generic_boundary => generic_mesh % boundary(b) &
               , mesh_boundary    => mesh % boundary(b)         )

        mesh_boundary % n_face = size(generic_boundary % face)
        allocate(mesh_boundary % face(mesh_boundary % n_face))

        do i = 1, mesh_boundary % n_face

          e = generic_boundary % face(i) % element_id   ! corresp element ID
          j = generic_boundary % face(i) % element_face ! corresp element face

          mesh_boundary % face(i) % element_id   = e
          mesh_boundary % face(i) % element_face = j

        end do

      end associate
    end do

  end subroutine ImportBoundaryFaces

  !-----------------------------------------------------------------------------
  !> Import of element domains

  subroutine ImportElementDomains(mesh, generic_mesh)
    class(Mesh_3D),        intent(inout) :: mesh         !< mesh partition
    class(GenericMesh_3D), intent(in)    :: generic_mesh !< generic mesh

    type(EmbeddedInterpolationOperator_1D), allocatable :: iop(:)
    real(RNP), allocatable :: xo(:), xi(:)
    integer :: e, i, po, po_max, po_min
    logical :: has_equidistant_nodes

    ! preliminaries ............................................................

    associate(p_geom => mesh % p_geom)

      has_equidistant_nodes = .false.
      po_min = huge(1)
      po_max = -1
      p_geom = -1

      do e = 1, size(generic_mesh % element)
        po = generic_mesh % element(e) % order
        if (generic_mesh % element(e) % basis == EQUIDISTANT_NODAL_BASIS) then
          has_equidistant_nodes = .true.
          po_min = min(po, po_min)
          po_max = max(po, po_max)
        else if (generic_mesh % element(e) % basis /= GAUSS_LOBATTO_BASIS) then
          call Error('ImportElementDomains', 'basis not supported')
        end if
        p_geom = max(p_geom, po)
      end do

      ! build interpolation operators
      if (has_equidistant_nodes) then
        allocate(iop(po_min:po_max))
        do po = po_min, po_max
          allocate(xo(0:po), xi(0:po))
          ! original points (equidistant)
          do i = 0, po
            xo(i) = TWO * i/po - ONE
          end do
          ! interpolation points (Lobatto)
          xi = LobattoPoints(po)
          iop(po) = EmbeddedInterpolationOperator_1D('N', xo, xi)
          deallocate(xo, xi)
        end do
      end if

    end associate

    ! element degree and points ................................................

    do e = 1, mesh % n_elem
      associate( geometry => mesh % element(e) % geometry  &
               , po  => generic_mesh % element(e) % order  &
               , x_g => generic_mesh % element(e) % x      )

        allocate(geometry % x_e(0:po, 0:po, 0:po, 1:3))

        geometry % po = po

        if (generic_mesh % element(e) % basis == EQUIDISTANT_NODAL_BASIS) then
          call TransformToLobattoBasis(po, iop(po) % A, x_g, geometry % x_e)
        else
          geometry % x_e(0:po,0:po,0:po,1:3) = reshape(x_g, [po+1,po+1,po+1,3])
        end if

      end associate
    end do

  end subroutine ImportElementDomains

  !-----------------------------------------------------------------------------
  !> Transformation to Gauss-Lobatto nodes

  pure subroutine TransformToLobattoBasis(po, A, x_g, x_e)
    integer,   intent(in)  :: po
    real(RNP), intent(in)  :: A(0:po,0:po)
    real(RNP), intent(in)  :: x_g(0:po,0:po,0:po,3)
    real(RNP), intent(out) :: x_e(0:po,0:po,0:po,3)

    real(RNP) :: z2(0:po,0:po,0:po), z3(0:po,0:po,0:po)
    real(RNP) :: tmp
    integer   :: d, i, j, k, p

    do d = 1, 3

      ! z3 = AxIxI x_g(:,:,:,d) ................................................

      do k = 0, po
      do j = 0, po
      do i = 0, po
        tmp = 0
        do p = 0, po
          tmp = tmp + A(k,p) * x_g(i,j,p,d)
        end do
        z3(i,j,k) = tmp
      end do
      end do
      end do

      ! z2 = IxAxI z3 ..........................................................

      do k = 0, po
      do j = 0, po
      do i = 0, po
        tmp = 0
        do p = 0, po
          tmp = tmp + A(j,p) * z3(i,p,k)
        end do
        z2(i,j,k) = tmp
      end do
      end do
      end do

      ! x_e(:,:,:,d) = IxIxA z2 ................................................

      do k = 0, po
      do j = 0, po
      do i = 0, po
        tmp = 0
        do p = 0, po
          tmp = tmp + A(i,p) * z2(p,j,k)
        end do
        x_e(i,j,k,d) = tmp
      end do
      end do
      end do

    end do

  end subroutine TransformToLobattoBasis

  !=============================================================================

end submodule MP_ImportGenericMesh
