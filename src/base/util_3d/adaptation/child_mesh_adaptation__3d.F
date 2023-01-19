module Child_Mesh_Adaptation__3D
  use Kind_Parameters
  use Constants
  use Execution_Control
  use XMPI
  use Mesh__3D
  use Mesh_Element__3D
  use Mesh_Element_Indexing__3D
  use Element_Transfer_Buffer__3D
  use Element_Distribution_Map__3D
  use Element_Distribution_Data__3D
  use Child_Distribution_Map__3D
  use Partitioner_Interface__3D

  implicit none
  private

  public :: PartitioningOptions_3D   ! from Partitioner_Interface__3D
  public :: ChildMeshAdaptation_3D

contains

  !=============================================================================
  ! Adaptation routines

  !-----------------------------------------------------------------------------
  !> Generation of an adapted child mesh

  subroutine ChildMeshAdaptation_3D(opt, parent, child)
    class(PartitioningOptions_3D), intent(in) :: opt
    class(Mesh_3D), intent(inout) :: parent
    class(Mesh_3D), intent(out)   :: child

    ! internal variables .......................................................

    type(MeshAttributes_3D)       :: attrib
    type(ChildDistributionMap_3D) :: child_map
    integer :: i

    !---------------------------------------------------------------------------
    !$omp master

    ! prepare new partitions ...................................................

    attrib = MeshAttributes_3D(parent)
    attrib % structured = .false.
    attrib % is_root = .false.
    attrib % dx = HALF * attrib % dx

    call child % Init_Mesh_3D(attrib, parent % comm_world)

    child % n_parts = opt % n_parts

    ! use identity for mapping parts to processes, so far
    if (child % proc < child % n_parts) then
      child % part = child % proc
    end if

    call child % BuildCommunicator()

    ! build child element map ..................................................

    if (parent % part >= 0) then
      call BuildChildMap(opt, parent, child_map)
    end if

    !---------------------------------------------------------------------------
    !$omp end master

  end subroutine ChildMeshAdaptation_3D

  !-----------------------------------------------------------------------------
  !> Computes the target family partitions and first child element IDs

  subroutine BuildChildMap(opt, parent, map)
    class(PartitioningOptions_3D), intent(in)  :: opt
    class(Mesh_3D),                intent(in)  :: parent
    type(ChildDistributionMap_3D), intent(out) :: map

    integer, allocatable, target :: tp_child(:)
    integer, contiguous, pointer :: tp_child_val(:,:,:,:)
    type(ElementTransferBuffer_3D), allocatable, asynchronous :: tp_child_buf

    allocate(tp_child(parent%n_elem + parent%n_ghost), source = -1)

    ! graph-based partitioning
    call ParMETIS_Partitioner(opt, parent, tp_child)

    ! transfer target partition IDs to ghosts
    if (parent % n_ghost > 0) then
      tp_child_val(1:1,1:1,1:1,1:parent%n_elem+parent%n_ghost) => tp_child
      tp_child_buf = ElementTransferBuffer_3D(parent, tp_elem_val)
      call tp_child_buf % Transfer(parent, tp_child_val, tag = 1000)
      call tp_child_buf % Merge(tp_child_val)
    end if

    map = ChildDistributionMap_3D(parent, opt%n_parts, tp_child)

  end subroutine BuildChildMap

  !=============================================================================
  ! Generation of child elements

  !-----------------------------------------------------------------------------
  !> Build child distribution data for sending from parent to child partitions

  subroutine BuildDistributionData(parent, child, map, dist_data)

    class(Mesh_3D), intent(in) :: parent
    class(Mesh_3D), intent(in) :: child
    class(ChildDistributionMap_3D),    intent(in)  :: map
    class(ElementDistributionData_3D), intent(out) :: dist_data(0:)

    !---------------------------------------------------------------------------
    ! internal variables

    integer, allocatable :: id_child_face(:,:,:) ! IDs of neighbor children
    integer, allocatable :: id_child_edge(:,:,:) ! at parent element faces,
    integer, allocatable :: id_child_vert(:,:)   ! edges and vertices

    integer, allocatable :: tp_child_face(:)     ! TPs of neighbor children
    integer, allocatable :: tp_child_edge(:,:)   ! at parent element faces,
    integer, allocatable :: tp_child_vert(:,:)   ! edges and vertices

    integer :: ce_part(0:map%n_parts-1)  ! element counters
    integer :: cn_part(0:map%n_parts-1)  ! neighbor counters
    integer :: cp_part(0:map%n_parts-1)  ! element point counters

    integer :: nc, nn, np, tp
    integer :: nn_max, nn_face, nn_edge, nn_vert
    integer :: e, i, j, k, l, m
    integer :: ce

    if (parent % part < 0) return

    associate(id_child => map % id_child, tp_child => map % tp_child)

      !-------------------------------------------------------------------------
      ! MPI and dimensions

      ! basic initialization
      do tp = 0, map%n_parts-1
        dist_data(tp) % proc       = child % proc_part(tp)
        dist_data(tp) % comm       = child % comm_world
        dist_data(tp) % n_element  = 0
        dist_data(tp) % n_neighbor = 0
        dist_data(tp) % n_point    = 0
      end do

      ! counters
      ce_part = 0  ! count of child elements per partition
      cn_part = 0  ! upper bound for count of child neighbors
      cp_part = 0  ! upper bound for count of child geometry points

      ! max dimensions
      nn_max  = 0  ! max number of neighbors per child element
      nn_face = 0  ! max number of neighbors per child element face
      nn_edge = 0  ! max number of neighbors per child element edge
      nn_vert = 0  ! max number of neighbors per child element vertex

      ! counts of elements, neighbors and geometry points per partition
      do e = 1, parent % n_elem
        associate(element => parent % element(e))
          tp = tp_child(e)
          if (tp >= 0) then

            ! number of child elements
            select case(element % adaptation % mark)
            case(1:6)
              nc = 4
            case(7:18)
              nc = 2
            case(19:26)
              nc = 1
            case default
              nc = 8
            end select

            ! max number of neighbors at faces, edges and vertices
            i = max( element % face   % n_neighbor, 1)
            j = max( element % edge   % n_neighbor, 1)
            k = max( element % vertex % n_neighbor, 1)

            ! max number of neighbors per child
            nn = 6 * i + 12 * j + 8 * k

            ! number of element points per child
            np = (element % geometry % po + 1) ** 3

            ! update counters
            ce_part(tp) = ce_part(tp) + nc
            cn_part(tp) = cn_part(tp) + nc * nn
            cp_part(tp) = cp_part(tp) + nc * np

            ! update dimensions
            nn_max  = max(nn, nn_max)
            nn_face = max(i, nn_face)
            nn_edge = max(j, nn_edge)
            nn_vert = max(k, nn_vert)

          end if
        end associate
      end do

      if (nn_face > 1) then
        call Error( 'BuildDistributionData'                             &
                  , 'adaptation of nonconforming mesh is not supported' &
                  , 'Child_Mesh_Adaptation__3D'                         )
      end if

      !-------------------------------------------------------------------------
      ! build components

      ! setup ..................................................................

      ! allocate components
      do tp = 0, map%n_parts-1
        allocate( dist_data(tp) % element              ( ce_part(tp)    ) )
        allocate( dist_data(tp) % neighbor_id          ( cn_part(tp)    ) )
        allocate( dist_data(tp) % neighbor_part        ( cn_part(tp)    ) )
        allocate( dist_data(tp) % neighbor_component   ( cn_part(tp)    ) )
        allocate( dist_data(tp) % neighbor_orientation ( cn_part(tp)    ) )
        allocate( dist_data(tp) % x_e                  ( cp_part(tp), 3 ) )
      end do

      ! auxiliary arrays
      allocate( id_child_face (2,2,6)        , tp_child_face (6)          )
      allocate( id_child_edge (2,nn_edge,12) , tp_child_edge (nn_edge,12) )
      allocate( id_child_vert (nn_vert,8)    , tp_child_vert (nn_vert,8)  )

      ! reset counters
      ce_part = 0
      cn_part = 0
      cp_part = 0

      ! generate child element data
      do e = 1, parent % n_elem

        tp = tp_child(e)
        if (tp < 0) cycle

        associate(element => parent % element(e), dd => dist_data(tp))

          ! child TPs and IDs of parent face neighbors .........................

          id_child_face =  0
          tp_child_face = -1
          do j = 1, 6
            if (element % face(j) % n_neighbor /= 1) cycle
            k = element % face(j) % i_neighbor
            l = element % neighbor(k) % id
            m = element % neighbor(k) % component
            select case(m)
            case(1)
              call element % AlignFromNeighborFace &
                       (j, k, id_child(1,:,:,l), id_child_face(:,:,j))
            case(2)
              call element % AlignFromNeighborFace &
                       (j, k, id_child(2,:,:,l), id_child_face(:,:,j))
            case(3)
              call element % AlignFromNeighborFace &
                       (j, k, id_child(:,1,:,l), id_child_face(:,:,j))
            case(4)
              call element % AlignFromNeighborFace &
                       (j, k, id_child(:,2,:,l), id_child_face(:,:,j))
            case(5)
              call element % AlignFromNeighborFace &
                       (j, k, id_child(:,:,1,l), id_child_face(:,:,j))
            case(6)
              call element % AlignFromNeighborFace &
                       (j, k, id_child(:,:,2,l), id_child_face(:,:,j))
            end select
            tp_child_face(j) = tp_child(l)
          end do

          ! child TPs and IDs of parent edge neighbors .........................

          id_child_edge =  0
          tp_child_edge = -1
          do j = 1, 12
            do i = 1, element % edge(j) % n_neighbor
              k = element % edge(j) % i_neighbor + i - 1
              l = element % neighbor(k) % id
              m = element % neighbor(k) % component
              select case(m)
              case(7) ! edge 1
                call element % AlignFromNeighborEdge &
                         (j, k, id_child(:,1,1,l), id_child_edge(:,i,j))
              case(8) ! edge 2
                call element % AlignFromNeighborEdge &
                         (j, k, id_child(:,2,1,l), id_child_edge(:,i,j))
              case(9) ! edge 3
                call element % AlignFromNeighborEdge &
                         (j, k, id_child(:,1,2,l), id_child_edge(:,i,j))
              case(10) ! edge 4
                call element % AlignFromNeighborEdge &
                         (j, k, id_child(:,2,2,l), id_child_edge(:,i,j))
              case(11) ! edge 5
                call element % AlignFromNeighborEdge &
                         (j, k, id_child(1,:,1,l), id_child_edge(:,i,j))
              case(12) ! edge 6
                call element % AlignFromNeighborEdge &
                         (j, k, id_child(2,:,1,l), id_child_edge(:,i,j))
              case(13) ! edge 7
                call element % AlignFromNeighborEdge &
                         (j, k, id_child(1,:,2,l), id_child_edge(:,i,j))
              case(14) ! edge 8
                call element % AlignFromNeighborEdge &
                         (j, k, id_child(2,:,2,l), id_child_edge(:,i,j))
              case(15) ! edge 9
                call element % AlignFromNeighborEdge &
                         (j, k, id_child(1,1,:,l), id_child_edge(:,i,j))
              case(16) ! edge 10
                call element % AlignFromNeighborEdge &
                         (j, k, id_child(2,1,:,l), id_child_edge(:,i,j))
              case(17) ! edge 11
                call element % AlignFromNeighborEdge &
                         (j, k, id_child(1,2,:,l), id_child_edge(:,i,j))
              case(18) ! edge 12
                call element % AlignFromNeighborEdge &
                         (j, k, id_child(2,2,:,l), id_child_edge(:,i,j))
              end select
              tp_child_edge(i,j) = tp_child(l)
            end do
          end do

          ! child TPs and IDs of parent vertex neighbors .......................

          id_child_vert =  0
          tp_child_vert = -1
          do j = 1, 8
            do i = 1, element % vertex(j) % n_neighbor
              k = element % vertex(j) % i_neighbor + i - 1
              l = element % neighbor(k) % id
              m = element % neighbor(k) % component
              select case(m - 18)
              case(19) ! vertex 1
                id_child_vert(i,j) = id_child(1,1,1,l)
              case(20) ! vertex 2
                id_child_vert(i,j) = id_child(2,1,1,l)
              case(21) ! vertex 3
                id_child_vert(i,j) = id_child(1,2,1,l)
              case(22) ! vertex 4
                id_child_vert(i,j) = id_child(2,2,1,l)
              case(23) ! vertex 5
                id_child_vert(i,j) = id_child(1,1,2,l)
              case(24) ! vertex 6
                id_child_vert(i,j) = id_child(2,1,2,l)
              case(25) ! vertex 7
                id_child_vert(i,j) = id_child(1,2,2,l)
              case(26) ! vertex 8
                id_child_vert(i,j) = id_child(2,2,2,l)
              end select
              tp_child_vert(i,j) = tp_child(l)
            end do
          end do

          ! child 1,1,1 ........................................................

          if (id_child(1,1,1,e) > 0) then

            ! rank in children contributed to partition tp
            ce = ce_part(tp) + 1

            ! element data
            dd % element(ce) % id = id_child(1,1,1,e)
            dd % element(ce) % adaptation % parent_proc = parent % proc
            dd % element(ce) % adaptation % parent_id   = e

            ! face data
            dd % element(ce) % face(1) % boundary = element % face(1) % boundary
            dd % element(ce) % face(3) % boundary = element % face(3) % boundary
            dd % element(ce) % face(5) % boundary = element % face(5) % boundary

            ! initialization of neighbor counters
            cn = cn_part(tp)
            nn = 0

            ! neighbors at faces 1:6, update cn and nn
            call FaceNeighbor_NF(1, id_child_face(1,1,1), tp_child_face(1))
            call FaceNeighbor_SI(2, id_child(2,1,1,e)   , tp              )
            call FaceNeighbor_NF(3, id_child_face(1,1,3), tp_child_face(3))
            call FaceNeighbor_SI(4, id_child(1,2,1,e)   , tp              )
            call FaceNeighbor_NF(5, id_child_face(1,1,5), tp_child_face(5))
            call FaceNeighbor_SI(6, id_child(1,1,2,e)   , tp              )

            ! neighbors at edges 1:12
            call EdgeNeighbor_NE( 1,    id_child_edge(1,:,1), tp_child_edge(:,1))
            call EdgeNeighbor_NF( 2, 5, id_child_face(1,2,5), tp_child_face(5)  )
            call EdgeNeighbor_NF( 3, 3, id_child_face(1,2,3), tp_child_face(3)  )
            call EdgeNeighbor_SI( 4,    id_child(1,2,2,e)   , tp                )
            call EdgeNeighbor_NE( 5,    id_child_edge(1,:,5), tp_child_edge(:,5))
            call EdgeNeighbor_NF( 6, 5, id_child_face(2,1,5), tp_child_face(5)  )
            call EdgeNeighbor_NF( 7, 1, id_child_face(1,2,1), tp_child_face(1)  )
            call EdgeNeighbor_SI( 8,    id_child(2,1,2,e)   , tp                )
            call EdgeNeighbor_NE( 9,    id_child_edge(1,:,9), tp_child_edge(:,9))
            call EdgeNeighbor_NF(10, 3, id_child_face(2,1,3), tp_child_face(3)  )
            call EdgeNeighbor_NF(11, 1, id_child_face(2,1,1), tp_child_face(1)  )
            call EdgeNeighbor_SI(12,    id_child(2,2,1,e)   , tp                )

            ! neighbors at vertices 1:8
            call VertNeighbor_NV(1,    id_child_vert(  :,1), tp_child_vert(:,1)   )
            call VertNeighbor_NE(2, 1, id_child_edge(2,:,1), tp_child_edge(:,1), 1)
            call VertNeighbor_NE(3, 5, id_child_edge(2,:,5), tp_child_edge(:,5), 1)
            call VertNeighbor_NF(4, 5, id_child_face(2,2,5), tp_child_face(5)     )
            call VertNeighbor_NE(5, 9, id_child_edge(2,:,9), tp_child_edge(:,9), 1)
            call VertNeighbor_NF(6, 3, id_child_face(2,2,3), tp_child_face(3)     )
            call VertNeighbor_NF(7, 1, id_child_face(2,2,1), tp_child_face(1)     )
            call VertNeighbor_SI(8,    id_child(2,2,2,e)   , tp                   )

            ! element geometry

            ce_part(tp) = ce
            cn_part(tp) = cn
          end if

        end associate
      end do

    end associate

  contains

    !---------------------------------------------------------------------------
    !> Face neighbor at parent neighbor face

    subroutine FaceNeighbor_NF(ef, nb_id, nb_tp)
      integer, intent(in) :: ef    !< element face
      integer, intent(in) :: nb_id !< neighbor ID
      integer, intent(in) :: nb_tp !< neighbor target partition

      integer :: ii

      if (nb_id <= 0) return

      cn = cn + 1
      nn = nn + 1

      dd % element(ce) % face(ef) % n_neighbor = 1
      dd % element(ce) % face(ef) % i_neighbor = nn

      ii = element % face(ef) % i_neighbor

      dd % neighbor_id          (cn) = nb_id
      dd % neighbor_part        (cn) = nb_tp
      dd % neighbor_component   (cn) = element % neighbor(ii) % component
      dd % neighbor_orientation (cn) = element % neighbor(ii) % orientation

    end subroutine FaceNeighbor_NF

    !---------------------------------------------------------------------------
    !> Sibling face neighbor

    subroutine FaceNeighbor_SI(ef, nb_id, nb_tp)
      integer, intent(in) :: ef    !< element face
      integer, intent(in) :: nb_id !< neighbor ID
      integer, intent(in) :: nb_tp !< neighbor target partition

      if (nb_id <= 0) return

      cn = cn + 1
      nn = nn + 1

      dd % element(ce) % face(ef) % n_neighbor = 1
      dd % element(ce) % face(ef) % i_neighbor = nn

      dd % neighbor_id   (cn) = nb_id
      dd % neighbor_part (cn) = nb_tp

      select case(ef)
      case(1)
        dd % neighbor_component (cn) = 2
      case(2)
        dd % neighbor_component (cn) = 1
      case(3)
        dd % neighbor_component (cn) = 4
      case(4)
        dd % neighbor_component (cn) = 3
      case(5)
        dd % neighbor_component (cn) = 6
      case(6)
        dd % neighbor_component (cn) = 5
      end select

      dd % neighbor_orientation (cn) = 12  ! always aligned !

    end subroutine FaceNeighbor_SI

    !---------------------------------------------------------------------------
    !> Edge neighbors at parent neighbor edge

    subroutine EdgeNeighbor_NE(ee, nb_id, nb_tp)
      integer, intent(in) :: ee       !< element edge
      integer, intent(in) :: nb_id(:) !< neighbor IDs
      integer, intent(in) :: nb_tp(:) !< neighbor target partitions

      integer :: ii, jj, nn_ce

      nn_ce = count(nb_id > 0)
      if (nn_ce == 0) return

      dd % element(ce) % edge(ee) % n_neighbor = nn_ce
      dd % element(ce) % edge(ee) % i_neighbor = nn + 1

      ii = element % edge(ee) % i_neighbor
      do jj = 1, element % edge(ee) % n_neighbor
        if (nb_id(jj) > 0) then
          cn = cn + 1
          dd % neighbor_id          (cn) = nb_id(jj)
          dd % neighbor_part        (cn) = nb_tp(jj)
          dd % neighbor_component   (cn) = element % neighbor(ii) % component
          dd % neighbor_orientation (cn) = element % neighbor(ii) % orientation
        end if
        ii = ii + 1
      end do
      nn = nn + nn_ce

    end subroutine EdgeNeighbor_NE

    !---------------------------------------------------------------------------
    !> Edge neighbor at parent neighbor face

    subroutine EdgeNeighbor_NF(ee, ef, nb_id, nb_tp)
      integer, intent(in) :: ee    !< element edge
      integer, intent(in) :: ef    !< element face
      integer, intent(in) :: nb_id !< neighbor ID
      integer, intent(in) :: nb_tp !< neighbor target partition

      integer(IXS) :: component, orientation
      integer :: xi_o(3), xi_t(3)
      integer :: ii, nee

      if (nb_id <= 0) return

      cn = cn + 1
      nn = nn + 1

      dd % element(ce) % edge(ee) % n_neighbor = 1
      dd % element(ce) % edge(ee) % i_neighbor = nn

      ! adopt orientation from parent face neighbor .............................

      ii = element % face(ef) % i_neighbor
      orientation = element % neighbor(ii) % orientation

      ! identify matching adjacent child edge ...................................

      ! 1: midpoint (ξ,η,ζ) coordinates of aligned neighbor edge → xi_o
      select case(ee)
      case( 1)
        xi_o = [  0,  1,  1 ]
      case( 2)
        xi_o = [  0, -1,  1 ]
      case( 3)
        xi_o = [  0,  1, -1 ]
      case( 4)
        xi_o = [  0, -1, -1 ]
      case( 5)
        xi_o = [  1,  0,  1 ]
      case( 6)
        xi_o = [ -1,  0,  1 ]
      case( 7)
        xi_o = [  1,  0, -1 ]
      case( 8)
        xi_o = [ -1,  0, -1 ]
      case( 9)
        xi_o = [  1,  1,  0 ]
      case(10)
        xi_o = [ -1,  1,  0 ]
      case(11)
        xi_o = [  1, -1,  0 ]
      case(12)
        xi_o = [ -1, -1,  0 ]
      end select

      ! 2: transformed neighbor edge coordinates → xi_t
      call TransformIndex(orientation, 3, o0 = -1, o = xi_o, t0 = -1, t = xi_t)

      ! 3: ID of neighbor element edge
      if (xi_o(1) == 0) then
        nee = 1 + (xi_o(2) + 1)/2 + (xi_o(3) + 1)
      else if (xi_o(2) == 0) then
        nee = 5 + (xi_o(1) + 1)/2 + (xi_o(3) + 1)
      else if (xi_o(3) == 0) then
        nee = 9 + (xi_o(1) + 1)/2 + (xi_o(2) + 1)
      else ! should not happen, trigger error by providing invalid number
        nee = -100
      end if

      ! 4: neighbor component ID
      component = int(6 + nee, IXS)

      ! assign neighbor data ...................................................

      dd % neighbor_id          (cn) = nb_id
      dd % neighbor_part        (cn) = nb_tp
      dd % neighbor_component   (cn) = component
      dd % neighbor_orientation (cn) = orientation

    end subroutine EdgeNeighbor_NF

    !---------------------------------------------------------------------------
    !> Sibling edge neighbor

    subroutine EdgeNeighbor_SI(ee, nb_id, nb_tp)
      integer, intent(in) :: ee    !< element edge
      integer, intent(in) :: nb_id !< neighbor ID
      integer, intent(in) :: nb_tp !< neighbor target partition

      if (nb_id <= 0) return

      cn = cn + 1
      nn = nn + 1

      dd % element(ce) % edge(ee) % n_neighbor = 1
      dd % element(ce) % edge(ee) % i_neighbor = nn

      dd % neighbor_id   (cn) = nb_id
      dd % neighbor_part (cn) = nb_tp

      ! ID of coupled component = ID of opposite element edge
      select case(ee)
      case(1:4)
        dd % neighbor_component (cn) = 6 +  5 - ee
      case(5:8)
        dd % neighbor_component (cn) = 6 + 13 - ee
      case(9:12)
        dd % neighbor_component (cn) = 6 + 21 - ee
      end select

      dd % neighbor_orientation (cn) = 12  ! always aligned !

    end subroutine EdgeNeighbor_SI

    !---------------------------------------------------------------------------
    !> Vertex neighbors at parent neighbor vertex

    subroutine VertNeighbor_NV(ev, nb_id, nb_tp)
      integer, intent(in) :: ev       !< element vertex
      integer, intent(in) :: nb_id(:) !< neighbor IDs
      integer, intent(in) :: nb_tp(:) !< neighbor target partitions

       integer :: ii, jj, nn_ce

      nn_ce = count(nb_id > 0)
      if (nn_ce== 0) return

      dd % element(ce) % vertex(ev) % n_neighbor = nn_ce
      dd % element(ce) % vertex(ev) % i_neighbor = nn + 1

      ii = element % vertex(ev) % i_neighbor
      do jj = 1, element % vertex(ev) % n_neighbor
        if (nb_id(jj) > 0) then
          cn = cn + 1
          dd % neighbor_id          (cn) = nb_id(jj)
          dd % neighbor_part        (cn) = nb_tp(jj)
          dd % neighbor_component   (cn) = element % neighbor(ii) % component
          dd % neighbor_orientation (cn) = element % neighbor(ii) % orientation
        end if
        ii = ii + 1
      end do
      nn = nn + nn_ce

    end subroutine VertNeighbor_NV

    !---------------------------------------------------------------------------
    !> Vertex neighbors at parent neighbor edge

    subroutine VertNeighbor_NE(ev, ee, nb_id, nb_tp, ec_rank)
      integer, intent(in) :: ev       !< element vertex
      integer, intent(in) :: ee       !< element edge
      integer, intent(in) :: nb_id(:) !< neighbor IDs
      integer, intent(in) :: nb_tp(:) !< neighbor target partitions
      integer, intent(in) :: ec_rank  !< edge child rank: first (1) or second (2)

      integer :: ii, jj, nee, nev, nn_ce

      nn_ce = count(nb_id > 0)
      if (nn_ce == 0) return

      dd % element(ce) % vertex(ev) % n_neighbor = nn_ce
      dd % element(ce) % vertex(ev) % i_neighbor = nn + 1

      ii = element % edge(ee) % i_neighbor
      do jj = 1, element % edge(ee) % n_neighbor
        if (nb_id(jj) > 0) then
          cn = cn + 1

          ! parent neighbor edge
          nee = ElementEdgeID(element % neighbor(ii) % component)

          ! matching neighbor vertex
          if (element % NeigborEdgeIsAligned(ee, ii)) then
            nev = V_EDGE(ec_rank, nee)
          else
            nev = V_EDGE(3 - ec_rank, nee)  ! edge vertex order reversed !
          end if

          ! assign neighbor data .................................................

          dd % neighbor_id          (cn) = nb_id(jj)
          dd % neighbor_part        (cn) = nb_tp(jj)
          dd % neighbor_component   (cn) = 18 + nev
          dd % neighbor_orientation (cn) = element % neighbor(ii) % orientation

        end if
        ii = ii + 1
      end do
      nn = nn + nn_ce

    end subroutine VertNeighbor_NE

    !---------------------------------------------------------------------------
    !> Vertex neighbor at parent neighbor face

    subroutine VertNeighbor_NF(ev, ef, nb_id, nb_tp)
      integer, intent(in) :: ev    !< element vertex
      integer, intent(in) :: ef    !< element face
      integer, intent(in) :: nb_id !< neighbor ID
      integer, intent(in) :: nb_tp !< neighbor target partition

      integer(IXS) :: orientation
      integer :: xi_o(3), xi_t(3)
      integer :: ii, nev

      if (nb_id <= 0) return

      cn = cn + 1
      nn = nn + 1

      dd % element(ce) % vertex(ev) % n_neighbor = 1
      dd % element(ce) % vertex(ev) % i_neighbor = nn

      ! adopt orientation from parent face neighbor .............................

      ii = element % face(ef) % i_neighbor
      orientation = element % neighbor(ii) % orientation

      ! identify matching adjacent child vertex ...............................

      ! 1: (ξ,η,ζ) coordinates of aligned neighbor vertex → xi_o
      select case(ev)
      case( 1)
        xi_o = [  1,  1,  1 ]
      case(2)
        xi_o = [ -1,  1,  1 ]
      case(3)
        xi_o = [  1, -1,  1 ]
      case(4)
        xi_o = [ -1, -1,  1 ]
      case(5)
        xi_o = [  1,  1, -1 ]
      case( 6)
        xi_o = [ -1,  1, -1 ]
      case(7)
        xi_o = [  1, -1, -1 ]
      case(8)
        xi_o = [ -1, -1, -1 ]
      end select

      ! 2: transformed neighbor vertex coordinates → xi_t
      call TransformIndex(orientation, 3, o0 = -1, o = xi_o, t0 = -1, t = xi_t)

      ! 3: ID of neighbor element vertex
      nev = 1 + (xi_o(1) + 1)/2 + (xi_o(2) + 1) + (xi_o(3) + 1)*2

      ! assign neighbor data ...................................................

      dd % neighbor_id          (cn) = nb_id
      dd % neighbor_part        (cn) = nb_tp
      dd % neighbor_component   (cn) = 18 + nev
      dd % neighbor_orientation (cn) = orientation

    end subroutine VertNeighbor_NF

    !---------------------------------------------------------------------------
    !> Sibling vertex neighbor

    subroutine VertNeighbor_SI(ev, nb_id, nb_tp)
      integer, intent(in) :: ev    !< element vertex
      integer, intent(in) :: nb_id !< neighbor ID
      integer, intent(in) :: nb_tp !< neighbor target partition

      if (nb_id <= 0) return

      cn = cn + 1
      nn = nn + 1

      dd % element(ce) % vertex(ev) % n_neighbor = 1
      dd % element(ce) % vertex(ev) % i_neighbor = nn

      dd % neighbor_id          (cn) = nb_id
      dd % neighbor_part        (cn) = nb_tp
      dd % neighbor_component   (cn) = 18 + 9 - ev ! diagonally opposite vertex
      dd % neighbor_orientation (cn) = 12          ! always aligned !

    end subroutine VertNeighbor_SI

    !---------------------------------------------------------------------------

  end subroutine BuildDistributionData


  !=============================================================================

end module Child_Mesh_Adaptation__3D
