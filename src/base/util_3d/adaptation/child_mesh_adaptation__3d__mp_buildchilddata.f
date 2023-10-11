submodule(Child_Mesh_Adaptation__3D) MP_BuildChildData
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Build child distribution data for sending from parent to child partitions

  module subroutine BuildChildData(parent, child, map, child_data)

    class(Mesh_3D), intent(in) :: parent
    class(Mesh_3D), intent(in) :: child
    class(ChildDistributionMap_3D),    intent(in)  :: map
    class(ElementDistributionData_3D), intent(out) :: child_data(0:)

    !---------------------------------------------------------------------------
    ! internal variables

    ! start position of element contributions to the respective target partition
    integer, allocatable, save :: start_child_element(:)  ! child elements
    integer, allocatable, save :: start_child_neighbor(:) ! child neighbor data
    integer, allocatable, save :: start_child_point(:)    ! child geometry data

    integer, allocatable, save :: id_child_face(:,:,:) ! IDs of neighbor children
    integer, allocatable, save :: id_child_edge(:,:,:) ! at parent element faces,
    integer, allocatable, save :: id_child_vert(:,:)   ! edges and vertices

    integer, allocatable, save :: tp_child_face(:)     ! TPs of neighbor children
    integer, allocatable, save :: tp_child_edge(:,:)   ! at parent element faces,
    integer, allocatable, save :: tp_child_vert(:,:)   ! edges and vertices

    integer, allocatable, save :: ce_part(:)  ! element counters
    integer, allocatable, save :: cn_part(:)  ! neighbor counters
    integer, allocatable, save :: cp_part(:)  ! element point counters

    type(StandardOperators_1D),          allocatable, save :: eop(:)
    type(ParentToChildInterpolation_1D), allocatable, save :: iop(:)

    integer, save :: nn_edge, nn_vert
    integer, save :: po_min, po_max

    integer :: nc, nn, np, po, tp
    integer :: e, i, j, k, l, m
    integer :: ce, cn, cp
    integer :: child_mark
!### CHECK
write(*,'(99(G0,1X))') '# BCD  0',', parent%proc',parent%proc
!### CHECK END

    if (parent % part < 0) return

    associate(id_child => map % id_child, tp_child => map % tp_child)

      !-------------------------------------------------------------------------
      ! initialization

      !$omp master

      allocate(start_child_element  ( parent%n_elem ), source = -1)
      allocate(start_child_neighbor ( parent%n_elem ), source = -1)
      allocate(start_child_point    ( parent%n_elem ), source = -1)

      allocate(ce_part( 0:map%n_parts-1 ))
      allocate(cn_part( 0:map%n_parts-1 ))
      allocate(cp_part( 0:map%n_parts-1 ))

      ! basic initialization ...................................................

      do tp = 0, map%n_parts-1
        child_data(tp) % proc       = child % proc_part(tp)
        child_data(tp) % comm       = child % comm_world
        child_data(tp) % n_elem     = 0
        child_data(tp) % n_neighbor = 0
        child_data(tp) % n_point    = 0
      end do

      ! max dimensions
      nn_edge = 0  ! max number of neighbors per child element edge
      nn_vert = 0  ! max number of neighbors per child element vertex

      ! polynomial degree of geometry
      po_min = parent % p_geom
      po_max = 1

      ! counters
      ce_part = 0  ! count of child elements per partition
      cn_part = 0  ! upper bound for count of child neighbors
      cp_part = 0  ! upper bound for count of child geometry points

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
            i = max( maxval( element % face   % n_neighbor ), 1)
            j = max( maxval( element % edge   % n_neighbor ), i)
            k = max( maxval( element % vertex % n_neighbor ), j)

            ! max number of neighbors per child
            nn = 6 * i + 12 * j + 8 * k

            ! number of element points per child
            np = (element % geometry % po + 1) ** 3

            ! start entries
            start_child_element  (e) = ce_part(tp) + 1
            start_child_neighbor (e) = cn_part(tp) + 1
            start_child_point    (e) = cp_part(tp) + 1

            ! update counters
            ce_part(tp) = ce_part(tp) + nc
            cn_part(tp) = cn_part(tp) + nc * nn
            cp_part(tp) = cp_part(tp) + nc * np

            ! update dimensions
            nn_edge = max(j, nn_edge)
            nn_vert = max(k, nn_vert)

            ! update bounds of geometry polynomial degree
            po_min = min(po_min, element % geometry % po)
            po_max = max(po_max, element % geometry % po)

          end if
        end associate
      end do
!### CHECK
write(*,'(99(G0,1X))') '# BCD  1',', parent%proc',parent%proc
!### CHECK END

      ! allocate components
      do tp = 0, map%n_parts-1

        nc = ce_part(tp)
        nn = cn_part(tp)
        np = cp_part(tp)

        child_data(tp) % n_elem     = nc
        child_data(tp) % n_neighbor = nn
        child_data(tp) % n_point    = np

        allocate( child_data(tp) % element              ( nc    ) )
        allocate( child_data(tp) % start_neighbor       ( nc    ) )
        allocate( child_data(tp) % neighbor_id          ( nn    ) )
        allocate( child_data(tp) % neighbor_part        ( nn    ) )
        allocate( child_data(tp) % neighbor_component   ( nn    ) )
        allocate( child_data(tp) % neighbor_orientation ( nn    ) )
        allocate( child_data(tp) % start_point          ( nc    ) )
        allocate( child_data(tp) % x_e                  ( np, 3 ) )

      end do
!### CHECK
write(*,'(99(G0,1X))') '# BCD  2',', parent%proc',parent%proc
!### CHECK END

      ! auxiliary arrays
      allocate( id_child_face (2,2,6)        , tp_child_face (6)          )
      allocate( id_child_edge (2,nn_edge,12) , tp_child_edge (nn_edge,12) )
      allocate( id_child_vert (nn_vert,8)    , tp_child_vert (nn_vert,8)  )

      ! interpolation operators
      allocate( eop(po_min:po_max) )
      allocate( iop(po_min:po_max) )
      do po = po_min, po_max
        eop(po) = StandardOperators_1D(po, basis = 'L', no_vdm = .true.)
        iop(po) = ParentToChildInterpolation_1D(eop(po))
      end do
!### CHECK
write(*,'(99(G0,1X))') '# BCD  3',', parent%proc',parent%proc
!### CHECK END

      !$omp end master
      !$omp end barrier
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      ! generate child element data

      !$omp do
      do e = 1, parent % n_elem

        tp = tp_child(e)
        if (tp < 0) cycle

        associate(element => parent % element(e), cd_tp => child_data(tp))

          ! preliminaries ......................................................

          ! set counter offsets to start - 1
          ce = start_child_element  (e) - 1
          cn = start_child_neighbor (e) - 1
          cp = start_child_point    (e) - 1

          po = element % geometry % po
          np = (po + 1)**3

          ! set child mark to old proc ID elements are retained
          if ( element % adaptation % refinement == 100 .and.   &
               element % adaptation % mark       == 100       ) &
          then
            child_mark = element % adaptation % child_proc
          else
            child_mark = -1
          end if
!### CHECK
!! if (parent%part == 1 .and. parent%n_elem == 40 .and. e == 25) then
!! print '(99(G0,1X))', 'L2 P1 E25','refinement =',element % adaptation % refinement
!! print '(99(G0,1X))', 'L2 P1 E25','mark       =',element % adaptation % mark
!! print '(99(G0,1X))', 'L2 P1 E25','child_mark =',child_mark
!! end if
!### CHECK END

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

            ce = ce + 1

            ! element data
            cd_tp % element(ce) % id = id_child(1,1,1,e)
            cd_tp % element(ce) % cluster_oct = 1
            cd_tp % element(ce) % frozen = element % adaptation % mark < 100
            cd_tp % element(ce) % adaptation % parent_proc = parent % proc
            cd_tp % element(ce) % adaptation % parent_id   = e
            cd_tp % element(ce) % adaptation % mark        = child_mark

            ! face data
            cd_tp % element(ce) % face(1) % boundary = element % face(1) % boundary
            cd_tp % element(ce) % face(3) % boundary = element % face(3) % boundary
            cd_tp % element(ce) % face(5) % boundary = element % face(5) % boundary

            ! start index of child neighbor data
            cd_tp % start_neighbor(ce) = cn + 1

            ! initialization of neighbor counter
            nn = 0

            ! neighbors at faces 1:6
            call FaceNeighbor_NF(cd_tp, element, ef=1, i1=1, i2=1)
            call FaceNeighbor_SI(cd_tp, ef=2, i1=2, i2=1, i3=1)
            call FaceNeighbor_NF(cd_tp, element, ef=3, i1=1, i2=1)
            call FaceNeighbor_SI(cd_tp, ef=4, i1=1, i2=2, i3=1)
            call FaceNeighbor_NF(cd_tp, element, ef=5, i1=1, i2=1)
            call FaceNeighbor_SI(cd_tp, ef=6, i1=1, i2=1, i3=2)

            ! neighbors at edges 1:12
            call EdgeNeighbor_NE(cd_tp, element, ee= 1, i1=1)
            call EdgeNeighbor_NF(cd_tp, element, ee= 2, ef=5, i1=1, i2=2)
            call EdgeNeighbor_NF(cd_tp, element, ee= 3, ef=3, i1=1, i2=2)
            call EdgeNeighbor_SI(cd_tp, ee= 4, i1=1, i2=2, i3=2)
            call EdgeNeighbor_NE(cd_tp, element, ee= 5, i1=1)
            call EdgeNeighbor_NF(cd_tp, element, ee= 6, ef=5, i1=2, i2=1)
            call EdgeNeighbor_NF(cd_tp, element, ee= 7, ef=1, i1=1, i2=2)
            call EdgeNeighbor_SI(cd_tp, ee= 8, i1=2, i2=1, i3=2)
            call EdgeNeighbor_NE(cd_tp, element, ee= 9, i1=1)
            call EdgeNeighbor_NF(cd_tp, element, ee=10, ef=3, i1=2, i2=1)
            call EdgeNeighbor_NF(cd_tp, element, ee=11, ef=1, i1=2, i2=1)
            call EdgeNeighbor_SI(cd_tp, ee=12, i1=2, i2=2, i3=1)

            ! neighbors at vertices 1:8
            call VertNeighbor_NV(cd_tp, element, ev=1)
            call VertNeighbor_NE(cd_tp, element, ev=2, ee=1, i1=2)
            call VertNeighbor_NE(cd_tp, element, ev=3, ee=5, i1=2)
            call VertNeighbor_NF(cd_tp, element, ev=4, ef=5, i1=2, i2=2)
            call VertNeighbor_NE(cd_tp, element, ev=5, ee=9, i1=2)
            call VertNeighbor_NF(cd_tp, element, ev=6, ef=3, i1=2, i2=2)
            call VertNeighbor_NF(cd_tp, element, ev=7, ef=1, i1=2, i2=2)
            call VertNeighbor_SI(cd_tp, ev=8, i1=2, i2=2, i3=2)

            ! geometry
            cd_tp % element(ce) % geometry % po = po
            cd_tp % start_point(ce) = cp + 1
            call InterpolateCoords( po = po                         &
                                  , A1  = iop(po) % A(:,:,1)        &
                                  , A2  = iop(po) % A(:,:,1)        &
                                  , A3  = iop(po) % A(:,:,1)        &
                                  , xp  = element % geometry % x_e  &
                                  , xc1 = cd_tp % x_e(cp+1:cp+np,1) &
                                  , xc2 = cd_tp % x_e(cp+1:cp+np,2) &
                                  , xc3 = cd_tp % x_e(cp+1:cp+np,3) )
            cp = cp + np

          end if

          ! child 2,1,1 ........................................................

          if (id_child(2,1,1,e) > 0) then

            ce = ce + 1

            ! element data
            cd_tp % element(ce) % id = id_child(2,1,1,e)
            cd_tp % element(ce) % cluster_oct = 2
            cd_tp % element(ce) % frozen = element % adaptation % mark < 100
            cd_tp % element(ce) % adaptation % parent_proc = parent % proc
            cd_tp % element(ce) % adaptation % parent_id   = e
            cd_tp % element(ce) % adaptation % mark        = child_mark

            ! face data
            cd_tp % element(ce) % face(2) % boundary = element % face(2) % boundary
            cd_tp % element(ce) % face(3) % boundary = element % face(3) % boundary
            cd_tp % element(ce) % face(5) % boundary = element % face(5) % boundary

            cd_tp % start_neighbor(ce) = cn + 1

            nn = 0

            ! neighbors at faces 1:6
            call FaceNeighbor_SI(cd_tp, ef=1, i1=1, i2=1, i3=1)
            call FaceNeighbor_NF(cd_tp, element, ef=2, i1=1, i2=1)
            call FaceNeighbor_NF(cd_tp, element, ef=3, i1=2, i2=1)
            call FaceNeighbor_SI(cd_tp, ef=4, i1=2, i2=2, i3=1)
            call FaceNeighbor_NF(cd_tp, element, ef=5, i1=2, i2=1)
            call FaceNeighbor_SI(cd_tp, ef=6, i1=2, i2=1, i3=2)

            ! neighbors at edges 1:12
            call EdgeNeighbor_NE(cd_tp, element, ee= 1, i1=2)
            call EdgeNeighbor_NF(cd_tp, element, ee= 2, ef=5, i1=2, i2=2)
            call EdgeNeighbor_NF(cd_tp, element, ee= 3, ef=3, i1=2, i2=2)
            call EdgeNeighbor_SI(cd_tp, ee= 4, i1=2, i2=2, i3=2)
            call EdgeNeighbor_NF(cd_tp, element, ee= 5, ef=5, i1=1, i2=1)
            call EdgeNeighbor_NE(cd_tp, element, ee= 6, i1=1)
            call EdgeNeighbor_SI(cd_tp, ee= 7, i1=1, i2=1, i3=2)
            call EdgeNeighbor_NF(cd_tp, element, ee= 8, ef=2, i1=1, i2=2)
            call EdgeNeighbor_NF(cd_tp, element, ee= 9, ef=3, i1=1, i2=1)
            call EdgeNeighbor_NE(cd_tp, element, ee=10, i1=1)
            call EdgeNeighbor_SI(cd_tp, ee=11, i1=1, i2=2, i3=1)
            call EdgeNeighbor_NF(cd_tp, element, ee=12, ef=2, i1=2, i2=1)

            ! neighbors at vertices 1:8
            call VertNeighbor_NE(cd_tp, element, ev=1, ee=1,  i1=1)
            call VertNeighbor_NV(cd_tp, element, ev=2)
            call VertNeighbor_NF(cd_tp, element, ev=3, ef=5,  i1=1, i2=2)
            call VertNeighbor_NE(cd_tp, element, ev=4, ee=6,  i1=2)
            call VertNeighbor_NF(cd_tp, element, ev=5, ef=3,  i1=1, i2=2)
            call VertNeighbor_NE(cd_tp, element, ev=6, ee=10, i1=2)
            call VertNeighbor_SI(cd_tp, ev=7, i1=1,  i2=2, i3=2)
            call VertNeighbor_NF(cd_tp, element, ev=8, ef=2,  i1=2, i2=2)

            ! geometry
            cd_tp % element(ce) % geometry % po = po
            cd_tp % start_point(ce) = cp + 1
            call InterpolateCoords( po = po                         &
                                  , A1  = iop(po) % A(:,:,2)        &
                                  , A2  = iop(po) % A(:,:,1)        &
                                  , A3  = iop(po) % A(:,:,1)        &
                                  , xp  = element % geometry % x_e  &
                                  , xc1 = cd_tp % x_e(cp+1:cp+np,1) &
                                  , xc2 = cd_tp % x_e(cp+1:cp+np,2) &
                                  , xc3 = cd_tp % x_e(cp+1:cp+np,3) )
            cp = cp + np

          end if

          ! child 1,2,1 ........................................................

          if (id_child(1,2,1,e) > 0) then

            ce = ce + 1

            ! element data
            cd_tp % element(ce) % id = id_child(1,2,1,e)
            cd_tp % element(ce) % cluster_oct = 3
            cd_tp % element(ce) % frozen = element % adaptation % mark < 100
            cd_tp % element(ce) % adaptation % parent_proc = parent % proc
            cd_tp % element(ce) % adaptation % parent_id   = e
            cd_tp % element(ce) % adaptation % mark        = child_mark

            ! face data
            cd_tp % element(ce) % face(1) % boundary = element % face(1) % boundary
            cd_tp % element(ce) % face(4) % boundary = element % face(4) % boundary
            cd_tp % element(ce) % face(5) % boundary = element % face(5) % boundary

            cd_tp % start_neighbor(ce) = cn + 1

            nn = 0

            ! neighbors at faces 1:6
            call FaceNeighbor_NF(cd_tp, element, ef=1, i1=2, i2=1)
            call FaceNeighbor_SI(cd_tp, ef=2, i1=2, i2=2, i3=1)
            call FaceNeighbor_SI(cd_tp, ef=3, i1=1, i2=1, i3=1)
            call FaceNeighbor_NF(cd_tp, element, ef=4, i1=1, i2=1)
            call FaceNeighbor_NF(cd_tp, element, ef=5, i1=1, i2=2)
            call FaceNeighbor_SI(cd_tp, ef=6, i1=1, i2=2, i3=2)

            ! neighbors at edges 1:12
            call EdgeNeighbor_NF(cd_tp, element, ee= 1, ef=5, i1=1, i2=1)
            call EdgeNeighbor_NE(cd_tp, element, ee= 2, i1=1)
            call EdgeNeighbor_SI(cd_tp, ee= 3, i1=1, i2=1, i3=2)
            call EdgeNeighbor_NF(cd_tp, element, ee= 4, ef=4, i1=1, i2=2)
            call EdgeNeighbor_NE(cd_tp, element, ee= 5, i1=2)
            call EdgeNeighbor_NF(cd_tp, element, ee= 6, ef=5, i1=2, i2=2)
            call EdgeNeighbor_NF(cd_tp, element, ee= 7, ef=1, i1=2, i2=2)
            call EdgeNeighbor_SI(cd_tp, ee= 8, i1=2, i2=2, i3=2)
            call EdgeNeighbor_NF(cd_tp, element, ee= 9, ef=1, i1=1, i2=1)
            call EdgeNeighbor_SI(cd_tp, ee=10, i1=2, i2=1, i3=1)
            call EdgeNeighbor_NE(cd_tp, element, ee=11, i1=1)
            call EdgeNeighbor_NF(cd_tp, element, ee=12, ef=4, i1=2, i2=1)

            ! neighbors at vertices 1:8
            call VertNeighbor_NE(cd_tp, element, ev=1, ee= 5, i1=1)
            call VertNeighbor_NF(cd_tp, element, ev=2, ef= 5, i1=2, i2=1)
            call VertNeighbor_NV(cd_tp, element, ev=3)
            call VertNeighbor_NE(cd_tp, element, ev=4, ee= 2, i1=2)
            call VertNeighbor_NF(cd_tp, element, ev=5, ef= 1, i1=1, i2=2)
            call VertNeighbor_SI(cd_tp, ev=6, i1= 2, i2=1, i3=2)
            call VertNeighbor_NE(cd_tp, element, ev=7, ee=11, i1=2)
            call VertNeighbor_NF(cd_tp, element, ev=8, ef= 4, i1=2, i2=2)

            ! geometry
            cd_tp % element(ce) % geometry % po = po
            cd_tp % start_point(ce) = cp + 1
            call InterpolateCoords( po = po                         &
                                  , A1  = iop(po) % A(:,:,1)        &
                                  , A2  = iop(po) % A(:,:,2)        &
                                  , A3  = iop(po) % A(:,:,1)        &
                                  , xp  = element % geometry % x_e  &
                                  , xc1 = cd_tp % x_e(cp+1:cp+np,1) &
                                  , xc2 = cd_tp % x_e(cp+1:cp+np,2) &
                                  , xc3 = cd_tp % x_e(cp+1:cp+np,3) )
            cp = cp + np

          end if

          ! child 2,2,1 ........................................................

          if (id_child(2,2,1,e) > 0) then

            ce = ce + 1

            ! element data
            cd_tp % element(ce) % id = id_child(2,2,1,e)
            cd_tp % element(ce) % cluster_oct = 4
            cd_tp % element(ce) % frozen = element % adaptation % mark < 100
            cd_tp % element(ce) % adaptation % parent_proc = parent % proc
            cd_tp % element(ce) % adaptation % parent_id   = e
            cd_tp % element(ce) % adaptation % mark        = child_mark

            ! face data
            cd_tp % element(ce) % face(2) % boundary = element % face(2) % boundary
            cd_tp % element(ce) % face(4) % boundary = element % face(4) % boundary
            cd_tp % element(ce) % face(5) % boundary = element % face(5) % boundary

            cd_tp % start_neighbor(ce) = cn + 1

            nn = 0

            ! neighbors at faces 1:6
            call FaceNeighbor_SI(cd_tp, ef=1, i1=1, i2=2, i3=1)
            call FaceNeighbor_NF(cd_tp, element, ef=2, i1=2, i2=1)
            call FaceNeighbor_SI(cd_tp, ef=3, i1=2, i2=1, i3=1)
            call FaceNeighbor_NF(cd_tp, element, ef=4, i1=2, i2=1)
            call FaceNeighbor_NF(cd_tp, element, ef=5, i1=2, i2=2)
            call FaceNeighbor_SI(cd_tp, ef=6, i1=2, i2=2, i3=2)

            ! neighbors at edges 1:12
            call EdgeNeighbor_NF(cd_tp, element, ee= 1, ef=5, i1=2, i2=1)
            call EdgeNeighbor_NE(cd_tp, element, ee= 2, i1=2)
            call EdgeNeighbor_SI(cd_tp, ee= 3, i1=2, i2=1, i3=2)
            call EdgeNeighbor_NF(cd_tp, element, ee= 4, ef=4, i1=2, i2=2)
            call EdgeNeighbor_NF(cd_tp, element, ee= 5, ef=5, i1=1, i2=2)
            call EdgeNeighbor_NE(cd_tp, element, ee= 6, i1=2)
            call EdgeNeighbor_SI(cd_tp, ee= 7, i1=1, i2=2, i3=2)
            call EdgeNeighbor_NF(cd_tp, element, ee= 8, ef=2, i1=2, i2=2)
            call EdgeNeighbor_SI(cd_tp, ee= 9, i1=1, i2=1, i3=1)
            call EdgeNeighbor_NF(cd_tp, element, ee=10, ef=2, i1=1, i2=1)
            call EdgeNeighbor_NF(cd_tp, element, ee=11, ef=4, i1=1, i2=1)
            call EdgeNeighbor_NE(cd_tp, element, ee=12, i1=1)

            ! neighbors at vertices 1:8
            call VertNeighbor_NF(cd_tp, element, ev=1, ef= 5, i1=1, i2=1)
            call VertNeighbor_NE(cd_tp, element, ev=2, ee= 6, i1=1)
            call VertNeighbor_NE(cd_tp, element, ev=3, ee= 2, i1=1)
            call VertNeighbor_NV(cd_tp, element, ev=4)
            call VertNeighbor_SI(cd_tp, ev=5, i1= 1, i2=1, i3=2)
            call VertNeighbor_NF(cd_tp, element, ev=6, ef= 2, i1=1, i2=2)
            call VertNeighbor_NF(cd_tp, element, ev=7, ef= 4, i1=1, i2=2)
            call VertNeighbor_NE(cd_tp, element, ev=8, ee=12, i1=2)

            ! geometry
            cd_tp % element(ce) % geometry % po = po
            cd_tp % start_point(ce) = cp + 1
            call InterpolateCoords( po = po                         &
                                  , A1  = iop(po) % A(:,:,2)        &
                                  , A2  = iop(po) % A(:,:,2)        &
                                  , A3  = iop(po) % A(:,:,1)        &
                                  , xp  = element % geometry % x_e  &
                                  , xc1 = cd_tp % x_e(cp+1:cp+np,1) &
                                  , xc2 = cd_tp % x_e(cp+1:cp+np,2) &
                                  , xc3 = cd_tp % x_e(cp+1:cp+np,3) )
            cp = cp + np

          end if

          ! child 1,1,2 ........................................................

          if (id_child(1,1,2,e) > 0) then

            ce = ce + 1

            ! element data
            cd_tp % element(ce) % id = id_child(1,1,2,e)
            cd_tp % element(ce) % cluster_oct = 5
            cd_tp % element(ce) % frozen = element % adaptation % mark < 100
            cd_tp % element(ce) % adaptation % parent_proc = parent % proc
            cd_tp % element(ce) % adaptation % parent_id   = e
            cd_tp % element(ce) % adaptation % mark        = child_mark

            ! face data
            cd_tp % element(ce) % face(1) % boundary = element % face(1) % boundary
            cd_tp % element(ce) % face(3) % boundary = element % face(3) % boundary
            cd_tp % element(ce) % face(6) % boundary = element % face(6) % boundary

            cd_tp % start_neighbor(ce) = cn + 1

            nn = 0

            ! neighbors at faces 1:6
            call FaceNeighbor_NF(cd_tp, element, ef=1, i1=1, i2=2)
            call FaceNeighbor_SI(cd_tp, ef=2, i1=2, i2=1, i3=2)
            call FaceNeighbor_NF(cd_tp, element, ef=3, i1=1, i2=2)
            call FaceNeighbor_SI(cd_tp, ef=4, i1=1, i2=2, i3=2)
            call FaceNeighbor_SI(cd_tp, ef=5, i1=1, i2=1, i3=1)
            call FaceNeighbor_NF(cd_tp, element, ef=6, i1=1, i2=1)

            ! neighbors at edges 1:12
            call EdgeNeighbor_NF(cd_tp, element, ee= 1, ef=3, i1=1, i2=1)
            call EdgeNeighbor_SI(cd_tp, ee= 2, i1=1, i2=2, i3=1)
            call EdgeNeighbor_NE(cd_tp, element, ee= 3, i1=1)
            call EdgeNeighbor_NF(cd_tp, element, ee= 4, ef=6, i1=1, i2=2)
            call EdgeNeighbor_NF(cd_tp, element, ee= 5, ef=1, i1=1, i2=1)
            call EdgeNeighbor_SI(cd_tp, ee= 6, i1=2, i2=1, i3=1)
            call EdgeNeighbor_NE(cd_tp, element, ee= 7, i1=1)
            call EdgeNeighbor_NF(cd_tp, element, ee= 8, ef=6, i1=2, i2=1)
            call EdgeNeighbor_NE(cd_tp, element, ee= 9, i1=2)
            call EdgeNeighbor_NF(cd_tp, element, ee=10, ef=3, i1=2, i2=2)
            call EdgeNeighbor_NF(cd_tp, element, ee=11, ef=1, i1=2, i2=2)
            call EdgeNeighbor_SI(cd_tp, ee=12, i1=2, i2=2, i3=2)

            ! neighbors at vertices 1:8
            call VertNeighbor_NE(cd_tp, element, ev=1, ee= 9, i1=1)
            call VertNeighbor_NF(cd_tp, element, ev=2, ef= 3, i1=2, i2=1)
            call VertNeighbor_NF(cd_tp, element, ev=3, ef= 1, i1=2, i2=1)
            call VertNeighbor_SI(cd_tp, ev=4, i1= 2, i2=2, i3=1)
            call VertNeighbor_NV(cd_tp, element, ev=5)
            call VertNeighbor_NE(cd_tp, element, ev=6, ee= 3, i1=2)
            call VertNeighbor_NE(cd_tp, element, ev=7, ee= 7, i1=2)
            call VertNeighbor_NF(cd_tp, element, ev=8, ef= 6, i1=2, i2=2)

            ! geometry
            cd_tp % element(ce) % geometry % po = po
            cd_tp % start_point(ce) = cp + 1
            call InterpolateCoords( po = po                         &
                                  , A1  = iop(po) % A(:,:,1)        &
                                  , A2  = iop(po) % A(:,:,1)        &
                                  , A3  = iop(po) % A(:,:,2)        &
                                  , xp  = element % geometry % x_e  &
                                  , xc1 = cd_tp % x_e(cp+1:cp+np,1) &
                                  , xc2 = cd_tp % x_e(cp+1:cp+np,2) &
                                  , xc3 = cd_tp % x_e(cp+1:cp+np,3) )
            cp = cp + np

          end if

          ! child 2,1,2 ........................................................

          if (id_child(2,1,2,e) > 0) then

            ce = ce + 1

            ! element data
            cd_tp % element(ce) % id = id_child(2,1,2,e)
            cd_tp % element(ce) % cluster_oct = 6
            cd_tp % element(ce) % frozen = element % adaptation % mark < 100
            cd_tp % element(ce) % adaptation % parent_proc = parent % proc
            cd_tp % element(ce) % adaptation % parent_id   = e
            cd_tp % element(ce) % adaptation % mark        = child_mark

            ! face data
            cd_tp % element(ce) % face(2) % boundary = element % face(2) % boundary
            cd_tp % element(ce) % face(3) % boundary = element % face(3) % boundary
            cd_tp % element(ce) % face(6) % boundary = element % face(6) % boundary

            cd_tp % start_neighbor(ce) = cn + 1

            nn = 0

            ! neighbors at faces 1:6
            call FaceNeighbor_SI(cd_tp, ef=1, i1=1, i2=1, i3=2)
            call FaceNeighbor_NF(cd_tp, element, ef=2, i1=1, i2=2)
            call FaceNeighbor_NF(cd_tp, element, ef=3, i1=2, i2=2)
            call FaceNeighbor_SI(cd_tp, ef=4, i1=2, i2=2, i3=2)
            call FaceNeighbor_SI(cd_tp, ef=5, i1=2, i2=1, i3=1)
            call FaceNeighbor_NF(cd_tp, element, ef=6, i1=2, i2=1)

            ! neighbors at edges 1:12
            call EdgeNeighbor_NF(cd_tp, element, ee= 1, ef=3, i1=2, i2=1)
            call EdgeNeighbor_SI(cd_tp, ee= 2, i1=2, i2=2, i3=1)
            call EdgeNeighbor_NE(cd_tp, element, ee= 3, i1=2)
            call EdgeNeighbor_NF(cd_tp, element, ee= 4, ef=6, i1=2, i2=2)
            call EdgeNeighbor_SI(cd_tp, ee= 5, i1=1, i2=1, i3=1)
            call EdgeNeighbor_NF(cd_tp, element, ee= 6, ef=2, i1=1, i2=1)
            call EdgeNeighbor_NF(cd_tp, element, ee= 7, ef=6, i1=1, i2=1)
            call EdgeNeighbor_NE(cd_tp, element, ee= 8, i1=1)
            call EdgeNeighbor_NF(cd_tp, element, ee= 9, ef=3, i1=1, i2=2)
            call EdgeNeighbor_NE(cd_tp, element, ee=10, i1=2)
            call EdgeNeighbor_SI(cd_tp, ee=11, i1=1, i2=2, i3=2)
            call EdgeNeighbor_NF(cd_tp, element, ee=12, ef=2, i1=2, i2=2)

            ! neighbors at vertices 1:8
            call VertNeighbor_NF(cd_tp, element, ev=1, ef= 3, i1=1, i2=1)
            call VertNeighbor_NE(cd_tp, element, ev=2, ee=10, i1=1)
            call VertNeighbor_SI(cd_tp, ev=3, i1= 1, i2=2, i3=1)
            call VertNeighbor_NF(cd_tp, element, ev=4, ef= 2, i1=2, i2=1)
            call VertNeighbor_NE(cd_tp, element, ev=5, ee= 3, i1=1)
            call VertNeighbor_NV(cd_tp, element, ev=6)
            call VertNeighbor_NF(cd_tp, element, ev=7, ef= 6, i1=1, i2=2)
            call VertNeighbor_NE(cd_tp, element, ev=8, ee= 8, i1=2)

            ! geometry
            cd_tp % element(ce) % geometry % po = po
            cd_tp % start_point(ce) = cp + 1
            call InterpolateCoords( po = po                         &
                                  , A1  = iop(po) % A(:,:,2)        &
                                  , A2  = iop(po) % A(:,:,1)        &
                                  , A3  = iop(po) % A(:,:,2)        &
                                  , xp  = element % geometry % x_e  &
                                  , xc1 = cd_tp % x_e(cp+1:cp+np,1) &
                                  , xc2 = cd_tp % x_e(cp+1:cp+np,2) &
                                  , xc3 = cd_tp % x_e(cp+1:cp+np,3) )
            cp = cp + np

          end if

          ! child 1,2,2 ........................................................

          if (id_child(1,2,2,e) > 0) then

            ce = ce + 1

            ! element data
            cd_tp % element(ce) % id = id_child(1,2,2,e)
            cd_tp % element(ce) % cluster_oct = 7
            cd_tp % element(ce) % frozen = element % adaptation % mark < 100
            cd_tp % element(ce) % adaptation % parent_proc = parent % proc
            cd_tp % element(ce) % adaptation % parent_id   = e
            cd_tp % element(ce) % adaptation % mark        = child_mark

            ! face data
            cd_tp % element(ce) % face(1) % boundary = element % face(1) % boundary
            cd_tp % element(ce) % face(4) % boundary = element % face(4) % boundary
            cd_tp % element(ce) % face(6) % boundary = element % face(6) % boundary

            cd_tp % start_neighbor(ce) = cn + 1

            nn = 0

            ! neighbors at faces 1:6
            call FaceNeighbor_NF(cd_tp, element, ef=1, i1=2, i2=2)
            call FaceNeighbor_SI(cd_tp, ef=2, i1=2, i2=2, i3=2)
            call FaceNeighbor_SI(cd_tp, ef=3, i1=1, i2=1, i3=2)
            call FaceNeighbor_NF(cd_tp, element, ef=4, i1=1, i2=2)
            call FaceNeighbor_SI(cd_tp, ef=5, i1=1, i2=2, i3=1)
            call FaceNeighbor_NF(cd_tp, element, ef=6, i1=1, i2=2)

            ! neighbors at edges 1:12
            call EdgeNeighbor_SI(cd_tp, ee= 1, i1=1, i2=1, i3=1)
            call EdgeNeighbor_NF(cd_tp, element, ee= 2, ef=4, i1=1, i2=1)
            call EdgeNeighbor_NF(cd_tp, element, ee= 3, ef=6, i1=1, i2=1)
            call EdgeNeighbor_NE(cd_tp, element, ee= 4, i1=1)
            call EdgeNeighbor_NF(cd_tp, element, ee= 5, ef=1, i1=2, i2=1)
            call EdgeNeighbor_SI(cd_tp, ee= 6, i1=2, i2=2, i3=1)
            call EdgeNeighbor_NE(cd_tp, element, ee= 7, i1=2)
            call EdgeNeighbor_NF(cd_tp, element, ee= 8, ef=6, i1=2, i2=2)
            call EdgeNeighbor_NF(cd_tp, element, ee= 9, ef=1, i1=1, i2=2)
            call EdgeNeighbor_SI(cd_tp, ee=10, i1=2, i2=1, i3=2)
            call EdgeNeighbor_NE(cd_tp, element, ee=11, i1=2)
            call EdgeNeighbor_NF(cd_tp, element, ee=12, ef=4, i1=2, i2=2)

            ! neighbors at vertices 1:8
            call VertNeighbor_NF(cd_tp, element, ev=1, ef= 1, i1=1, i2=1)
            call VertNeighbor_SI(cd_tp, ev=2, i1= 2, i2=1, i3=1)
            call VertNeighbor_NE(cd_tp, element, ev=3, ee=11, i1=1)
            call VertNeighbor_NF(cd_tp, element, ev=4, ef= 4, i1=2, i2=1)
            call VertNeighbor_NE(cd_tp, element, ev=5, ee= 7, i1=1)
            call VertNeighbor_NF(cd_tp, element, ev=6, ef= 6, i1=2, i2=1)
            call VertNeighbor_NV(cd_tp, element, ev=7)
            call VertNeighbor_NE(cd_tp, element, ev=8, ee= 4, i1=2)

            ! geometry
            cd_tp % element(ce) % geometry % po = po
            cd_tp % start_point(ce) = cp + 1
            call InterpolateCoords( po = po                         &
                                  , A1  = iop(po) % A(:,:,1)        &
                                  , A2  = iop(po) % A(:,:,2)        &
                                  , A3  = iop(po) % A(:,:,2)        &
                                  , xp  = element % geometry % x_e  &
                                  , xc1 = cd_tp % x_e(cp+1:cp+np,1) &
                                  , xc2 = cd_tp % x_e(cp+1:cp+np,2) &
                                  , xc3 = cd_tp % x_e(cp+1:cp+np,3) )
            cp = cp + np

          end if

          ! child 2,2,2 ........................................................

          if (id_child(2,2,2,e) > 0) then

            ce = ce + 1

            ! element data
            cd_tp % element(ce) % id = id_child(2,2,2,e)
            cd_tp % element(ce) % cluster_oct = 8
            cd_tp % element(ce) % frozen = element % adaptation % mark < 100
            cd_tp % element(ce) % adaptation % parent_proc = parent % proc
            cd_tp % element(ce) % adaptation % parent_id   = e
            cd_tp % element(ce) % adaptation % mark        = child_mark

            ! face data
            cd_tp % element(ce) % face(2) % boundary = element % face(2) % boundary
            cd_tp % element(ce) % face(4) % boundary = element % face(4) % boundary
            cd_tp % element(ce) % face(6) % boundary = element % face(6) % boundary

            cd_tp % start_neighbor(ce) = cn + 1

            nn = 0

            ! neighbors at faces 1:6
            call FaceNeighbor_SI(cd_tp, ef=1, i1=1, i2=2, i3=2)
            call FaceNeighbor_NF(cd_tp, element, ef=2, i1=2, i2=2)
            call FaceNeighbor_SI(cd_tp, ef=3, i1=2, i2=1, i3=2)
            call FaceNeighbor_NF(cd_tp, element, ef=4, i1=2, i2=2)
            call FaceNeighbor_SI(cd_tp, ef=5, i1=2, i2=2, i3=1)
            call FaceNeighbor_NF(cd_tp, element, ef=6, i1=2, i2=2)

            ! neighbors at edges 1:12
            call EdgeNeighbor_SI(cd_tp, ee= 1, i1=2, i2=1, i3=1)
            call EdgeNeighbor_NF(cd_tp, element, ee= 2, ef=4, i1=2, i2=1)
            call EdgeNeighbor_NF(cd_tp, element, ee= 3, ef=6, i1=2, i2=1)
            call EdgeNeighbor_NE(cd_tp, element, ee= 4, i1=2)
            call EdgeNeighbor_SI(cd_tp, ee= 5, i1=1, i2=2, i3=1)
            call EdgeNeighbor_NF(cd_tp, element, ee= 6, ef=2, i1=2, i2=1)
            call EdgeNeighbor_NF(cd_tp, element, ee= 7, ef=6, i1=1, i2=2)
            call EdgeNeighbor_NE(cd_tp, element, ee= 8, i1=2)
            call EdgeNeighbor_SI(cd_tp, ee= 9, i1=1, i2=1, i3=2)
            call EdgeNeighbor_NF(cd_tp, element, ee=10, ef=2, i1=1, i2=2)
            call EdgeNeighbor_NF(cd_tp, element, ee=11, ef=4, i1=1, i2=2)
            call EdgeNeighbor_NE(cd_tp, element, ee=12, i1=2)

            ! neighbors at vertices 1:8
            call VertNeighbor_SI(cd_tp, ev=1, i1= 1, i2=1, i3=1)
            call VertNeighbor_NF(cd_tp, element, ev=2, ef= 2, i1=1, i2=1)
            call VertNeighbor_NF(cd_tp, element, ev=3, ef= 4, i1=1, i2=1)
            call VertNeighbor_NE(cd_tp, element, ev=4, ee=12, i1=1)
            call VertNeighbor_NF(cd_tp, element, ev=5, ef= 6, i1=1, i2=1)
            call VertNeighbor_NE(cd_tp, element, ev=6, ee= 8, i1=1)
            call VertNeighbor_NE(cd_tp, element, ev=7, ee= 4, i1=1)
            call VertNeighbor_NV(cd_tp, element, ev=8)

            ! geometry
            cd_tp % element(ce) % geometry % po = po
            cd_tp % start_point(ce) = cp + 1
            call InterpolateCoords( po = po                         &
                                  , A1  = iop(po) % A(:,:,2)        &
                                  , A2  = iop(po) % A(:,:,2)        &
                                  , A3  = iop(po) % A(:,:,2)        &
                                  , xp  = element % geometry % x_e  &
                                  , xc1 = cd_tp % x_e(cp+1:cp+np,1) &
                                  , xc2 = cd_tp % x_e(cp+1:cp+np,2) &
                                  , xc3 = cd_tp % x_e(cp+1:cp+np,3) )
            cp = cp + np

          end if

        end associate
      end do
!### CHECK
write(*,'(99(G0,1X))') '# BCD  4',', parent%proc',parent%proc
!### CHECK END

      !-------------------------------------------------------------------------
      ! finalization

      !$omp master
      deallocate( start_child_element, start_child_neighbor, start_child_point )
      deallocate( ce_part, cn_part, cp_part )
      deallocate( id_child_face, tp_child_face )
      deallocate( id_child_edge, tp_child_edge )
      deallocate( id_child_vert, tp_child_vert )
      deallocate( eop, iop )
      !$omp end master
      !$omp end barrier

    end associate
!### CHECK
write(*,'(99(G0,1X))') '# BCD  X',', parent%proc',parent%proc
!### CHECK END

  contains

    !---------------------------------------------------------------------------
    !> Face neighbor at parent neighbor face

    subroutine FaceNeighbor_NF(cd_tp, element, ef, i1, i2)
      class(ElementDistributionData_3D), intent(inout) :: cd_tp
      class(MeshElement_3D), intent(in) :: element
      integer, intent(in) :: ef     !< element face
      integer, intent(in) :: i1, i2 !< neighbor face child double index

      integer :: ii

      associate(nb_id => id_child_face(i1,i2,ef), nb_tp => tp_child_face(ef))

        if (nb_id <= 0) return

        cn = cn + 1
        nn = nn + 1

        cd_tp % element(ce) % face(ef) % n_neighbor = int( 1 , IXS )
        cd_tp % element(ce) % face(ef) % i_neighbor = int( nn, IXS )

        ii = element % face(ef) % i_neighbor

        cd_tp % neighbor_id          (cn) = nb_id
        cd_tp % neighbor_part        (cn) = nb_tp
        cd_tp % neighbor_component   (cn) = element % neighbor(ii) % component
        cd_tp % neighbor_orientation (cn) = element % neighbor(ii) % orientation

      end associate

    end subroutine FaceNeighbor_NF

    !---------------------------------------------------------------------------
    !> Sibling face neighbor

    subroutine FaceNeighbor_SI(cd_tp, ef, i1, i2, i3)
      class(ElementDistributionData_3D), intent(inout) :: cd_tp
      integer, intent(in) :: ef         !< element face
      integer, intent(in) :: i1, i2, i3 !< element child triple index

      associate(nb_id => map%id_child(i1,i2,i3,e), nb_tp => tp)

        if (nb_id <= 0) return

        cn = cn + 1
        nn = nn + 1

        cd_tp % element(ce) % face(ef) % n_neighbor = int( 1 , IXS )
        cd_tp % element(ce) % face(ef) % i_neighbor = int( nn, IXS )

        cd_tp % neighbor_id   (cn) = nb_id
        cd_tp % neighbor_part (cn) = nb_tp

        select case(ef)
        case(1)
          cd_tp % neighbor_component (cn) = 2
        case(2)
          cd_tp % neighbor_component (cn) = 1
        case(3)
          cd_tp % neighbor_component (cn) = 4
        case(4)
          cd_tp % neighbor_component (cn) = 3
        case(5)
          cd_tp % neighbor_component (cn) = 6
        case(6)
          cd_tp % neighbor_component (cn) = 5
        end select

        cd_tp % neighbor_orientation (cn) = 12  ! always aligned !

      end associate

    end subroutine FaceNeighbor_SI

    !---------------------------------------------------------------------------
    !> Edge neighbors at parent neighbor edge

    subroutine EdgeNeighbor_NE(cd_tp, element, ee, i1)
      class(ElementDistributionData_3D), intent(inout) :: cd_tp
      class(MeshElement_3D), intent(in) :: element
      integer, intent(in) :: ee !< element edge
      integer, intent(in) :: i1 !< neighbor edge child index

      integer :: ii, jj, nn_ce

      associate(nb_id => id_child_edge(i1,:,ee), nb_tp => tp_child_edge(:,ee))

        nn_ce = count(nb_id > 0)
        if (nn_ce == 0) return

        cd_tp % element(ce) % edge(ee) % n_neighbor = int( nn_ce , IXS )
        cd_tp % element(ce) % edge(ee) % i_neighbor = int( nn + 1, IXS )

        ii = element % edge(ee) % i_neighbor
        do jj = 1, element % edge(ee) % n_neighbor
          if (nb_id(jj) > 0) then
            cn = cn + 1
            cd_tp % neighbor_id          (cn) = nb_id(jj)
            cd_tp % neighbor_part        (cn) = nb_tp(jj)
            cd_tp % neighbor_component   (cn) = element % neighbor(ii) % component
            cd_tp % neighbor_orientation (cn) = element % neighbor(ii) % orientation
          end if
          ii = ii + 1
        end do
        nn = nn + nn_ce

      end associate

    end subroutine EdgeNeighbor_NE

    !---------------------------------------------------------------------------
    !> Edge neighbor at parent neighbor face

    subroutine EdgeNeighbor_NF(cd_tp, element, ee, ef, i1, i2)
      class(ElementDistributionData_3D), intent(inout) :: cd_tp
      class(MeshElement_3D), intent(in) :: element
      integer, intent(in) :: ee     !< element edge
      integer, intent(in) :: ef     !< element face
      integer, intent(in) :: i1, i2 !< neighbor face child double index

      integer(IXS) :: component, orientation
      integer :: xi_o(3), xi_t(3)
      integer :: ii, nee

      associate(nb_id => id_child_face(i1,i2,ef), nb_tp => tp_child_face(ef))

        if (nb_id <= 0) return

        cn = cn + 1
        nn = nn + 1

        cd_tp % element(ce) % edge(ee) % n_neighbor = int( 1 , IXS )
        cd_tp % element(ce) % edge(ee) % i_neighbor = int( nn, IXS )

        ! adopt orientation from parent face neighbor ...........................

        ii = element % face(ef) % i_neighbor
        orientation = element % neighbor(ii) % orientation

        ! identify matching adjacent child edge .................................

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
        if (xi_t(1) == 0) then
          nee = 1 + (xi_t(2) + 1)/2 + (xi_t(3) + 1)
        else if (xi_t(2) == 0) then
          nee = 5 + (xi_t(1) + 1)/2 + (xi_t(3) + 1)
        else if (xi_t(3) == 0) then
          nee = 9 + (xi_t(1) + 1)/2 + (xi_t(2) + 1)
        else ! should not happen, trigger error by providing invalid number
          nee = -100
        end if

        ! 4: neighbor component ID
        component = int(6 + nee, IXS)

        ! assign neighbor data .................................................

        cd_tp % neighbor_id          (cn) = nb_id
        cd_tp % neighbor_part        (cn) = nb_tp
        cd_tp % neighbor_component   (cn) = component
        cd_tp % neighbor_orientation (cn) = orientation

      end associate

    end subroutine EdgeNeighbor_NF

    !---------------------------------------------------------------------------
    !> Sibling edge neighbor

    subroutine EdgeNeighbor_SI(cd_tp, ee, i1, i2, i3)
      class(ElementDistributionData_3D), intent(inout) :: cd_tp
      integer, intent(in) :: ee         !< element edge
      integer, intent(in) :: i1, i2, i3 !< element child triple index

      associate(nb_id => map%id_child(i1,i2,i3,e), nb_tp => tp)

        if (nb_id <= 0) return

        cn = cn + 1
        nn = nn + 1

        cd_tp % element(ce) % edge(ee) % n_neighbor = int( 1 , IXS )
        cd_tp % element(ce) % edge(ee) % i_neighbor = int( nn, IXS )

        cd_tp % neighbor_id   (cn) = nb_id
        cd_tp % neighbor_part (cn) = nb_tp

        ! ID of coupled component = ID of opposite element edge
        select case(ee)
        case(1:4)
          cd_tp % neighbor_component (cn) = int(6 +  5 - ee, IXS)
        case(5:8)
          cd_tp % neighbor_component (cn) = int(6 + 13 - ee, IXS)
        case(9:12)
          cd_tp % neighbor_component (cn) = int(6 + 21 - ee, IXS)
        end select

        cd_tp % neighbor_orientation (cn) = 12  ! always aligned !

      end associate

    end subroutine EdgeNeighbor_SI

    !---------------------------------------------------------------------------
    !> Vertex neighbors at parent neighbor vertex

    subroutine VertNeighbor_NV(cd_tp, element, ev)
      class(ElementDistributionData_3D), intent(inout) :: cd_tp
      class(MeshElement_3D), intent(in) :: element
      integer, intent(in) :: ev !< element vertex

      integer :: ii, jj, nn_ce

      associate(nb_id => id_child_vert(:,ev), nb_tp => tp_child_vert(:,ev))

        nn_ce = count(nb_id > 0)
        if (nn_ce== 0) return

        cd_tp % element(ce) % vertex(ev) % n_neighbor = int( nn_ce , IXS )
        cd_tp % element(ce) % vertex(ev) % i_neighbor = int( nn + 1, IXS )

        ii = element % vertex(ev) % i_neighbor
        do jj = 1, element % vertex(ev) % n_neighbor
          if (nb_id(jj) > 0) then
            cn = cn + 1
            cd_tp % neighbor_id          (cn) = nb_id(jj)
            cd_tp % neighbor_part        (cn) = nb_tp(jj)
            cd_tp % neighbor_component   (cn) = element % neighbor(ii) % component
            cd_tp % neighbor_orientation (cn) = element % neighbor(ii) % orientation
          end if
          ii = ii + 1
        end do
        nn = nn + nn_ce

      end associate

    end subroutine VertNeighbor_NV

    !---------------------------------------------------------------------------
    !> Vertex neighbors at parent neighbor edge

    subroutine VertNeighbor_NE(cd_tp, element, ev, ee, i1)
      class(ElementDistributionData_3D), intent(inout) :: cd_tp
      class(MeshElement_3D), intent(in) :: element
      integer, intent(in) :: ev !< element vertex
      integer, intent(in) :: ee !< element edge
      integer, intent(in) :: i1 !< neighbor edge child index

      integer :: ii, jj, nee, nev, nn_ce

      associate(nb_id => id_child_edge(i1,:,ee), nb_tp => tp_child_edge(:,ee))

        nn_ce = count(nb_id > 0)
        if (nn_ce == 0) return

        cd_tp % element(ce) % vertex(ev) % n_neighbor = int( nn_ce , IXS )
        cd_tp % element(ce) % vertex(ev) % i_neighbor = int( nn + 1, IXS )

        ii = element % edge(ee) % i_neighbor
        do jj = 1, element % edge(ee) % n_neighbor
          if (nb_id(jj) > 0) then
            cn = cn + 1

            ! parent neighbor edge
            nee = ElementEdgeID(element % neighbor(ii) % component)

            ! matching neighbor vertex
            if (element % NeigborEdgeIsAligned(ee, ii)) then
              nev = V_EDGE(3 - i1, nee)
            else
              nev = V_EDGE(i1, nee)
            end if

            cd_tp % neighbor_id          (cn) = nb_id(jj)
            cd_tp % neighbor_part        (cn) = nb_tp(jj)
            cd_tp % neighbor_component   (cn) = int(18 + nev, IXS)
            cd_tp % neighbor_orientation (cn) = element % neighbor(ii) % orientation

          end if
          ii = ii + 1
        end do
        nn = nn + nn_ce

      end associate

    end subroutine VertNeighbor_NE

    !---------------------------------------------------------------------------
    !> Vertex neighbor at parent neighbor face

    subroutine VertNeighbor_NF(cd_tp, element, ev, ef, i1, i2)
      class(ElementDistributionData_3D), intent(inout) :: cd_tp
      class(MeshElement_3D), intent(in) :: element
      integer, intent(in) :: ev     !< element vertex
      integer, intent(in) :: ef     !< element face
      integer, intent(in) :: i1, i2 !< neighbor face child double index

      integer(IXS) :: orientation
      integer :: xi_o(3), xi_t(3)
      integer :: ii, nev

      associate(nb_id => id_child_face(i1,i2,ef), nb_tp => tp_child_face(ef))

        if (nb_id <= 0) return

        cn = cn + 1
        nn = nn + 1

        cd_tp % element(ce) % vertex(ev) % n_neighbor = int( 1 , IXS )
        cd_tp % element(ce) % vertex(ev) % i_neighbor = int( nn, IXS )

        ! adopt orientation from parent face neighbor ..........................

        ii = element % face(ef) % i_neighbor
        orientation = element % neighbor(ii) % orientation

        ! identify matching adjacent child vertex ..............................

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
        nev = 1 + (xi_t(1) + 1)/2 + (xi_t(2) + 1) + (xi_t(3) + 1)*2

        ! assign neighbor data .................................................

        cd_tp % neighbor_id          (cn) = nb_id
        cd_tp % neighbor_part        (cn) = nb_tp
        cd_tp % neighbor_component   (cn) = int(18 + nev, IXS)
        cd_tp % neighbor_orientation (cn) = orientation

      end associate

    end subroutine VertNeighbor_NF

    !---------------------------------------------------------------------------
    !> Sibling vertex neighbor

    subroutine VertNeighbor_SI(cd_tp, ev, i1, i2, i3)
      class(ElementDistributionData_3D), intent(inout) :: cd_tp
      integer, intent(in) :: ev         !< element vertex
      integer, intent(in) :: i1, i2, i3 !< element child triple index

      associate(nb_id => map%id_child(i1,i2,i3,e), nb_tp => tp)

        if (nb_id <= 0) return

        cn = cn + 1
        nn = nn + 1

        cd_tp % element(ce) % vertex(ev) % n_neighbor = int( 1 , IXS)
        cd_tp % element(ce) % vertex(ev) % i_neighbor = int( nn, IXS)

        cd_tp % neighbor_id          (cn) = nb_id
        cd_tp % neighbor_part        (cn) = nb_tp
        cd_tp % neighbor_component   (cn) = int(18 + 9 - ev, IXS)
        cd_tp % neighbor_orientation (cn) = 12
        ! Note
        ! – component:   is the diagonally opposite vertex
        ! – orientation: sibling neighbor is always aligned

      end associate

    end subroutine VertNeighbor_SI

    !---------------------------------------------------------------------------

  end subroutine BuildChildData

  !-----------------------------------------------------------------------------
  !> Interpolates parent element points to child

  subroutine InterpolateCoords(po, A1, A2, A3, xp, xc1, xc2, xc3)
    integer,   intent(in)  :: po               !< polynomial order
    real(RNP), intent(in)  :: A1  (0:po, 0:po) !< x₁ interpolation operator
    real(RNP), intent(in)  :: A2  (0:po, 0:po) !< x₂ interpolation operator
    real(RNP), intent(in)  :: A3  (0:po, 0:po) !< x₃ interpolation operator
    real(RNP), intent(in)  :: xp  (0:po, 0:po, 0:po, 3) !< parent coordinates
    real(RNP), intent(out) :: xc1 (0:po, 0:po, 0:po)    !< child x₁ coordinates
    real(RNP), intent(out) :: xc2 (0:po, 0:po, 0:po)    !< child x₂ coordinates
    real(RNP), intent(out) :: xc3 (0:po, 0:po, 0:po)    !< child x₃ coordinates

    real(RNP) :: z1(0:po, 0:po, 0:po, 3)
    real(RNP) :: z2(0:po, 0:po, 0:po, 3)

    integer :: i, j, k, p

    do k = 0, po
    do j = 0, po
    do i = 0, po
      z1(i,j,k,1) = 0
      z1(i,j,k,2) = 0
      z1(i,j,k,3) = 0
      do p = 0, po
        z1(i,j,k,1) = z1(i,j,k,1) + A1(i,p) * xp(p,j,k,1)
        z1(i,j,k,2) = z1(i,j,k,2) + A1(i,p) * xp(p,j,k,2)
        z1(i,j,k,3) = z1(i,j,k,3) + A1(i,p) * xp(p,j,k,3)
      end do
    end do
    end do
    end do

    do k = 0, po
    do j = 0, po
    do i = 0, po
      z2(i,j,k,1) = 0
      z2(i,j,k,2) = 0
      z2(i,j,k,3) = 0
      do p = 0, po
        z2(i,j,k,1) = z2(i,j,k,1) + A2(j,p) * z1(i,p,k,1)
        z2(i,j,k,2) = z2(i,j,k,2) + A2(j,p) * z1(i,p,k,2)
        z2(i,j,k,3) = z2(i,j,k,3) + A2(j,p) * z1(i,p,k,3)
      end do
    end do
    end do
    end do

    do k = 0, po
    do j = 0, po
    do i = 0, po
      xc1(i,j,k) = 0
      xc2(i,j,k) = 0
      xc3(i,j,k) = 0
      do p = 0, po
        xc1(i,j,k) = xc1(i,j,k) + A3(k,p) * z2(i,j,p,1)
        xc2(i,j,k) = xc2(i,j,k) + A3(k,p) * z2(i,j,p,2)
        xc3(i,j,k) = xc3(i,j,k) + A3(k,p) * z2(i,j,p,3)
      end do
    end do
    end do
    end do

  end subroutine InterpolateCoords

  !=============================================================================

end submodule MP_BuildChildData
