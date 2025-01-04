module Root_Mesh_Partitioning__3D
  use Logging_Levels
  use XMPI
  use Mesh__3D
  use Element_Transfer_Buffer__3D
  use Element_Distribution_Map__3D
  use Element_Distribution_Data__3D
  use Partitioner_Interface__3D
  use Data_Exchange__3D

  implicit none
  private

  public :: PartitioningOptions_3D   ! from Partitioner_Interface__3D
  public :: RootMeshPartitioning_3D

contains

  !=============================================================================
  ! Partitioning routines

  !-----------------------------------------------------------------------------
  !> Redistribution of root mesh data

  subroutine RootMeshPartitioning_3D(opt, old_mesh, new_mesh, child, x_plan)

    ! arguments ................................................................

    class(PartitioningOptions_3D),        intent(in)    :: opt
    class(Mesh_3D),                       intent(in)    :: old_mesh
    class(Mesh_3D),                       intent(out)   :: new_mesh
    class(Mesh_3D),             optional, intent(inout) :: child
    class(DataExchangePlan_3D), optional, intent(out)   :: x_plan

    ! internal variables .......................................................

    type(ElementDistributionMap_3D) :: dist_map
    type(ElementDistributionData_3D), allocatable :: dist_send_data(:)
    type(ElementDistributionData_3D), allocatable :: dist_recv_data(:)

    type(MeshAttributes_3D) :: attrib
    character(len=:), allocatable :: prefix
    logical :: logging
    integer :: i

    !---------------------------------------------------------------------------
    !$omp master

    ! initialize logging .......................................................

    logging = old_mesh%proc == 0 .and. log_level > 0 .or. &
              old_mesh%proc  > 0 .and. log_level > 1
    prefix  = LoggingPrefix('RootMeshPartitioning_3D', old_mesh%proc)

    if (logging) then
      print '(2A)', prefix, 'start'
    end if

    ! prepare new partitions ...................................................

    attrib = MeshAttributes_3D(old_mesh)  ! adopt attributes from old mesh
    attrib % structured = .false.         ! new mesh is always unstructured

    call new_mesh % Init_Mesh_3D(attrib, old_mesh % comm_world)

    new_mesh % n_parts = opt % n_parts

    ! use identity for mapping parts to processes, so far
    if (new_mesh % proc < new_mesh % n_parts) then
      new_mesh % part = new_mesh % proc
    end if

    call new_mesh % BuildCommunicator()

    ! build new element distribution map .......................................

    if (old_mesh % part >= 0) then
      call BuildElementDistributionMap(opt, old_mesh, dist_map)
    end if

    ! distribute element data ..................................................

    ! prepare send data
    if (old_mesh % part >= 0) then
      allocate(dist_send_data(0:new_mesh%n_parts-1))
      call BuildElementDistributionData(old_mesh, new_mesh, dist_map, dist_send_data)
    end if

    ! prepare receive data
    if (new_mesh % part >= 0) then
      allocate(dist_recv_data(0:old_mesh%n_parts-1))
      do i = 0, old_mesh%n_parts-1
        dist_recv_data(i) % proc = old_mesh % proc_part(i)
        dist_recv_data(i) % comm = old_mesh % comm_world
      end do
    end if

    ! send dimensions (skip local part)
    if (old_mesh % part >= 0) then
      do i = 0, new_mesh%n_parts-1
        if (dist_send_data(i) % proc == old_mesh % proc) cycle
        call dist_send_data(i) % Send_Dimensions()
      end do
    end if

    ! receive dimensions (skip local part)
    if (new_mesh % part >= 0) then
      do i = 0, old_mesh%n_parts-1
        if (dist_recv_data(i) % proc == old_mesh % proc) cycle
        call dist_recv_data(i) % Recv_Dimensions()
      end do
    end if

    ! send data (skip local part)
    if (old_mesh % part >= 0) then
      do i = 0, new_mesh%n_parts-1
        if (dist_send_data(i) % proc == old_mesh % proc) cycle
        call dist_send_data(i) % Send_Data()
      end do
    end if

    ! receive data (skip local part)
    if (new_mesh % part >= 0) then
      do i = 0, old_mesh%n_parts-1
        if (dist_recv_data(i) % proc == old_mesh % proc) cycle
        call dist_recv_data(i) % Recv_Data()
      end do
    end if

    ! finalize transfer
    if (old_mesh % part >= 0) then
      do i = 0, new_mesh%n_parts-1
        if (dist_send_data(i) % proc == old_mesh % proc) cycle
        call dist_send_data(i) % Send_Finish()
      end do
    end if
    if (new_mesh % part >= 0) then
      do i = 0, old_mesh%n_parts-1
        if (dist_recv_data(i) % proc == old_mesh % proc) cycle
        call dist_recv_data(i) % Recv_Finish()
      end do
    end if

    ! assign element data ......................................................

    if (new_mesh % part >= 0) then

      ! number of new elements
      if (old_mesh % part >= 0) then
        new_mesh % n_elem = sum(dist_recv_data % n_elem) &
                          + dist_send_data(old_mesh % part) % n_elem
      else
        new_mesh % n_elem = sum(dist_recv_data % n_elem)
      end if
      new_mesh % n_elem_active = new_mesh % n_elem
      new_mesh % n_elem_frozen = 0

      allocate(new_mesh % element(new_mesh % n_elem))

      ! assign received data
      if (new_mesh % part >= 0) then
        do i = 0, old_mesh%n_parts-1
          if (dist_recv_data(i) % proc /= old_mesh % proc) then
            call dist_recv_data(i) % Assign_Data(new_mesh)
          else
            call dist_send_data(old_mesh%part) % Assign_Data(new_mesh)
          end if
        end do
      end if

    end if

    ! complete mesh.............................................................

    if (new_mesh % part >= 0) then

      call new_mesh % IdentifyVertices()
      call new_mesh % IdentifyEdges()
      call new_mesh % BuildFaces()
      call new_mesh % BuildBoundaryFaces()
      call new_mesh % BuildLinks()
      call new_mesh % BuildGhosts()
      call new_mesh % BuildCuboids()
      call new_mesh % IdentifyRanks()
      call new_mesh % BuildMapToChild()

    end if

    ! redistribution maps for retained data ....................................

    if (present(x_plan)) then

      if (new_mesh % part >= 0) then

        allocate(x_plan % send_map( 0 : new_mesh % n_parts - 1 ))
        allocate(x_plan % recv_map( 0 : old_mesh % n_parts - 1 ))

        do i = 0, new_mesh % n_parts - 1
          x_plan % send_map(i) % comm = new_mesh % comm_world
          x_plan % send_map(i) % proc = new_mesh % proc_part(i)
        end do
        call BuildDataExchangeSendMap(dist_map, x_plan % send_map)

        call BuildDataExchangeRecvMap(dist_recv_data, x_plan % recv_map)

      else

        allocate(x_plan % send_map(0))
        allocate(x_plan % recv_map(0))

      end if
    end if

    ! update child partitions ..................................................

    if (present(child)) then
      call UpdateChild(old_mesh, new_mesh, dist_map, child)
    end if

    if (logging) then
      print '(2A)', prefix, 'exit'
    end if

    !---------------------------------------------------------------------------
    !$omp end master

  end subroutine RootMeshPartitioning_3D

  !-----------------------------------------------------------------------------
  !> Computes the new target partitions and element IDs

  subroutine BuildElementDistributionMap(opt, old_mesh, dist_map)
    class(PartitioningOptions_3D),   intent(in)  :: opt
    class(Mesh_3D),                  intent(in)  :: old_mesh
    type(ElementDistributionMap_3D), intent(out) :: dist_map

    integer, allocatable, target :: tp_elem(:)
    integer, contiguous, pointer :: tp_elem_val(:,:,:,:)
    type(ElementTransferBuffer_3D), allocatable, asynchronous :: tp_elem_buf

    allocate(tp_elem(old_mesh%n_elem + old_mesh%n_ghost), source = -1)

    ! graph-based partitioning
    call ParMETIS_Partitioner_3D(opt, old_mesh, tp_elem)

    ! transfer target partition IDs to ghosts
    if (old_mesh % n_ghost > 0) then
      tp_elem_val(1:1,1:1,1:1,1:old_mesh%n_elem+old_mesh%n_ghost) => tp_elem
      tp_elem_buf = ElementTransferBuffer_3D(old_mesh, tp_elem_val)
      call tp_elem_buf % Transfer(old_mesh, tp_elem_val, tag = 1000)
      call tp_elem_buf % Merge(tp_elem_val)
    end if

    ! new distribution map, assuming identity for proc_part
    dist_map = ElementDistributionMap_3D(old_mesh, opt%n_parts, tp_elem)

  end subroutine BuildElementDistributionMap

  !=============================================================================
  ! Distribution routines

  !-----------------------------------------------------------------------------
  !> Build element distribution data for sending from old to new partitions

  subroutine BuildElementDistributionData(old_mesh, new_mesh, dist_map, dist_data)

    class(Mesh_3D), intent(in) :: old_mesh
    class(Mesh_3D), intent(in) :: new_mesh
    class(ElementDistributionMap_3D) , intent(in)  :: dist_map
    class(ElementDistributionData_3D), intent(out) :: dist_data(0:)

    ! internal variables .......................................................

    integer, allocatable, save :: rk_elem(:) ! rank among local elements in tp
    integer, allocatable, save :: ce_part(:) ! element counters
    integer, allocatable, save :: cn_part(:) ! neighbor counters
    integer, allocatable, save :: cp_part(:) ! element point counters

    integer :: ce, cn, cp, ne, nn, np, po, tp
    integer :: i, j, k, l

    if (old_mesh % part < 0) return

    !---------------------------------------------------------------------------
    ! initialization

    !!$omp master

    allocate(rk_elem( old_mesh%n_elem ))
    allocate(ce_part( 0:dist_map%n_parts-1 ))
    allocate(cn_part( 0:dist_map%n_parts-1 ))
    allocate(cp_part( 0:dist_map%n_parts-1 ))

    ! basic initialization .....................................................

    do tp = 0, dist_map%n_parts-1

      ne = dist_map % ne_part(tp)

      dist_data(tp) % proc       = new_mesh % proc_part(tp)
      dist_data(tp) % comm       = old_mesh % comm_world
      dist_data(tp) % n_elem  = ne
      dist_data(tp) % n_neighbor = 0
      dist_data(tp) % n_point    = 0

      allocate(dist_data(tp) % element        (ne) )
      allocate(dist_data(tp) % start_neighbor (ne) )
      allocate(dist_data(tp) % start_point    (ne) )

    end do

    ! counts and start indices .................................................

    ce_part = 0
    cn_part = 0
    cp_part = 0

    do l = 1, old_mesh % n_elem

      tp = dist_map % tp_elem(l)
      ce = ce_part(tp) + 1
      nn = size(old_mesh % element(l) % neighbor)
      np = (old_mesh % element(l) % geometry % po + 1) ** 3

      rk_elem(l) = ce

      dist_data(tp) % start_neighbor (ce) = cn_part(tp) + 1
      dist_data(tp) % start_point    (ce) = cp_part(tp) + 1

      ce_part(tp) = ce_part(tp) + 1
      cn_part(tp) = cn_part(tp) + nn
      cp_part(tp) = cp_part(tp) + np

    end do

    ! remaining components .....................................................

    do tp = 0, dist_map%n_parts-1

      nn = cn_part(tp)
      np = cp_part(tp)

      dist_data(tp) % n_neighbor  =  nn
      dist_data(tp) % n_point     =  np

      allocate(dist_data(tp) % neighbor_id          ( nn    ) )
      allocate(dist_data(tp) % neighbor_part        ( nn    ) )
      allocate(dist_data(tp) % neighbor_component   ( nn    ) )
      allocate(dist_data(tp) % neighbor_orientation ( nn    ) )
      allocate(dist_data(tp) % x_e                  ( np, 3 ) )

    end do

    !!$omp end master
    !!$omp end barrier
    !---------------------------------------------------------------------------

    ! build components .........................................................

    ! extract and transform element data
    !!$omp do
    do l = 1, old_mesh % n_elem
      tp = dist_map % tp_elem(l)
      associate(dd => dist_data(tp))

        ! static element components
        associate(element => old_mesh % element(l))
          ce = rk_elem(l)
          dd % element(ce) % id              = dist_map % id_elem(l)
          dd % element(ce) % vertex          = element % vertex
          dd % element(ce) % edge            = element % edge
          dd % element(ce) % face            = element % face
          dd % element(ce) % adaptation      = element % adaptation
          dd % element(ce) % geometry % po   = element % geometry % po
          dd % element(ce) % geometry % x_c  = element % geometry % x_c
          dd % element(ce) % geometry % dx_m = element % geometry % dx_m
        end associate

        ! neighbor data
        associate(neighbor => old_mesh % element(l) % neighbor)
          cn = dd % start_neighbor(ce)
          nn = size(neighbor)
          do i = 1, nn
            dd % neighbor_id          (cn) = dist_map % id_elem( neighbor(i)%id )
            dd % neighbor_part        (cn) = dist_map % tp_elem( neighbor(i)%id )
            dd % neighbor_component   (cn) = neighbor(i) % component
            dd % neighbor_orientation (cn) = neighbor(i) % orientation
            cn = cn + 1
          end do
        end associate

        ! element points
        associate(x_e => old_mesh % element(l) % geometry % x_e)
          po = old_mesh % element(l) % geometry % po
          cp = dd % start_point(ce)
          do k = 0, po
          do j = 0, po
          do i = 0, po
            dd % x_e(cp,1) = x_e(i,j,k,1)
            dd % x_e(cp,2) = x_e(i,j,k,2)
            dd % x_e(cp,3) = x_e(i,j,k,3)
            cp = cp + 1
          end do
          end do
          end do
        end associate

      end associate
    end do

    !!$omp master
    deallocate(rk_elem, ce_part, cn_part, cp_part)
    !!$omp end master
    !!$omp end barrier

  end subroutine BuildElementDistributionData

  !=============================================================================
  ! Retained data redistribution routines

  !-----------------------------------------------------------------------------
  !> Build map indicating the new location of retained elements

  subroutine BuildDataExchangeSendMap(dist_map, send_map)
    class(ElementDistributionMap_3D), intent(in)    :: dist_map
    class(DataExchangeMap_3D),        intent(inout) :: send_map(0:)

    integer, allocatable, save :: e(:)
    integer :: i, p

    associate(n_parts => dist_map % n_parts)

      do i = 0, n_parts-1
        allocate(send_map(i) % id_elem( dist_map % ne_part(i)), source = -1)
      end do

      allocate(e(0:n_parts-1), source = 0)

      do i = 1, dist_map % n_elem
        p = dist_map % tp_elem(i)
        if (p < 0) cycle
        e(p) = e(p) + 1
        send_map(p) % id_elem(e(p)) = i
      end do

      deallocate(e)

    end associate

  end subroutine BuildDataExchangeSendMap

  !-----------------------------------------------------------------------------
  !> Build map indicating the old location of retained elements

  subroutine BuildDataExchangeRecvMap(dist_recv_data, recv_map)
    class(ElementDistributionData_3D), intent(in)  :: dist_recv_data(0:)
    class(DataExchangeMap_3D),         intent(out) :: recv_map(0:)

    integer :: i, e

    do i = 0, ubound(dist_recv_data,1)
      recv_map(i) % comm = dist_recv_data(i) % comm
      recv_map(i) % proc = dist_recv_data(i) % proc
      allocate(recv_map(i) % id_elem( dist_recv_data(i) % n_elem ))
      do e = 1, dist_recv_data(i) % n_elem
        recv_map(i) % id_elem(e) = dist_recv_data(i) % element(e) % id
      end do
    end do

  end subroutine BuildDataExchangeRecvMap

  !=============================================================================
  ! Child update

  !-----------------------------------------------------------------------------
  !> Update parent info and maps in child level partitions

  subroutine UpdateChild(old_mesh, new_mesh, dist_map, child)
    class(Mesh_3D), intent(in) :: old_mesh
    class(Mesh_3D), intent(in) :: new_mesh
    class(ElementDistributionMap_3D), intent(in) :: dist_map
    class(Mesh_3D), intent(inout) :: child

    type(DataExchangeMap_3D), allocatable :: send_map(:)
    type(DataExchangeMap_3D), allocatable :: recv_map(:)

    type(DataExchangeSendBuf_3D), allocatable :: send_buf(:)
    type(DataExchangeRecvBuf_3D), allocatable :: recv_buf(:)

    integer, allocatable :: send_attrib(:,:)
    integer, allocatable :: recv_attrib(:,:)
    integer :: n_recv, n_send
    integer :: i

    ! initialization ...........................................................

    n_send = old_mesh % n_child
    allocate(send_map(n_send), send_buf(n_send))
    do i = 1, n_send
      send_map(i) = DataExchangeMap_3D(old_mesh % map_child(i))
    end do

    n_recv = child % n_parent
    allocate(recv_map(n_recv), recv_buf(n_recv))
    do i = 1, n_recv
      recv_map(i) = DataExchangeMap_3D(child % map_parent(i))
    end do

    allocate(send_attrib(2, old_mesh % n_elem   ), source = -1)
    allocate(recv_attrib(2, child    % n_cluster), source = -1)

    ! prepare send data ........................................................

    do i = 1, old_mesh % n_elem
      send_attrib(1,i) = new_mesh % proc_part(dist_map % tp_elem(i)) ! new process
      send_attrib(2,i) = dist_map % id_elem(i)                       ! new ID
    end do

    ! pass parent attributes to children .......................................

    do i = 1, n_send
      call send_buf(i) % Extract_Data(send_map(i), send_attrib)
      call send_buf(i) % Send_Start()
    end do

    do i = 1, n_recv
      call recv_buf(i) % Init(recv_map(i), recv_attrib)
      call recv_buf(i) % Recv_Start()
    end do

    do i = 1, n_send
      call send_buf(i) % Send_Finish()
    end do

    do i = 1, n_recv
      call recv_buf(i) % Recv_Finish()
      call recv_buf(i) % Assign_Data(recv_attrib)
    end do

    ! extract received attributes ..............................................

    do i = 1, child % n_elem
      associate(element => child % element(i))
        element % adaptation % parent_proc = recv_attrib(1, element%cluster_id)
        element % adaptation % parent_id   = recv_attrib(2, element%cluster_id)
      end associate
    end do

    ! rebuild map from child to parent .........................................

    call child % BuildMapToParent()

  end subroutine UpdateChild

  !=============================================================================

end module Root_Mesh_Partitioning__3D
