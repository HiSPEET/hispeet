module Partition_Root_Mesh__3D
  use XMPI
  use Mesh__3D
  use Element_Transfer_Buffer__3D
  use Element_Distribution_Map__3D
  use Element_Distribution_Data__3D

  implicit none
  private

  public :: RootMeshPartitioningOptions_3D
  public :: RootMeshPartitioning_3D

  type RootMeshPartitioningOptions_3D
    integer :: n_parts   = 1         !< number of new partitions
    integer :: method    = 1         !< weighting method
    integer :: weight(4) = [1,0,0,0] !< element and connectivity weights
  contains
    procedure :: Bcast => Bcast_PartitioningOptions
  end type RootMeshPartitioningOptions_3D

contains

  !=============================================================================
  ! RootMeshPartitioningOptions_3D: type-bound procedures

  subroutine Bcast_PartitioningOptions(opt, root, comm)
    class(RootMeshPartitioningOptions_3D), intent(inout) :: opt !< options
    integer       , intent(in) :: root !< rank root process
    type(MPI_Comm), intent(in) :: comm !< MPI communicator

    call XMPI_Bcast(opt % n_parts , root, comm)
    call XMPI_Bcast(opt % method  , root, comm)
    call XMPI_Bcast(opt % weight  , root, comm)

  end subroutine Bcast_PartitioningOptions

  !=============================================================================
  ! Partitioning routines

  !-----------------------------------------------------------------------------
  !>

  subroutine RootMeshPartitioning_3D(opt, old_mesh, new_mesh, map)
    class(RootMeshPartitioningOptions_3D), intent(in) :: opt
    class(Mesh_3D), intent(in)  :: old_mesh
    class(Mesh_3D), intent(out) :: new_mesh
    class(ElementDistributionMap_3D), intent(out) :: map

    ! internal variables .......................................................

    type(ElementDistributionData_3D), allocatable :: send_data(:)
    type(ElementDistributionData_3D), allocatable :: recv_data(:)

    type(MeshAttributes_3D) :: attrib
    integer :: i

    !---------------------------------------------------------------------------
    !$omp master

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

    ! build distribution map ...................................................

    if (old_mesh % part >= 0) then
      call BuildDistributionMap(opt, old_mesh, map)
    end if

    ! distribute element data ..................................................

    ! prepare send data
    if (old_mesh % part >= 0) then
      allocate(send_data(0:new_mesh%n_parts-1))
      call BuildDistributionData(old_mesh, new_mesh, map, send_data)
    end if

    ! prepare receive data
    if (new_mesh % part >= 0) then
      allocate(recv_data(0:old_mesh%n_parts-1))
      do i = 0, old_mesh%n_parts-1
        recv_data(i) % proc = old_mesh % proc_active(i)
        recv_data(i) % comm = old_mesh % comm_world
      end do
    end if

    ! send dimensions (skip local part)
    if (old_mesh % part >= 0) then
      do i = 0, new_mesh%n_parts-1
        if (send_data(i) % proc == old_mesh % proc) cycle
        call send_data(i) % Send_Dimensions()
      end do
    end if

    ! receive dimensions (skip local part)
    if (new_mesh % part >= 0) then
      do i = 0, old_mesh%n_parts-1
        if (recv_data(i) % proc == old_mesh % proc) cycle
        call recv_data(i) % Recv_Dimensions()
      end do
    end if

    ! send data (skip local part)
    if (old_mesh % part >= 0) then
      do i = 0, new_mesh%n_parts-1
        if (send_data(i) % proc == old_mesh % proc) cycle
        call send_data(i) % Send_Data()
      end do
    end if

    ! receive data (skip local part)
    if (new_mesh % part >= 0) then
      do i = 0, old_mesh%n_parts-1
        if (recv_data(i) % proc == old_mesh % proc) cycle
        call recv_data(i) % Recv_Data()
      end do
    end if

    ! finalize transfer
    if (old_mesh % part >= 0) then
      do i = 0, new_mesh%n_parts-1
        if (send_data(i) % proc == old_mesh % proc) cycle
        call send_data(i) % Send_Finish()
      end do
    end if
    if (new_mesh % part >= 0) then
      do i = 0, old_mesh%n_parts-1
        if (recv_data(i) % proc == old_mesh % proc) cycle
        call recv_data(i) % Recv_Finish()
      end do
    end if

    ! nothing more to do if new partition is empty
    if (new_mesh % part < 0) return

    ! assign element data ......................................................

    ! number of new elements
    if (old_mesh % part >= 0) then
      new_mesh % n_elem = sum(recv_data % n_element) &
                        + send_data(old_mesh % part) % n_element
    else
      new_mesh % n_elem = sum(recv_data % n_element)
    end if

    allocate(new_mesh % element(new_mesh % n_elem))

    ! assign received data
    if (new_mesh % part >= 0) then
      do i = 0, old_mesh%n_parts-1
        if (recv_data(i) % proc /= old_mesh % proc) then
          call AssignDistributionData(new_mesh, recv_data(i))
        else
          call AssignDistributionData(new_mesh, send_data(old_mesh % part))
        end if
      end do
    end if

    ! complete mesh.............................................................

    call new_mesh % IdentifyVertices()
    call new_mesh % IdentifyEdges()
    call new_mesh % BuildFaces()
    call new_mesh % BuildBoundaryFaces()
    call new_mesh % BuildLinks()
    call new_mesh % BuildGhosts()
    call new_mesh % BuildCuboids()
    call new_mesh % IdentifyRanks()

    !---------------------------------------------------------------------------
    !$omp end master

  end subroutine RootMeshPartitioning_3D

  !-----------------------------------------------------------------------------
  !> Computes the new target partitions and element IDs

  subroutine BuildDistributionMap(opt, old_mesh, map)
    class(RootMeshPartitioningOptions_3D), intent(in) :: opt
    class(Mesh_3D), intent(in) :: old_mesh
    type(ElementDistributionMap_3D), intent(out) :: map

    integer, allocatable, target :: tp_elem(:)
    integer, contiguous, pointer :: tp_elem_val(:,:,:,:)
    type(ElementTransferBuffer_3D), allocatable, asynchronous :: tp_elem_buf

    allocate(tp_elem(old_mesh%n_elem + old_mesh%n_ghost), source = -1)

    ! graph-based partitioning
    call ParMETIS_Partitioner(opt, old_mesh, tp_elem)

    ! transfer target partitio ID to ghosts
    if (old_mesh % n_ghost > 0) then
      tp_elem_val(1:1,1:1,1:1,1:old_mesh%n_elem+old_mesh%n_ghost) => tp_elem
      tp_elem_buf = ElementTransferBuffer_3D(old_mesh, tp_elem_val)
      call tp_elem_buf % Transfer(old_mesh, tp_elem_val, tag = 1000)
      call tp_elem_buf % Merge(tp_elem_val)
    end if

    ! distribution map, assuming identity for proc_part
    map = ElementDistributionMap_3D(old_mesh, opt%n_parts, tp_elem)

  end subroutine BuildDistributionMap

  !-----------------------------------------------------------------------------
  !> Partitioning with ParMETIS

  subroutine ParMETIS_Partitioner(opt, old_mesh, tp_elem)
    use ParMETIS_Binding

    ! arguments ................................................................

    class(RootMeshPartitioningOptions_3D), intent(in) :: opt
    class(Mesh_3D), intent(in) :: old_mesh
    integer, intent(inout) :: tp_elem(:) !< element target partitions

    ! internal variables .......................................................

    ! ParMETIS arguments
    integer(METIS_IDX_T) :: wgtflag, numflag, ncon, nparts, edgecut
    integer(METIS_IDX_T) :: options(3) = 0
    integer(METIS_IDX_T), allocatable :: vtxdist(:), xadj(:), adjncy(:)
    integer(METIS_IDX_T), allocatable :: vwgt(:,:)
    integer(METIS_IDX_T), allocatable :: adjwgt(:)
    real(METIS_REAL_T),   allocatable :: tpwgts(:,:)
    real(METIS_REAL_T),   allocatable :: ubvec(:)
    integer(METIS_IDX_T), allocatable :: part(:)

    integer, allocatable :: nvtx_proc(:)        ! num graph vertices per process
    integer, allocatable, target :: vtx_elem(:) ! graph vertex IDs
    integer, pointer :: var_vtx_elem(:,:,:,:)   ! maps o element variable

    type(ElementTransferBuffer_3D), asynchronous :: buf_vtx_elem

    integer :: e, i, j, k, m, n

    !---------------------------------------------------------------------------
    ! Body, using ParMETIS-style naming

    associate( comm   => old_mesh % comm_parts  &
             , proc   => old_mesh % part        &
             , nproc  => old_mesh % n_parts     &
             , nvtx   => old_mesh % n_elem      &
             , weight => opt % weight           )

      ! initialization .........................................................

      ! use C-style numbering for ParMETIS
      numflag = 0

      ! set default to unweighted graph
      ncon    = 0
      wgtflag = 0

      ! element weighting
      if (weight(1) > 0) then
        wgtflag = 2 ! activate vertex constraints
        select case(opt % method)
        case(1)
          ncon = 1  ! consider only child level
        case(2)
          ncon = 2  ! consider also higher level descendants
        end select
      end if

      ! adjacency weighting
      if (any(weight(1:3) > 0)) then
        wgtflag = wgtflag + 1
      end if

      ! number of new partitions
      nparts = opt % n_parts

      ! ParMETIS array arguments, using C-style numbering
      allocate( vtxdist ( 0:nproc              ) )
      allocate( vwgt    ( 0:ncon-1, 0:nvtx-1   ) )
      allocate( tpwgts  ( 0:ncon-1, 0:nparts-1 ) )
      allocate( xadj    ( 0:nvtx               ) )
      allocate( ubvec   ( 0:ncon-1             ) )
      allocate( part    ( 0:nvtx-1             ) )

      ! auxiliary arrays
      allocate( nvtx_proc ( 0:nproc-1               ) )
      allocate( vtx_elem  ( 1:nvtx+old_mesh%n_ghost ) )

      ! graph vertex ID element variable and transfer buffer
      var_vtx_elem(1:1, 1:1, 1:1, 1:size(vtx_elem)) => vtx_elem
      buf_vtx_elem = ElementTransferBuffer_3D(old_mesh, var_vtx_elem)

      if (ncon > 0) then
        ubvec  = 1.05          ! tolerances for multi-constraint weighting
        tpwgts = 1.00 / nparts ! fractions of vertex weight per partition
      end if

      ! graph vertex distribution ..............................................

      ! build list of element counts per partition
      m = max(nvtx, 0)
      call MPI_Allgather(m, 1, MPI_INTEGER, nvtx_proc, 1, MPI_INTEGER, comm)

      ! compute graph vertex offsets
      vtxdist(0) = 0
      do i = 1, nproc
        vtxdist(i) = vtxdist(i-1) + nvtx_proc(i-1)
      end do

      ! graph vertex IDs by element index (local+ghost) ........................

      ! graph vertex ID of local elements
      m = vtxdist(proc) - 1
      do e = 1, nvtx
        vtx_elem(e) = e + m
      end do

      ! transfer graph vertex IDs to ghosts
      call buf_vtx_elem % Transfer(old_mesh, var_vtx_elem, tag=1000)
      call buf_vtx_elem % Merge(var_vtx_elem)

      ! graph vertex weights and adjacency offsets .............................

      xadj(0) = 0

      do e = 1, nvtx
        associate(element => old_mesh % element(e))

          i = e - 1

          select case(ncon)
          case(1)
            vwgt(0,i) = weight(1)
          case(2)
            select case(element % adaptation % mark)
            case(1)
              m = 8     ! = 8^1
            case(2)
              m = 64    ! = 8^2
            case(3)
              m = 576   ! = 8^2 + 8^3
            case(4:)
              m = 4672  ! = 8^2 + 8^3 + 8^4
            case default
              m = 0
            end select
            vwgt(0,i) = weight(1)
            vwgt(1,i) = weight(1) * m
          end select

          i = i + 1

          xadj(i) = xadj(i-1)
          do k = 1, 6
            xadj(i) = xadj(i) + element % face(k) % n_neighbor
          end do

          if (weight(3) > 0) then ! include element-edge neighbors
            do k = 1, 12
              xadj(i) = xadj(i) + element % edge(k) % n_neighbor
            end do
          end if

          if (weight(4) > 0) then ! include element-vertex neighbors
            do k = 1, 8
              xadj(i) = xadj(i) + element % vertex(k) % n_neighbor
            end do
          end if

        end associate
      end do

      ! adjacency and adjacency weights ........................................

      m = xadj(nvtx) - 1 ! number of graph edges
      allocate(adjncy(0:m))
      if (any(weight(1:3) > 0)) then
        allocate(adjwgt(0:m))
      else
        allocate(adjwgt(0:0))
      end if

      m = 0

      do e = 1, nvtx
        associate(element => old_mesh % element(e))

          ! element faces
          do k = 1, 6
            n = element % face(k) % n_neighbor - 1
            if (n < 0) cycle
            i = element % face(k) % i_neighbor
            do j = i, i+n
              adjncy(m) = vtx_elem(element % neighbor(j) % id)
              if (weight(2) > 0) then
                adjwgt(m) = weight(2)
              end if
              m = m + 1
            end do
          end do

          ! element edges
          if (weight(3) > 0) then
            do k = 1, 12
              n = element % edge(k) % n_neighbor - 1
              if (n < 0) cycle
              i = element % edge(k) % i_neighbor
              do j = i, i+n
                adjncy(m) = vtx_elem(element % neighbor(j) % id)
                adjwgt(m) = weight(3)
                m = m + 1
              end do
            end do
          end if

          ! element vertices
          if (weight(4) > 0) then
            do k = 1, 8
              n = element % vertex(k) % n_neighbor - 1
              if (n < 0) cycle
              i = element % vertex(k) % i_neighbor
              do j = i, i+n
                adjncy(m) = vtx_elem(element % neighbor(j) % id)
                adjwgt(m) = weight(4)
                m = m + 1
              end do
            end do
          end if

        end associate
      end do

      ! ParMETIS ...............................................................

      call ParMETIS_V3_PartKway( vtxdist, xadj, adjncy, vwgt, adjwgt, wgtflag  &
                               , numflag, ncon, nparts, tpwgts, ubvec, options &
                               , edgecut, part, comm % MPI_VAL                 )

      ! result .................................................................

      tp_elem(1:nvtx) = part

    end associate

    !---------------------------------------------------------------------------

  end subroutine ParMETIS_Partitioner

  !=============================================================================
  ! Distribution routines

  !-----------------------------------------------------------------------------
  !> Build element distribution data fro sending from old to new partitions

  subroutine BuildDistributionData(old_mesh, new_mesh, map, dist_data)

    class(Mesh_3D), intent(in) :: old_mesh
    class(Mesh_3D), intent(in) :: new_mesh
    class(ElementDistributionMap_3D) , intent(in)  :: map
    class(ElementDistributionData_3D), intent(out) :: dist_data(0:)

    ! internal variables .......................................................

    integer :: ce_part(0:map%n_parts-1)  ! element counters
    integer :: cn_part(0:map%n_parts-1)  ! neighbor counters
    integer :: cp_part(0:map%n_parts-1)  ! element point counters

    integer :: ce, cn, cp, nn, np, po, tp
    integer :: i, j, k, l

    if (old_mesh % part < 0) return

    ! MPI and dimensions .......................................................

    ! basic initialization
    do tp = 0, map%n_parts-1
      dist_data(tp) % proc       = new_mesh % proc_part(tp)
      dist_data(tp) % comm       = old_mesh % comm_world
      dist_data(tp) % n_element  = map % ne_part(tp)
      dist_data(tp) % n_neighbor = 0
      dist_data(tp) % n_point    = 0
    end do

    ! counting neighbors and element points
    do l = 1, old_mesh % n_elem

      tp = map % tp_elem(l)
      nn = size(old_mesh % element(l) % neighbor)
      np = (old_mesh % element(l) % geometry % po + 1) ** 3

      dist_data(tp) % n_neighbor  =  dist_data(tp) % n_neighbor  +  nn
      dist_data(tp) % n_point     =  dist_data(tp) % n_point     +  np

    end do

    ! build components .........................................................

    ! allocate components
    do tp = 0, map%n_parts-1
      associate(dd => dist_data(tp))
        allocate( dd % element              ( dd % n_element  ) )
        allocate( dd % neighbor_id          ( dd % n_neighbor ) )
        allocate( dd % neighbor_part        ( dd % n_neighbor ) )
        allocate( dd % neighbor_component   ( dd % n_neighbor ) )
        allocate( dd % neighbor_orientation ( dd % n_neighbor ) )
        allocate( dd % x_e                  ( dd % n_point, 3 ) )
      end associate
    end do

    ! initialize counters
    ce_part = 0
    cn_part = 0
    cp_part = 0

    ! extract and transform element data
    do l = 1, old_mesh % n_elem
      tp = map % tp_elem(l)
      associate(dd => dist_data(tp))

        ! static element components
        associate(element => old_mesh % element(l))
          ce = ce_part(tp) + 1
          dd % element(ce) % id              = map % id_elem(l)
          dd % element(ce) % vertex          = element % vertex
          dd % element(ce) % edge            = element % edge
          dd % element(ce) % face            = element % face
          dd % element(ce) % adaptation      = element % adaptation
          dd % element(ce) % geometry % po   = element % geometry % po
          dd % element(ce) % geometry % x_c  = element % geometry % x_c
          dd % element(ce) % geometry % dx_m = element % geometry % dx_m
          ce_part(tp) = ce
        end associate

        ! neighbor data
        associate(neighbor => old_mesh % element(l) % neighbor)
          cn = cn_part(tp)
          nn = size(neighbor)
          do i = 1, nn
            cn = cn + 1
            dd % neighbor_id          (cn) = map % id_elem( neighbor(i)% id )
            dd % neighbor_part        (cn) = map % tp_elem( neighbor(i)% id )
            dd % neighbor_component   (cn) = neighbor(i) % component
            dd % neighbor_orientation (cn) = neighbor(i) % orientation
          end do
          cn_part(tp) = cn
        end associate

        ! element points
        associate(x_e => old_mesh % element(l) % geometry % x_e)
          po = old_mesh % element(l) % geometry % po
          cp = cp_part(tp)
          do k = 0, po
          do j = 0, po
          do i = 0, po
            cp = cp + 1
            dd % x_e(cp,1) = x_e(i,j,k,1)
            dd % x_e(cp,2) = x_e(i,j,k,2)
            dd % x_e(cp,3) = x_e(i,j,k,3)
          end do
          end do
          end do
          cp_part(tp) = cp
        end associate

      end associate
    end do

  end subroutine BuildDistributionData

  !-----------------------------------------------------------------------------
  !> Assign element distribution data to new mesh partition

  subroutine AssignDistributionData(new_mesh, dist_data)
    class(Mesh_3D), intent(inout) :: new_mesh
    type(ElementDistributionData_3D), intent(in) :: dist_data

    integer :: cn, cp, nn, po
    integer :: e, i, j, k, l

    ! initialization ...........................................................

    cn = 0 ! neighbor counter
    cp = 0 ! point counter

    do l = 1, dist_data % n_element

      e = dist_data % element(l) % id

      ! static element components ..............................................

      new_mesh % element(e) = dist_data % element(l)

      ! neighbors ..............................................................

      nn = sum( new_mesh % element(e) % face   % n_neighbor ) &
         + sum( new_mesh % element(e) % edge   % n_neighbor ) &
         + sum( new_mesh % element(e) % vertex % n_neighbor )

      allocate(new_mesh % element(e) % neighbor(nn))

      associate(neighbor => new_mesh % element(e) % neighbor)
        do i = 1, nn
          cn = cn + 1
          neighbor(i) % id          = dist_data % neighbor_id          (cn)
          neighbor(i) % part        = dist_data % neighbor_part        (cn)
          neighbor(i) % component   = dist_data % neighbor_component   (cn)
          neighbor(i) % orientation = dist_data % neighbor_orientation (cn)
        end do
      end associate

      ! element points .........................................................

      po = new_mesh % element(e) % geometry % po

      allocate(new_mesh % element(e) % geometry % x_e(0:po,0:po,0:po,3))

      associate(x_e => new_mesh % element(e) % geometry % x_e)
        do k = 0, po
        do j = 0, po
        do i = 0, po
          cp = cp + 1
          x_e(i,j,k,1) = dist_data % x_e(cp, 1)
          x_e(i,j,k,2) = dist_data % x_e(cp, 2)
          x_e(i,j,k,3) = dist_data % x_e(cp, 3)
        end do
        end do
        end do
      end associate

    end do

  end subroutine AssignDistributionData

  !=============================================================================

end module Partition_Root_Mesh__3D
