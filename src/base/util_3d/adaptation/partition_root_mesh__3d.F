module Partition_Root_Mesh__3D
  use XMPI
  use Mesh__3D
  use Element_Transfer_Buffer__3D
  use Element_Distribution_Data__3D

  implicit none
  private

  public :: PartitioningOptions_3D

!## CHECK
public :: ComputePartitions
!## CHECK END

  type PartitioningOptions_3D
    integer :: method    = 1         !< weighting method
    integer :: weight(4) = [1,0,0,0] !< element and connectivity weights
  contains
    procedure :: Bcast => Bcast_PartitioningOptions
  end type PartitioningOptions_3D

contains

  !=============================================================================
  ! PartitioningOptions_3D: type-bound procedures

  subroutine Bcast_PartitioningOptions(opt, root, comm)
    class(PartitioningOptions_3D), intent(inout) :: opt !< options
    integer       , intent(in) :: root !< rank root process
    type(MPI_Comm), intent(in) :: comm !< MPI communicator

    call XMPI_Bcast(opt % method , root, comm)
    call XMPI_Bcast(opt % weight , root, comm)

  end subroutine Bcast_PartitioningOptions

  !=============================================================================
  ! Partitioning routines

  !-----------------------------------------------------------------------------
  !> Computes the new target partitions and element IDs
  !>
  !> @note
  !>   Runs only on the master thread with OpenMP

  subroutine ComputePartitions(mesh, opt, n_parts, tp_elem, id_elem, ne_part)

    ! arguments ................................................................

    class(Mesh_3D), intent(in) :: mesh
    !< local mesh partition

    class(PartitioningOptions_3D), intent(in) :: opt
    !< partitioning options

    integer, intent(in) :: n_parts
    !< number of target partitions

    integer, intent(out) :: tp_elem(1:mesh%n_elem+mesh%n_ghost)
    !< element target partitions

    integer, intent(out) :: id_elem(1:mesh%n_elem+mesh%n_ghost)
    !< element IDs in target partitions

    integer, asynchronous, intent(out) :: ne_part(0:n_parts-1)
    !< number of elements contributed to target partitions

    ! internal variables .......................................................

    integer, contiguous, pointer :: tp_elem_(:)
    integer, contiguous, pointer :: id_elem_(:)
    integer, allocatable, target :: elem_val(:,:,:,:,:)
    type(ElementTransferBuffer_3D), asynchronous :: elem_buf

    !---------------------------------------------------------------------------
    !$omp master

    ! workspace ................................................................

    allocate(elem_val(1, 1, 1, mesh%n_elem + mesh%n_ghost, 2), source = -1)

    elem_buf =  ElementTransferBuffer_3D(mesh, elem_val)
    tp_elem_ => elem_val(1,1,1,:,1)
    id_elem_ => elem_val(1,1,1,:,2)

    ! target partitions ........................................................

    call ParMETIS_Partitioner(mesh, opt, n_parts, tp_elem_)

    ! element distribution .....................................................

    call ComputeElementDistribution(mesh, n_parts, tp_elem_, id_elem_, ne_part)

    call elem_buf % Transfer(mesh, elem_val, tag = 1000)
    call elem_buf % Merge(elem_val)

    ! assign results ...........................................................

    tp_elem = tp_elem_
    id_elem = id_elem_

    !---------------------------------------------------------------------------
    !$omp end master

  end subroutine ComputePartitions

  !-----------------------------------------------------------------------------
  !> Partitioning with ParMETIS

  subroutine ParMETIS_Partitioner(mesh, opt, n_parts, tp_elem)
    use ParMETIS_Binding

    ! arguments ................................................................

    class(Mesh_3D), intent(in) :: mesh
    !< local mesh partition

    class(PartitioningOptions_3D), intent(in) :: opt
    !< partitioning options

    integer, intent(in)  :: n_parts
    !< number of target partitions

    integer, intent(out) :: tp_elem(:)
    !< element target partitions

    ! local variables ..........................................................

    ! ParMETIS arguments
    integer(METIS_IDX_T) :: wgtflag, numflag, ncon, nparts, edgecut
    integer(METIS_IDX_T) :: options(3) = 0
    integer(METIS_IDX_T), allocatable :: vtxdist(:), xadj(:), adjncy(:)
    integer(METIS_IDX_T), allocatable :: vwgt(:,:)
    integer(METIS_IDX_T), allocatable :: adjwgt(:)
    real(METIS_REAL_T),   allocatable :: tpwgts(:,:)
    real(METIS_REAL_T),   allocatable :: ubvec(:)
    integer(METIS_IDX_T), allocatable :: part(:)

    integer, allocatable :: n_elem_proc(:)   ! num elements per process
    integer, allocatable, target :: vtxid(:) ! graph vertex IDs
    integer, pointer :: var_vtxid(:,:,:,:)   ! maps vtxid to element variable

    type(ElementTransferBuffer_3D), asynchronous :: buf_vtxid

    integer :: e, i, j, k, m, n

    !$omp master

    associate( comm   => mesh % comm_parts &
             , i_proc => mesh % part       & ! i_proc = present partition ID
             , n_proc => mesh % n_parts    & ! n_proc = present num partitions
             , n_elem => mesh % n_elem     &
             , weight => opt  % weight     )

      ! initialization .........................................................

      ! use C-style numbering for ParMETIS
      numflag = 0

      ! set default to unweighted graph
      ncon    = 0
      wgtflag = 0

      ! element weighting
      if (weight(1) > 0) then
        select case(opt % method)
        case(2)
          ncon = 2
        case(3)
          ncon = 3
        case default
          ncon = 1
        end select
        wgtflag = 2
      end if

      ! adjacency weighting
      if (any(weight(1:3) > 0)) then
        wgtflag = wgtflag + 1
      end if

      ! number of partitions
      nparts = n_parts

      ! ParMETIS array arguments, using C-style numbering
      allocate( vtxdist     ( 0:n_proc                ) )
      allocate( vwgt        ( 0:ncon-1, 0:n_elem -1   ) )
      allocate( tpwgts      ( 0:ncon-1, 0:n_parts-1   ) )
      allocate( xadj        ( 0:n_elem                ) )
      allocate( ubvec       ( 0:ncon-1                ) )
      allocate( part        ( 0:n_elem-1              ) )

      ! auxiliary arrays
      allocate( n_elem_proc ( 0:n_proc-1              ) )
      allocate( vtxid       ( 1:n_elem + mesh%n_ghost ) )

      ! graph vertex ID element variable and transfer buffer
      var_vtxid(1:1, 1:1, 1:1, 1:size(vtxid)) => vtxid
      buf_vtxid = ElementTransferBuffer_3D(mesh, var_vtxid)

      if (ncon > 0) then
        ubvec  = 1.05          ! tolerances for multi-constraint weighting
        tpwgts = 1.00 / nparts ! fractions of vertex weight per partition
      end if

      ! graph vertex distribution ..............................................

      ! build list of element counts per partition
      m = max(n_elem, 0)
      call MPI_Allgather(m, 1, MPI_INTEGER, n_elem_proc, 1, MPI_INTEGER, comm)

      ! compute graph vertex offsets
      vtxdist(0) = 0
      do i = 1, n_proc
        vtxdist(i) = vtxdist(i-1) + n_elem_proc(i-1)
      end do

      ! graph vertex IDs by element index (local+ghost) ........................

      ! graph vertex ID of local elements
      m = vtxdist(i_proc) - 1
      do e = 1, n_elem
        vtxid(e) = e + m
      end do

      ! transfer graph vertex IDs to ghosts
      call buf_vtxid % Transfer(mesh, var_vtxid, tag=1000)
      call buf_vtxid % Merge(var_vtxid)

      ! graph vertex weights and adjacency offsets .............................

      xadj(0) = 0

      do e = 1, mesh % n_elem
        associate(element => mesh % element(e))

          i = e - 1

          if (ncon > 0) then
            vwgt(0, i) = weight(1)
            if (ncon > 1) then
              select case(element % adaptation % mark)
              case(0)
                vwgt(1, i) = weight(1)
              case(1)
                vwgt(1, i) = weight(1) * 8
              case default
                vwgt(1, i) = 0
              end select
              if (ncon > 2) then
                select case(element % adaptation % levels)
                case(2)
                  i = 64    ! = 8^2
                case(3)
                  i = 576   ! = 8^2 + 8^3
                case(4:)
                  i = 4672  ! = 8^2 + 8^3 + 8^4
                case default
                  i = 0
                end select
                vwgt(2, i) = weight(1) * i
              end if
            end if
          end if

          xadj(e) = xadj(e-1)
          do k = 1, 6
            xadj(e) = xadj(e) + element % face(k) % n_neighbor
          end do

          if (weight(3) > 0) then ! include element-edge neighbors
            do k = 1, 12
              xadj(e) = xadj(e) + element % edge(k) % n_neighbor
            end do
          end if

          if (weight(4) > 0) then ! include element-vertex neighbors
            do k = 1, 8
              xadj(e) = xadj(e) + element % vertex(k) % n_neighbor
            end do
          end if

        end associate
      end do

      ! adjacency and adjacency weights ........................................

      m = xadj(n_elem) - 1 ! number of graph edges
      allocate(adjncy(0:m))
      if (any(weight(1:3) > 0)) then
        allocate(adjwgt(0:m))
      else
        allocate(adjwgt(0:0))
      end if

      m = 0

      do e = 1, mesh % n_elem
        associate(element => mesh % element(e))

          ! element faces
          do k = 1, 6
            n = element % face(k) % n_neighbor - 1
            if (n < 0) cycle
            i = element % face(k) % i_neighbor
            do j = i, i+n
              adjncy(m) = vtxid(element % neighbor(j) % id)
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
                adjncy(m) = vtxid(element % neighbor(j) % id)
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
                adjncy(m) = vtxid(element % neighbor(j) % id)
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

      tp_elem = part

    end associate

    !$omp end master

  end subroutine ParMETIS_Partitioner

  !-----------------------------------------------------------------------------
  !> Compute element IDs in their target partitions

  subroutine ComputeElementDistribution(mesh, n_parts, tp_elem, id_elem, ne_part)
    class(Mesh_3D), intent(in) :: mesh  !< local mesh partition
    integer, intent(in)  :: n_parts     !< number of target partitions
    integer, intent(in)  :: tp_elem(1:) !< element target partitions
    integer, intent(out) :: id_elem(1:) !< element IDs in target partition
    integer, intent(out) :: ne_part(0:) !< element IDs in target partition

    integer, allocatable, asynchronous :: id_part(:)

    integer :: e, p

    !$omp master

    allocate(id_part(0:n_parts-1))

    ne_part = 0
    do e = 1, mesh % n_elem
      ne_part(tp_elem(e)) = ne_part(tp_elem(e)) + 1
    end do

    call ComputeElementOffsets(mesh, ne_part, id_part)

    do e = 1, mesh % n_elem
      p = tp_elem(e)               ! fetch target partition
      id_part(p) = id_part(p) + 1  ! increment element counter
      id_elem(e) = id_part(p)      ! set element ID in target partition
    end do

    !$omp end master

  end subroutine ComputeElementDistribution

  !-----------------------------------------------------------------------------
  !> Computes the offsets for numbering the redistributed elements

  subroutine ComputeElementOffsets(mesh, ne_part, id_part)
    class(Mesh_3D), intent(in) :: mesh
    integer, asynchronous, intent(in)  :: ne_part(0:)
    integer, asynchronous, intent(out) :: id_part(0:)

    type(MPI_Win) :: window
    integer(MPI_ADDRESS_KIND) :: integer_extent, lb
    integer(MPI_ADDRESS_KIND) :: buf_size
    integer(MPI_ADDRESS_KIND) :: target_disp = 0
    integer :: disp_unit
    integer :: i, np

    ! preliminaries ............................................................

    np = size(ne_part)

    ! create MPI window ........................................................

    call MPI_Type_get_extent(MPI_INTEGER, lb, integer_extent)
    buf_size  = integer_extent * np
    disp_unit = int(integer_extent)

    call MPI_Win_create( base       =  id_part            &
                       , size       =  buf_size           &
                       , disp_unit  =  disp_unit          &
                       , info       =  MPI_INFO_NULL      &
                       , comm       =  mesh % comm_parts  &
                       , win        =  window             )

    ! compute offsets via accumulation .........................................

    id_part = 0

    call MPI_Win_fence(MPI_MODE_NOSTORE + MPI_MODE_NOPRECEDE, window)

    do i = mesh%part + 1, mesh%n_parts - 1
      call MPI_Accumulate( origin_addr     = ne_part       &
                         , origin_count    = np            &
                         , origin_datatype = MPI_INTEGER   &
                         , target_rank     = i             &
                         , target_disp     = target_disp   &
                         , target_count    = np            &
                         , target_datatype = MPI_INTEGER   &
                         , op              = MPI_SUM       &
                         , win             = window        )
    end do

    call MPI_Win_fence( MPI_MODE_NOSTORE + MPI_MODE_NOPUT + MPI_MODE_NOSUCCEED &
                      , window )

    call MPI_Win_free(window)

  end subroutine ComputeElementOffsets

  !=============================================================================
  ! Distribution routines

  !-----------------------------------------------------------------------------
  !> Build element distribution data

  subroutine BuildDistributionData( mesh, n_parts, tp_elem, id_elem, ne_part, &
                                    proc_part, dist_data )

    ! arguments ................................................................

    class(Mesh_3D), intent(in) :: mesh
    !< local mesh partition

    integer, intent(in) :: n_parts
    !< number of target partitions

    integer, intent(in) :: tp_elem(1:mesh%n_elem+mesh%n_ghost)
    !< element target partitions

    integer, intent(in) :: id_elem(1:mesh%n_elem+mesh%n_ghost)
    !< element IDs in target partitions

    integer, intent(in) :: ne_part(0:n_parts-1)
    !< number of elements contributed to target partitions

    integer, optional, intent(in) :: proc_part(0:n_parts-1)
    !< maps partition ID to precess ID in `mesh%comm_world` [identity]

    type(ElementDistributionData_3D), intent(out) :: dist_data(0:n_parts-1)
    !< local element distribution data being prepared for sending

    ! internal variables .......................................................

    integer :: ie_part(0:n_parts-1)  ! element counters
    integer :: in_part(0:n_parts-1)  ! neighbor counters
    integer :: ip_part(0:n_parts-1)  ! element point counters

    integer :: e, i, j, nn, np, tp

    if (mesh % part < 0) return

    !---------------------------------------------------------------------------
    !$omp master

    ! MPI and dimensions .......................................................

    ! basic initialization
    do tp = 0, n_parts-1
      dist_data(tp) % proc       = tp
      dist_data(tp) % comm       = mesh % comm_world
      dist_data(tp) % n_element  = ne_part(tp)
      dist_data(tp) % n_neighbor = 0
      dist_data(tp) % n_point    = 0
    end do

    ! optional transformation of part IDs to process IDs
    if (present(proc_part)) then
      dist_data % proc = proc_part( dist_data % proc )
    end if

    ! counting neighbors and element points
    do e = 1, mesh % n_elem

      tp = tp_elem(e)
      nn = size(mesh % element(e) % neighbor)
      np = size(mesh % element(e) % geometry % x_e, 1)

      dist_data(tp) % n_neighbor  =  dist_data(tp) % n_neighbor  +  nn
      dist_data(tp) % n_point     =  dist_data(tp) % n_point     +  np

    end do

    ! built components .........................................................

    ! allocate components
    do tp = 0, n_parts-1
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
    ie_part = 0
    in_part = 0
    ip_part = 0

    ! extract and transform element data
    do e = 1, mesh % n_elem
      tp = tp_elem(e)
      associate(element => mesh % element(e), dd => dist_data(tp))

        ! static element components
        i = ie_part(tp) + 1
        dd % element(i) % id              = id_elem(e)
        dd % element(i) % vertex          = element % vertex
        dd % element(i) % edge            = element % edge
        dd % element(i) % face            = element % face
        dd % element(i) % adaptation      = element % adaptation
        dd % element(i) % geometry % po   = element % geometry % po
        dd % element(i) % geometry % x_c  = element % geometry % x_c
        dd % element(i) % geometry % dx_m = element % geometry % dx_m
        ie_part(tp) = i

        ! neighbor data
        i  = in_part(tp)
        nn = size(element % neighbor)
        do j = 1, nn
          i = i + 1
          dd % neighbor_id          (i) = id_elem( element%neighbor(j)% id   )
          dd % neighbor_part        (i) = tp_elem( element%neighbor(j)% part )
          dd % neighbor_component   (i) = element % neighbor(j) % component
          dd % neighbor_orientation (i) = element % neighbor(j) % orientation
        end do
        in_part(tp) = i

        ! element points
        i  = ip_part(tp)
        np = size(element % geometry % x_e, 1)
        dd % x_e(i+1:i+np, 1) = reshape(element % geometry % x_e(:,:,:,1), [np])
        dd % x_e(i+1:i+np, 2) = reshape(element % geometry % x_e(:,:,:,2), [np])
        dd % x_e(i+1:i+np, 3) = reshape(element % geometry % x_e(:,:,:,3), [np])
        ip_part(tp) = i + np

      end associate
    end do

    !---------------------------------------------------------------------------
    !$omp end master

  end subroutine BuildDistributionData

  !=============================================================================

end module Partition_Root_Mesh__3D
