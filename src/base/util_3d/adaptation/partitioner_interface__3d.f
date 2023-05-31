module Partitioner_Interface__3D
  use XMPI
  use Mesh__3D
  use Element_Transfer_Buffer__3D

  implicit none
  private

  public :: PartitioningOptions_3D
  public :: ParMETIS_Partitioner_3D

  !-----------------------------------------------------------------------------
  !> Partitioning options
  !>
  !> For partitioning, an abstract graph is formed. The vertices of this graph
  !> represent the mesh elements, while its edges correspond to connections
  !> between elements. Depending to the adjacency type, the graph edges are
  !> classified into three groups: element face, element edge and element
  !> vertex connections.
  !>
  !> The number of partitions is set via `n_parts`.
  !>
  !> The `mode` option activates one of the following partitioning modes:
  !>
  !>   - `1` single-level
  !>   - `2` parent
  !>   - `3` child
  !>   - `4` weighted
  !>
  !> In single-level mode, all graph vertices are weighted equally. In parent
  !> mode, a second constraint is derived from `mesh%element%adaptation%mark`.
  !> The child mode is intended for generating and distributing a new child
  !> mesh. In this case the constraint is identical to the second constraint
  !> in the parent mode. In the weighted mode, precomputed weights are passed
  !> by an additional argument.
  !>
  !> The `weight` option assigns relative weights to the graph components:
  !>
  !>   - `weight(1) ≥ 1` graph vertices, i.e. mesh elements
  !>   - `weight(2) ≥ 0` graph edges emerging from element face connections
  !>   - `weight(3) ≥ 0` graph edges emerging from element edge connections
  !>   - `weight(4) ≥ 0` graph edges emerging from element vertex connections
  !>
  !> The vertex weight is multiplied by the workload determined by the selected
  !> partitioning mode. Vertices with no workload are excluded from the graph.
  !> Connections that have no weight can be cut without penalty.
  !> Note that the lower bounds of the weights are enforced for regularity.

  type PartitioningOptions_3D
    integer :: n_parts   = 1         !< number of new partitions
    integer :: mode      = 1         !< partitioning mode
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

    call XMPI_Bcast(opt % n_parts   , root, comm)
    call XMPI_Bcast(opt % mode      , root, comm)
    call XMPI_Bcast(opt % weight    , root, comm)

  end subroutine Bcast_PartitioningOptions

  !=============================================================================
  ! Partitioner interfaces

  !-----------------------------------------------------------------------------
  !> ParMETIS interface

  subroutine ParMETIS_Partitioner_3D(opt, mesh, tp_elem)
    use ParMETIS_Binding

    ! arguments ................................................................

    class(PartitioningOptions_3D), intent(in) :: opt
    class(Mesh_3D), intent(in) :: mesh
    integer, intent(inout) :: tp_elem(:) !< child element target partitions

    ! internal variables .......................................................

    type(MPI_Comm) :: comm
    integer :: proc, nproc

    ! ParMETIS arguments
    integer(METIS_IDX_T) :: wgtflag, numflag, ncon, nparts, edgecut
    integer(METIS_IDX_T) :: options(3) = 0
    integer(METIS_IDX_T), allocatable :: vtxdist(:), xadj(:), adjncy(:)
    integer(METIS_IDX_T), allocatable :: vwgt(:,:)
    integer(METIS_IDX_T), allocatable :: adjwgt(:)
    real(METIS_REAL_T),   allocatable :: tpwgts(:,:)
    real(METIS_REAL_T),   allocatable :: ubvec(:)
    integer(METIS_IDX_T), allocatable :: part(:)

    integer :: nvtx
    integer, allocatable :: nvtx_proc(:)        ! num graph vertices per process
    integer, allocatable, target :: vtx_elem(:) ! graph vertex IDs
    integer, pointer :: var_vtx_elem(:,:,:,:)   ! map to element variable

    type(ElementTransferBuffer_3D), asynchronous :: buf_vtx_elem

    integer :: e, i, j, k, l, m, n

    !---------------------------------------------------------------------------
    ! Body

    associate( n_elem  => mesh % n_elem   &
             , n_ghost => mesh % n_ghost  &
             , weight  => opt % weight           )

      ! initialization .........................................................

      ! vertices = contributing elements
      allocate(vtx_elem(n_elem + n_ghost), source = -1)
      i = 0
      do e = 1, n_elem
        if (opt%mode == 3 .and. mesh%element(e)%adaptation%mark < 1) cycle
        vtx_elem(e) = i
        i = i + 1
      end do
      nvtx = i

      ! create MPI communicator comprising all parent meshes with nvtx > 0
      if (nvtx > 0) then
        m = 1
      else
        m = 0
      end if
      call MPI_Comm_split(mesh%comm_parts, m, mesh%part, comm)

      call MPI_Comm_rank(comm, proc)
      call MPI_Comm_size(comm, nproc)

      ! graph vertex ID element variable and transfer buffer
      var_vtx_elem(1:1, 1:1, 1:1, 1:size(vtx_elem)) => vtx_elem
      buf_vtx_elem = ElementTransferBuffer_3D(mesh, var_vtx_elem)

      ! ParMetis input arguments ...............................................

      ! use C-style numbering for ParMETIS
      numflag = 0

      ! number of new partitions
      nparts = opt % n_parts

      ! number of constraints
      select case(opt % mode)
      case(1)
        ncon = 1
      case(2)
        ncon = 2
      case(3)
        ncon = 1
      end select

      ! weight flag
      if (any(weight(1:3) > 0)) then
        ! vertex and edge constraints
        wgtflag = 3
      else
        ! vertex constraints only
        wgtflag = 2
      end if

      ! ParMETIS array arguments, using C-style numbering
      allocate( vtxdist ( 0:nproc              ) )
      allocate( vwgt    ( 0:ncon-1, 0:nvtx-1   ) )
      allocate( tpwgts  ( 0:ncon-1, 0:nparts-1 ) )
      allocate( xadj    ( 0:nvtx               ) )
      allocate( ubvec   ( 0:ncon-1             ) )
      allocate( part    ( 0:nvtx-1             ) )

      ! tolerance for multi-constraint weighting
      ubvec  = 1.05

      ! fractions of vertex weight per partition
      tpwgts = 1.00 / nparts

      ! graph vertex distribution ..............................................

      ! build list of graph vertex counts per partition
      m = max(nvtx, 0)
      allocate(nvtx_proc(0:nproc-1))
      call MPI_Allgather(m, 1, MPI_INTEGER, nvtx_proc, 1, MPI_INTEGER, comm)

      ! compute graph vertex offsets
      vtxdist(0) = 0
      do i = 1, nproc
        vtxdist(i) = vtxdist(i-1) + nvtx_proc(i-1)
      end do

      ! graph vertex IDs by element index (local+ghost) ........................

      ! graph vertex ID of local parent elements
      where(vtx_elem >= 0)
        vtx_elem = vtx_elem + vtxdist(proc)
      end where

      ! transfer graph vertex IDs to ghosts
      call buf_vtx_elem % Transfer(mesh, var_vtx_elem, tag=1000)
      call buf_vtx_elem % Merge(var_vtx_elem)

      ! graph vertex weights and adjacency offsets .............................

      xadj(0) = 0

      i = 0
      do e = 1, n_elem
        if (vtx_elem(e) < 0) cycle
        associate(element => mesh % element(e))

          ! number of children
          select case(element % adaptation % mark)
          case(1:6)
            m = 4      ! face refined
          case(7:18)
            m = 2      ! edge refined
          case(19:26)
            m = 1      ! vertex refined
          case default
            m = 8      ! regular refinement
          end select

          ! vertex weights
          select case(opt % mode)
          case(1)
            vwgt(0,i) = weight(1)
          case(1)
            vwgt(0,i) = weight(1)
            vwgt(1,i) = weight(1) * m
          case(2)
            vwgt(0,i) = weight(1) * m
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

      do e = 1, n_elem
        if (vtx_elem(e) < 0) cycle
        associate(element => mesh % element(e))

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

      i = 0
      do e = 1, n_elem
        if (vtx_elem(e) >= 0) then
          tp_elem(e) = part(i)
          i = i + 1
        else
          tp_elem(e) = -1
        end if
      end do

      ! clean-up ...............................................................

      call MPI_Comm_free(comm)

    end associate

    !---------------------------------------------------------------------------

  end subroutine ParMETIS_Partitioner_3D

  !=============================================================================

end module Partitioner_Interface__3D
