module Partitioner_Interface__3D
  use XMPI
  use Logging_Levels
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
  !> Set `child` to partition the child mesh. Otherwise the given `mesh` will be
  !> repartitioned.
  !>
  !> Set `split` to obtain child partitions by subdividing parent partitions.
  !> Only used with `child = T`.
  !>
  !> Set `lazy` to adopt the parent partitioning if the number of parts remains
  !> the same.
  !>
  !> Set `n_parts` to the requested number of partitions
  !>
  !> The partitioning is subject to `n_con` constraints which depend on the
  !> target mesh (given or child) and the number of sublevels considered.
  !> For any element, let be `c` the number of children and `g` the number of
  !> grand children and further descendants in sublevels 2 up to `n_sub`. Then
  !> the constraints are defined as follows:
  !>
  !>   | n_con | given mesh  | child mesh |
  !>   | ----- | ----------- | ---------- |
  !>   |   1   |  `1`        |   `c+g`    |
  !>   |   2   |  `1, c+g`   |   `c, g`   |
  !>   |   3   |  `1, c, g`  |     -      |
  !>
  !> The option `w_adj` determines the weights of the graph edges resulting
  !> from element adjacency at the faces, edges and vertices. While the faces
  !> are always taken into account, edges and vertex connections are included
  !> only if the corresponding weights do not vanish.

  type PartitioningOptions_3D
    logical :: child    = .false. !< set T/F to partition child or given level
    logical :: split    = .false. !< use subdivision of parent partitions
    logical :: lazy     = .false. !< adopt parent partitioning if possible
    integer :: n_parts  = 1       !< number of partitions requested
    integer :: n_con    = 1       !< number of constraints (1 .. 3)
    integer :: n_sub    = huge(1) !< max number of sublevels to be weighted
    integer :: c_active = 2       !< cost of active child elements
    integer :: c_frozen = 1       !< cost of frozen child elements
    integer :: w_adj(3) = [1,0,0] !< element face/edge/vertex adjacency weights
  contains
    procedure :: Bcast => Bcast_PartitioningOptions
  end type PartitioningOptions_3D

  integer, parameter :: parent_mode = 1
  integer, parameter :: child_mode  = 2

contains

  !=============================================================================
  ! PartitioningOptions_3D: type-bound procedures

  subroutine Bcast_PartitioningOptions(opt, root, comm)
    class(PartitioningOptions_3D), intent(inout) :: opt !< options
    integer       , intent(in) :: root !< rank root process
    type(MPI_Comm), intent(in) :: comm !< MPI communicator

    call XMPI_Bcast(opt % child    , root, comm)
    call XMPI_Bcast(opt % split    , root, comm)
    call XMPI_Bcast(opt % lazy     , root, comm)
    call XMPI_Bcast(opt % n_parts  , root, comm)
    call XMPI_Bcast(opt % n_con    , root, comm)
    call XMPI_Bcast(opt % n_sub    , root, comm)
    call XMPI_Bcast(opt % c_active , root, comm)
    call XMPI_Bcast(opt % c_frozen , root, comm)
    call XMPI_Bcast(opt % w_adj    , root, comm)

  end subroutine Bcast_PartitioningOptions

  !=============================================================================
  ! Partitioner interfaces

  !-----------------------------------------------------------------------------
  !> ParMETIS interface

  subroutine ParMETIS_Partitioner_3D(opt, mesh, tp_elem, n_parts)
    use ParMETIS_Binding

    ! arguments ................................................................

    class(PartitioningOptions_3D), intent(in) :: opt
    class(Mesh_3D), intent(in) :: mesh
    integer, intent(inout) :: tp_elem(:) !< child element target partitions
    integer, intent(out)   :: n_parts    !< actual number of partitions

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

    character(len=:), allocatable :: prefix
    logical :: logging
    integer :: c(0:3)
    integer :: e, i, j, k, l, m, n

    !---------------------------------------------------------------------------
    ! Body

    logging = mesh%proc == 0 .and. log_level > 0 .or. &
              mesh%proc  > 0 .and. log_level > 1
    prefix  = LoggingPrefix('ParMETIS_Partitioner_3D', mesh%proc)

    if (logging) then
      print '(2A)', prefix, 'start'
    end if

    associate( n_elem   => mesh % n_elem   &
             , n_ghost  => mesh % n_ghost  &
             , c_active => opt % c_active  &
             , c_frozen => opt % c_frozen  &
             , w_adj    => opt % w_adj     )

      ! initialization .........................................................

      ! vertices = contributing elements
      allocate(vtx_elem(n_elem + n_ghost), source = -1)
      n = 0
      do e = 1, n_elem
        if (opt%child) then
          if (mesh%element(e)%adaptation%mark < 100) cycle
        end if
        vtx_elem(e) = n
        n = n + 1
      end do
      nvtx = n







      ! create MPI communicator comprising all parent meshes with nvtx > 0
      if (nvtx > 0) then
        m = 1
      else
        m = 0
      end if
      call MPI_Comm_split(mesh%comm_parts, m, mesh%part, comm)

      call MPI_Comm_rank(comm, proc)
      call MPI_Comm_size(comm, nproc)

      ! current vertex distribution over processes
      allocate(nvtx_proc(0:nproc-1))
      call MPI_Allgather(n, 1, MPI_INTEGER, nvtx_proc, 1, MPI_INTEGER, comm)

      ! graph vertex ID element variable and transfer buffer
      var_vtx_elem(1:1, 1:1, 1:1, 1:size(vtx_elem)) => vtx_elem
      buf_vtx_elem = ElementTransferBuffer_3D(mesh, var_vtx_elem)

      ! ParMetis input arguments ...............................................

      ! use C-style numbering for ParMETIS
      numflag = 0

      ! number of new partitions
      nparts = min(opt%n_parts, sum(nvtx_proc))

      ! number of constraints, 1 ≤ ncon ≤ 3
      ncon = max(min(opt%n_con,3), 1)

      ! weight flag
      if (any(w_adj > 0)) then
        ! graph vertex and edge constraints
        wgtflag = 3
      else
        ! graph vertex constraints only
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

      ! graph vertex offsets ...................................................

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

      m = 0
      do e = 1, n_elem
        if (vtx_elem(e) < 0) cycle
        associate(element => mesh % element(e))

          m = m + 1

          xadj(m) = xadj(m-1)

          ! count graph edges contributed by element faces
          do k = 1, 6
            n = element % face(k) % n_neighbor - 1
            if (n < 0) cycle
            i = element % face(k) % i_neighbor
            do j = i, i+n
              if (vtx_elem(element % neighbor(j) % id) < 0) cycle
              xadj(m) = xadj(m) + 1
            end do
          end do

          if (w_adj(2) > 0) then
            ! count graph edges contributed by element edges
            do k = 1, 12
              n = element % edge(k) % n_neighbor - 1
              if (n < 0) cycle
              i = element % edge(k) % i_neighbor
              do j = i, i+n
                if (vtx_elem(element % neighbor(j) % id) < 0) cycle
                xadj(m) = xadj(m) + 1
              end do
            end do
          end if

          if (w_adj(3) > 0) then
            ! count graph edges contributed by element vertices
            do k = 1, 8
              n = element % vertex(k) % n_neighbor - 1
              if (n < 0) cycle
              i = element % vertex(k) % i_neighbor
              do j = i, i+n
                if (vtx_elem(element % neighbor(j) % id) < 0) cycle
                xadj(m) = xadj(m) + 1
              end do
            end do
          end if

        end associate
      end do

      ! adjacency and adjacency weights ........................................

      m = xadj(nvtx) - 1 ! number of graph edges
      allocate(adjncy(0:m))
      if (any(w_adj > 0)) then
        allocate(adjwgt(0:m), source = 0)
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
              l = vtx_elem(element % neighbor(j) % id)
              if (l < 0) cycle
              adjncy(m) = l
              if (w_adj(1) > 0) then
                adjwgt(m) = w_adj(1)
              end if
              m = m + 1
            end do
          end do

          ! element edges
          if (w_adj(2) > 0) then
            do k = 1, 12
              n = element % edge(k) % n_neighbor - 1
              if (n < 0) cycle
              i = element % edge(k) % i_neighbor
              do j = i, i+n
                l = vtx_elem(element % neighbor(j) % id)
                if (l < 0) cycle
                adjncy(m) = l
                adjwgt(m) = w_adj(2)
                m = m + 1
              end do
            end do
          end if

          ! element vertices
          if (w_adj(3) > 0) then
            do k = 1, 8
              n = element % vertex(k) % n_neighbor - 1
              if (n < 0) cycle
              i = element % vertex(k) % i_neighbor
              do j = i, i+n
                l = vtx_elem(element % neighbor(j) % id)
                if (l < 0) cycle
                adjncy(m) = l
                adjwgt(m) = w_adj(3)
                m = m + 1
              end do
            end do
          end if

        end associate
      end do

      if (logging) then
        print '(2A)', prefix, 'parameters'
        print '(A,2X,99(G0,1X))', prefix, 'wgtflag      =', wgtflag
        print '(A,2X,99(G0,1X))', prefix, 'numflag      =', numflag
        print '(A,2X,99(G0,1X))', prefix, 'ncon         =', ncon
        print '(A,2X,99(G0,1X))', prefix, 'nparts       =', nparts
        print '(A,2X,99(G0,1X))', prefix, 'ubvec        =', ubvec
        print '(A,2X,99(G0,1X))', prefix, 'options      =', options
        print '(A,2X,99(G0,1X))', prefix, 'size(vtxdist)=', size(vtxdist)
        print '(A,2X,99(G0,1X))', prefix, 'size(xadj)   =', size(xadj)
        print '(A,2X,99(G0,1X))', prefix, 'size(adjncy) =', size(adjncy)
        print '(A,2X,99(G0,1X))', prefix, 'size(vwgt)   =', size(vwgt)
        print '(A,2X,99(G0,1X))', prefix, 'size(adjwgt) =', size(adjwgt)
        print '(A,2X,99(G0,1X))', prefix, 'min(vtxdist) =', minval(vtxdist)
        print '(A,2X,99(G0,1X))', prefix, 'max(vtxdist) =', maxval(vtxdist)
        print '(A,2X,99(G0,1X))', prefix, 'min(xadj)    =', minval(xadj)
        print '(A,2X,99(G0,1X))', prefix, 'max(xadj)    =', maxval(xadj)
        print '(A,2X,99(G0,1X))', prefix, 'min(adjncy)  =', minval(adjncy)
        print '(A,2X,99(G0,1X))', prefix, 'max(adjncy)  =', maxval(adjncy)
        print '(A,2X,99(G0,1X))', prefix, 'min(vwgt)    =', minval(vwgt)
        print '(A,2X,99(G0,1X))', prefix, 'max(vwgt)    =', maxval(vwgt)
        print '(A,2X,99(G0,1X))', prefix, 'min(adjwgt)  =', minval(adjwgt)
        print '(A,2X,99(G0,1X))', prefix, 'max(adjwgt)  =', maxval(adjwgt)
        print '(A,2X,99(G0,1X))', prefix, 'min(tpwgts)  =', minval(tpwgts)
        print '(A,2X,99(G0,1X))', prefix, 'max(tpwgts)  =', maxval(tpwgts)
      end if

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

      n_parts = nparts

      ! clean-up ...............................................................

      call MPI_Comm_free(comm)

    end associate

    !---------------------------------------------------------------------------

    if (logging) then
      print '(2A)', prefix, 'exit'
    end if

  end subroutine ParMETIS_Partitioner_3D

  !-----------------------------------------------------------------------------
  !> Identification of graph vertices

  subroutine IdentifyGraphVertices(opt, mesh, vtx_elem, nvtx)
    class(PartitioningOptions_3D), intent(in) :: opt
    class(Mesh_3D), intent(in) :: mesh
    integer, intent(inout) :: vtx_elem(:) !< map from elements to graph vertices
    integer, intent(out)   :: nvtx        !< number of graph vertices

    integer :: e, n

    n = 0
    do e = 1, mesh%n_elem
      if (opt%child .and. mesh%element(e)%adaptation%mark < 100) cycle
      vtx_elem(e) = n
      n = n + 1
    end do
    nvtx = n

  end subroutine IdentifyGraphVertices

  !-----------------------------------------------------------------------------
  !> Determination of graph vertex weights

  subroutine GetGraphVertexWeights(opt, mesh, vtx_elem, vwgt)
    class(PartitioningOptions_3D), intent(in) :: opt
    class(Mesh_3D),       intent(in)  :: mesh
    integer,              intent(in)  :: vtx_elem(:)
    integer(METIS_IDX_T), intent(out) :: vwgt(0:,0:) !< graph vertex weights

    integer :: c_active, c_frozen
    integer :: c, e, g, i, l, w

    n_con = size(vwgt,1)

    if (opt%child .and. n_con > 1 .or. n_con > 2) then
      c_active = opt % c_active
      c_frozen = opt % c_frozen
    else
      c_active = 1
      c_frozen = 1
    end if

    do e = 1, mesh%n_elem
      associate(element => mesh % element(e))

        i = vtx_elem(e)
        if (i < 0) cycle

        ! children
        select case(element % adaptation % mark)
        case(100:108)
          c = 1 * c_frozen  ! 1 frozen child @ vertex or cloned
        case(201:212)
          c = 2 * c_frozen  ! 2 frozen children @ edge
        case(401:406)
          c = 4 * c_frozen  ! 4 frozen children @ face
        case(800)
          c = 8 * c_frozen  ! 8 frozen children @ element
        case(1000,8000)
          c = 8 * c_active  ! 8 active children @ element or cloned
        end select

        ! grandchildren and further descendants
        g = 0
        w = 1
        do l = 2, min(element%adaptation%sublevels, opt%n_sub)
          g = g + w
          w = 8 * w
        end do

        if (opt%children) then
          select case(n_con)
          case(1)
            vwgt(0,i) = c + 64 * g
          case(2)
            vwgt(0,i) = c
            vwgt(1,i) = g
          end select
        else
          vwgt(0,i) = 1
          select case(n_con)
          case(2)
            vwgt(1,i) = c + 64 * g
          case(3)
            vwgt(1,i) = c
            vwgt(2,i) = g
          end select
        end if

      end associate
    end do

  end subroutine GetGraphVertexWeights

  !=============================================================================

end module Partitioner_Interface__3D
