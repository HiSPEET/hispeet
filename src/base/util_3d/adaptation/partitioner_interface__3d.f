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
  !>   - `1` parent
  !>   - `2` child
  !>
  !> In parent mode, the given mesh itself is (re)partitioned and in child mode,
  !> the child mesh is partitioned. This choice determines the constraints that
  !> are considered as follows:
  !>
  !>   | constraint | parent | child |
  !>   | ---------- | ------ | ----- |
  !>   |     1      |  `1`   |  `r`  |
  !>   |     2      |  `r`   |  `c`  |
  !>   |     3      |  `c`   |  `g`  |
  !>
  !> where
  !>
  !>   - `1`  unit workload assigned to all elements
  !>   - `r`  workload derived from planned refinement of given (parent) mesh
  !>   - `c`  workload derived from planned refinement of child mesh
  !>   - `g`  workload derived from planned refinement of grandchild mesh
  !>
  !> For the parent the refinement is evaluated from `element%adaptation%mark`
  !> and for the child and grandchild from `element%adaptation%sublevels`.
  !>
  !> The option `n_const` determines the number of constraints to be considered,
  !> i.e. `1` for the first one only, `2` for the first two, and `3` for all.
  !>
  !> The `w_comp` option assigns relative weights to the graph components:
  !>
  !>   - `w_comp(1) ≥ 1` graph vertices, i.e. mesh elements
  !>   - `w_comp(2) ≥ 0` graph edges emerging from element face connections
  !>   - `w_comp(3) ≥ 0` graph edges emerging from element edge connections
  !>   - `w_comp(4) ≥ 0` graph edges emerging from element vertex connections
  !>
  !> The vertex weight is multiplied with the costs determined by the selected
  !> partitioning mode. Vertices with no cost are excluded from the graph.
  !> Connections that have no weight can be cut without penalty.
  !> Note that the lower bounds of the weights are enforced for regularity.

  type PartitioningOptions_3D
    integer :: mode      = 1         !< partitioning mode
    integer :: n_parts   = 1         !< number of new partitions
    integer :: n_const   = 1         !< number of constraints (1 .. 3)
    integer :: c_active  = 5         !< cost of active elements
    integer :: c_frozen  = 1         !< cost of frozen elements
    integer :: w_comp(4) = [1,0,0,0] !< element and connectivity weights
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

    call XMPI_Bcast(opt % mode      , root, comm)
    call XMPI_Bcast(opt % n_parts   , root, comm)
    call XMPI_Bcast(opt % n_const   , root, comm)
    call XMPI_Bcast(opt % c_active  , root, comm)
    call XMPI_Bcast(opt % c_frozen  , root, comm)
    call XMPI_Bcast(opt % w_comp    , root, comm)

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

    integer :: c(0:3)
    integer :: e, i, j, k, l, m, n
!### CHECK
!! print '(99(G0,1X))', 'PMP 0, proc',mesh%proc
!### CHECK END

    !---------------------------------------------------------------------------
    ! Body

    associate( n_elem   => mesh % n_elem   &
             , n_ghost  => mesh % n_ghost  &
             , c_active => opt % c_active  &
             , c_frozen => opt % c_frozen  &
             , w_comp   => opt % w_comp    )

      ! initialization .........................................................

      ! vertices = contributing elements
      allocate(vtx_elem(n_elem + n_ghost), source = -1)
      i = 0
      do e = 1, n_elem
        if (opt%mode == child_mode) then
          if (mesh%element(e)%adaptation%mark < 100) cycle
        end if
        vtx_elem(e) = i
        i = i + 1
      end do
      nvtx = i
!### CHECK
!! print '(99(G0,1X))', 'PMP 1, proc',mesh%proc
!! print '(99(G0,1X))', 'PMP 1, proc',mesh%proc,'w_comp =', w_comp
!! print '(99(G0,1X))', 'PMP 1, proc',mesh%proc,'min/max mark =', &
!! minval(mesh%element%adaptation%mark),maxval(mesh%element%adaptation%mark)
!### CHECK END

      ! create MPI communicator comprising all parent meshes with nvtx > 0
      if (nvtx > 0) then
        m = 1
      else
        m = 0
      end if
      call MPI_Comm_split(mesh%comm_parts, m, mesh%part, comm)

      call MPI_Comm_rank(comm, proc)
      call MPI_Comm_size(comm, nproc)
!### CHECK
!! print '(99(G0,1X))', 'PMP 2, proc',mesh%proc
!! print '(99(G0,1X))', 'PMP 2, proc',mesh%proc,'m =',m
!! print '(99(G0,1X))', 'PMP 2, proc',mesh%proc,'ParMetis proc  =',proc
!! print '(99(G0,1X))', 'PMP 2, proc',mesh%proc,'ParMetis nproc =',nproc
!### CHECK END

      ! graph vertex ID element variable and transfer buffer
      var_vtx_elem(1:1, 1:1, 1:1, 1:size(vtx_elem)) => vtx_elem
      buf_vtx_elem = ElementTransferBuffer_3D(mesh, var_vtx_elem)
!### CHECK
!! print '(99(G0,1X))', 'PMP 3, proc',mesh%proc
!### CHECK END

      ! ParMetis input arguments ...............................................

      ! use C-style numbering for ParMETIS
      numflag = 0

      ! number of new partitions
      nparts = opt % n_parts

      ! number of constraints, 1 ≤ ncon ≤ 3
      ncon = max(min(opt % n_const,3), 1)

      ! w_comp flag
      if (any(w_comp(1:3) > 0)) then
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
!### CHECK
!! print '(99(G0,1X))', 'PMP 4, proc',mesh%proc
!### CHECK END

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
!### CHECK
!! print '(99(G0,1X))', 'PMP 5, proc',mesh%proc
!### CHECK END

      ! graph vertex IDs by element index (local+ghost) ........................

      ! graph vertex ID of local parent elements
      where(vtx_elem >= 0)
        vtx_elem = vtx_elem + vtxdist(proc)
      end where
!### CHECK
!! print '(99(G0,1X))', 'PMP 6, proc',mesh%proc
!### CHECK END

      ! transfer graph vertex IDs to ghosts
      call buf_vtx_elem % Transfer(mesh, var_vtx_elem, tag=1000)
      call buf_vtx_elem % Merge(var_vtx_elem)
!### CHECK
!! print '(99(G0,1X))', 'PMP 7, proc',mesh%proc
!### CHECK END

      ! graph vertex weights and adjacency offsets .............................

      xadj(0) = 0

      m = 0
      do e = 1, n_elem
        if (vtx_elem(e) < 0) cycle
        associate(element => mesh % element(e))

          ! initialize costs
          c = 0

          ! cost associated with parent
          c(0) = 1

          ! cost associated with children
          select case(element % adaptation % mark)
          case(100:108)
            c(1) = 1 * c_frozen  ! 1 frozen child @ vertex or cloned
          case(201:212)
            c(1) = 2 * c_frozen  ! 2 frozen children @ edge
          case(401:406)
            c(1) = 4 * c_frozen  ! 4 frozen children @ face
          case(800)
            c(1) = 8 * c_frozen  ! 8 frozen children @ element
          case(1000)
            c(1) = 1 * c_active  ! 1 active child    @ element
          case(8000)
            c(1) = 8 * c_active  ! 8 active children @ element
          end select

          ! cost associated with grandchildren and great-grandchildren
          select case(element % adaptation % sublevels)
          case(2)
            c(2) = 1
          case(3:)
            c(2) = 1
            c(3) = 1
          end select

          ! vertex weights
          select case(opt % mode)
          case(parent_mode)
            vwgt(0:ncon-1,m) = c(0:ncon-1)
          case(child_mode)
            vwgt(0:ncon-1,m) = c(1:ncon)
          end select

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

          if (w_comp(3) > 0) then
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

          if (w_comp(4) > 0) then
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
!### CHECK
!! print '(99(G0,1X))', 'PMP 8, proc',mesh%proc
!### CHECK END

      ! adjacency and adjacency weights ........................................

      m = xadj(nvtx) - 1 ! number of graph edges
      allocate(adjncy(0:m))
      if (any(w_comp(1:3) > 0)) then
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
              if (w_comp(2) > 0) then
                adjwgt(m) = w_comp(2)
              end if
              m = m + 1
            end do
          end do

          ! element edges
          if (w_comp(3) > 0) then
            do k = 1, 12
              n = element % edge(k) % n_neighbor - 1
              if (n < 0) cycle
              i = element % edge(k) % i_neighbor
              do j = i, i+n
                l = vtx_elem(element % neighbor(j) % id)
                if (l < 0) cycle
                adjncy(m) = l
                adjwgt(m) = w_comp(3)
                m = m + 1
              end do
            end do
          end if

          ! element vertices
          if (w_comp(4) > 0) then
            do k = 1, 8
              n = element % vertex(k) % n_neighbor - 1
              if (n < 0) cycle
              i = element % vertex(k) % i_neighbor
              do j = i, i+n
                l = vtx_elem(element % neighbor(j) % id)
                if (l < 0) cycle
                adjncy(m) = l
                adjwgt(m) = w_comp(4)
                m = m + 1
              end do
            end do
          end if

        end associate
      end do
!### CHECK
!! print '(99(G0,1X))', 'PMP 9, proc',mesh%proc
!! print '(99(G0,1X))', 'PMP 9, proc',mesh%proc,', wgtflag          =',wgtflag
!! print '(99(G0,1X))', 'PMP 9, proc',mesh%proc,', numflag          =',numflag
!! print '(99(G0,1X))', 'PMP 9, proc',mesh%proc,', ncon             =',ncon
!! print '(99(G0,1X))', 'PMP 9, proc',mesh%proc,', nparts           =',nparts
!! print '(99(G0,1X))', 'PMP 9, proc',mesh%proc,', ubvec            =',ubvec
!! print '(99(G0,1X))', 'PMP 9, proc',mesh%proc,', options          =',options
!! print '(99(G0,1X))', 'PMP 9, proc',mesh%proc,', size(vtxdist)    =',size(vtxdist)
!! print '(99(G0,1X))', 'PMP 9, proc',mesh%proc,', size(xadj)       =',size(xadj)
!! print '(99(G0,1X))', 'PMP 9, proc',mesh%proc,', size(adjncy)     =',size(adjncy)
!! print '(99(G0,1X))', 'PMP 9, proc',mesh%proc,', size(vwgt)       =',size(vwgt)
!! print '(99(G0,1X))', 'PMP 9, proc',mesh%proc,', size(adjwgt)     =',size(adjwgt)
!! print '(99(G0,1X))', 'PMP 9, proc',mesh%proc,', min/max(vtxdist) =',minval(vtxdist),maxval(vtxdist)
!! print '(99(G0,1X))', 'PMP 9, proc',mesh%proc,', min/max(xadj)    =',minval(xadj),maxval(xadj)
!! print '(99(G0,1X))', 'PMP 9, proc',mesh%proc,', min/max(adjncy)  =',minval(adjncy),maxval(adjncy)
!! print '(99(G0,1X))', 'PMP 9, proc',mesh%proc,', min/max(vwgt)    =',minval(vwgt),maxval(vwgt)
!! print '(99(G0,1X))', 'PMP 9, proc',mesh%proc,', min/max(adjwgt)  =',minval(adjwgt),maxval(adjwgt)
!! print '(99(G0,1X))', 'PMP 9, proc',mesh%proc,', min/max(tpwgts)  =',minval(tpwgts),maxval(tpwgts)
!### CHECK END

      ! ParMETIS ...............................................................

      call ParMETIS_V3_PartKway( vtxdist, xadj, adjncy, vwgt, adjwgt, wgtflag  &
                               , numflag, ncon, nparts, tpwgts, ubvec, options &
                               , edgecut, part, comm % MPI_VAL                 )

      ! result .................................................................
!### CHECK
!! print '(99(G0,1X))', 'PMP 10, proc',mesh%proc
!### CHECK END

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
!### CHECK
!! print '(99(G0,1X))', 'PMP X, proc',mesh%proc
!### CHECK END

  end subroutine ParMETIS_Partitioner_3D

  !=============================================================================

end module Partitioner_Interface__3D

