module Child_Mesh_Adaptation
  use XMPI
  use Mesh__3D
  use Element_Transfer_Buffer__3D
  use Element_Distribution_Map__3D

  implicit none
  private

  public :: ChildMeshAdaptationOptions_3D
  public :: ChildMeshAdaptation_3D

  type ChildMeshAdaptationOptions_3D
    integer :: n_parts   = 1         !< number of new partitions
    integer :: method    = 1         !< weighting method
    integer :: weight(4) = [1,0,0,0] !< element and connectivity weights
  contains
    procedure :: Bcast => Bcast_AdaptationOptions
  end type ChildMeshAdaptationOptions_3D

contains

  !=============================================================================
  ! ChildMeshAdaptationOptions_3D: type-bound procedures

  subroutine Bcast_AdaptationOptions(opt, root, comm)
    class(ChildMeshAdaptationOptions_3D), intent(inout) :: opt !< options
    integer       , intent(in) :: root !< rank root process
    type(MPI_Comm), intent(in) :: comm !< MPI communicator

    call XMPI_Bcast(opt % n_parts , root, comm)
    call XMPI_Bcast(opt % method  , root, comm)
    call XMPI_Bcast(opt % weight  , root, comm)

  end subroutine Bcast_AdaptationOptions

  !=============================================================================
  ! Adaptation routines

  !-----------------------------------------------------------------------------
  !> Computes the target family partitions and first child element IDs

  subroutine BuildFamilyDistributionMap(opt, parent_mesh, family_map)
    class(ChildMeshAdaptationOptions_3D), intent(in) :: opt
    class(Mesh_3D), intent(in) :: parent_mesh
    type(ElementDistributionMap_3D), intent(out) :: family_map

    integer, allocatable, target :: tp_elem(:)
    integer, contiguous, pointer :: tp_elem_val(:,:,:,:)
    type(ElementTransferBuffer_3D), allocatable, asynchronous :: tp_elem_buf

    allocate(tp_elem(parent_mesh%n_elem + parent_mesh%n_ghost), source = -1)

    ! graph-based partitioning
    call ParMETIS_Partitioner(opt, parent_mesh, tp_elem)

  end subroutine BuildFamilyDistributionMa

  !-----------------------------------------------------------------------------
  !> Partitioning with ParMETIS

  subroutine ParMETIS_Partitioner(opt, parent_mesh, tp_elem)
    use ParMETIS_Binding

    ! arguments ................................................................

    class(ChildMeshAdaptationOptions_3D), intent(in) :: opt
    class(Mesh_3D), intent(in) :: parent_mesh
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

    integer :: e, i, j, k, m, n

    !---------------------------------------------------------------------------
    ! Body

    associate( n_elem  => parent_mesh % n_elem   &
             , n_ghost => parent_mesh % n_ghost  &
             , weight  => opt % weight           )

      ! initialization .........................................................

      ! vertices = elements with descendants
      allocate(vtx_elem(n_elem + n_ghost), source = -1)
      i = 0
      do e = 1, n_elem
        if (parent_mesh % element(e) % adaptation % mark >= 0) then
          vtx_elem(e) = i
          i = i + 1
        end if
      end do
      nvtx = i

      ! create MPI communicator comprising all parent meshes with nvtx > 0
      if (nvtx > 0) then
        m = 1
      else
        m = 0
      end if
      call MPI_Comm_split(parent_mesh%comm_parts, m, parent_mesh%part, comm)

      call MPI_Comm_rank(comm, proc)
      call MPI_Comm_size(comm, nproc)

      ! ParMetis input arguments ...............................................

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

      ! graph vertex ID element variable and transfer buffer
      var_vtx_elem(1:1, 1:1, 1:1, 1:size(vtx_elem)) => vtx_elem
      buf_vtx_elem = ElementTransferBuffer_3D(parent_mesh, var_vtx_elem)

      if (ncon > 0) then
        ubvec  = 1.05          ! tolerances for multi-constraint weighting
        tpwgts = 1.00 / nparts ! fractions of vertex weight per partition
      end if

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
      end do

      ! transfer graph vertex IDs to ghosts
      call buf_vtx_elem % Transfer(parent_mesh, var_vtx_elem, tag=1000)
      call buf_vtx_elem % Merge(var_vtx_elem)

      ! graph vertex weights and adjacency offsets .............................

      xadj(0) = 0

      i = 0
      do e = 1, n_elem
        if (vtx_elem(e) < 0) cycle
        associate(element => parent_mesh % element(e))

          ! primary weight set proportional to designated number of children
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
          vwgt(0,i) = weight(1) * m

          ! secondary weight according to number of sublevels
          if (ncon >= 2) then
            select case(element % adaptation % sublevels)
            case(2)
              m = 8     ! = 8^1
            case(3)
              m = 64    ! = 8^2
            case(4)
              m = 576   ! = 8^2 + 8^3
            case(5:)
              m = 4672  ! = 8^2 + 8^3 + 8^4
            case default
              m = 0
            end select
            vwgt(1,i) = weight(1) * m
          end if

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
        associate(element => parent_mesh % element(e))

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

  end subroutine ParMETIS_Partitioner

  !=============================================================================

end module Child_Mesh_Adaptation
