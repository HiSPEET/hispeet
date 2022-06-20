module Partition_Root_Mesh
  use XMPI
  use Mesh__3D
  use Element_Transfer_Buffer__3D

contains

  ! module
  subroutine ComputeTargetPartitions(mesh, n_parts, weight, method, tp)
    class(Mesh_3D), intent(in) :: mesh  !< local mesh partition
    integer, intent(in)  :: n_parts     !< requested number of partitions
    integer, intent(in)  :: weight(4)   !< element and connectivity weights
    integer, intent(in)  :: method      !< element weighting method
    integer, intent(out) :: tp(mesh%n_elem) !< element target partitions

    ! so far
    call ParMETIS_Partitioner(mesh, n_parts, weight, method, tp)

  end subroutine ComputeTargetPartitions

  !-----------------------------------------------------------------------------
  !> Partitioning with ParMETIS
  !>
  !> @note
  !>   Runs only on the master thread with OpenMP

  subroutine ParMETIS_Partitioner(mesh, n_parts, weight, method, tp)
    use ParMETIS_Binding

    ! arguments ................................................................

    class(Mesh_3D), intent(in) :: mesh  !< local mesh partition
    integer, intent(in)  :: n_parts     !< requested number of partitions
    integer, intent(in)  :: weight(4)   !< element and connectivity weights
    integer, intent(in)  :: method      !< element weighting method
    integer, intent(out) :: tp(mesh%n_elem) !< element target partitions

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
    integer, pointer :: var_vtxid(:,:,:,:) ! maps vtxid to element variable

    type(ElementTransferBuffer_3D), asynchronous :: buf_vtxid

    integer :: e, i, j, k, m

    !$omp master

    associate( comm   => mesh % comm_parts &
             , i_proc => mesh % part       & ! i_proc = present partition ID
             , n_proc => mesh % n_parts    & ! n_proc = present num partitions
             , n_elem => mesh % n_elem     )

      ! initialization .........................................................

      ! use C-style numbering for ParMETIS
      numflag = 0

      ! set default to unweighted graph
      ncon    = 0
      wgtflag = 0

      ! element weighting
      if (weight(1) > 0) then
        select case(method)
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

      tp = part

    end associate

    !$omp end master

  end subroutine ParMETIS_Partitioner

end module Partition_Root_Mesh
