module Child_Distribution_Map__3D
  use XMPI
  use Mesh__3D
  use Element_Transfer_Buffer__3D

  implicit none
  private

  public :: ChildDistributionMap_3D

  !-----------------------------------------------------------------------------
  !> Map defining the distribution of child element data
  !>
  !> The map refers to the elements and ghosts of the `parent` mesh that was
  !> given at instantiation. Hence
  !>
  !>     n_elem  = parent % n_elem
  !>     n_ghost = parent % n_ghost
  !>
  !> `n_parts` is the number of new child partitions and generally differs from
  !> `parent%n_parts`. The dimensions and the meaning of the array components
  !> are as follows:
  !>
  !>  – `mark(1:n_elem + n_ghost)`
  !>     refinement marks of local elements and ghosts
  !>
  !>  – `tp_child(1:n_elem + n_ghost)`
  !>     child target partitions of local elements and ghosts, or `-1` if none
  !>
  !>  – `id_child(2,2,2, 1:n_elem + n_ghost)`
  !>     child IDs in the target partition or `0` if none
  !>       +  in the case of subdividing, id_child(i,j,k,*) is the ID of the
  !>          child at position i,j,k
  !>       +  in the case of cloning, all entries ar set to the child ID
  !>
  !>  – `nc_part(2, 0:n_parts -1)`
  !>     number of active (1) and frozen (2) children contributed to partitions

  type ChildDistributionMap_3D
    integer :: n_parts = 0                    !< num target partitions
    integer :: n_elem  = 0                    !< num elements
    integer :: n_ghost = 0                    !< num ghosts
    integer, allocatable :: mark(:)           !< parent refinement marks
    integer, allocatable :: tp_child(:)       !< child target partitions
    integer, allocatable :: id_child(:,:,:,:) !< child IDs in target partitions
    integer, allocatable :: nc_part(:,:)      !< num children per partition
  end type ChildDistributionMap_3D

  ! constructor
  interface ChildDistributionMap_3D
    procedure New_Map
  end interface

contains

  !-----------------------------------------------------------------------------
  !> New distribution map

  function New_Map(parent, n_parts, tp_child) result(this)

    type(ChildDistributionMap_3D) :: this

    class(Mesh_3D), intent(in) :: parent
    integer, intent(in) :: n_parts
    integer, intent(in) :: tp_child(:)

    call BuildMap(this, parent, n_parts, tp_child)

  end function New_Map

  !-----------------------------------------------------------------------------
  !> Build distribution map including ghost contributions

  subroutine BuildMap(this, parent, n_parts, tp_child)
    class(ChildDistributionMap_3D), target, intent(inout) :: this
    class(Mesh_3D), intent(in) :: parent
    integer, intent(in) :: n_parts
    integer, intent(in) :: tp_child(parent%n_elem + parent%n_ghost)

    ! internal variables .......................................................

    integer, parameter :: &
        ic_regular (2,2,2) = reshape( [1,2,3,4,5,6,7,8], [2,2,2] ), &
        ic_face    (2,2)   = reshape( [1,2,3,4],         [2,2]   ), &
        ic_edge    (2)     = reshape( [1,2],             [2]     )

    type(ElementTransferBuffer_3D), allocatable, asynchronous :: mark_buf
    type(ElementTransferBuffer_3D), allocatable, asynchronous :: id_child_buf

    integer, allocatable :: mark(:,:,:,:), id_child(:,:,:,:), nc_part(:,:)

    integer :: oc_part(2,0:n_parts-1)
    integer :: i, m, p

    ! basic initialization .....................................................

    this % n_parts  = n_parts
    this % n_elem   = parent % n_elem
    this % n_ghost  = parent % n_ghost
    this % tp_child = tp_child
    this % tp_child = max(this % tp_child, -1)

    allocate(mark    (1,1,1, 1 : this%n_elem + this%n_ghost), source = 0)
    allocate(id_child(2,2,2, 1 : this%n_elem + this%n_ghost), source = 0)
    allocate(nc_part (2    , 0 : n_parts - 1               ), source = 0)

    ! number of active and frozen children per target partition ................
    do i = 1, this % n_elem
      mark(1,1,1,i) = parent % element(i) % adaptation % mark
      p = tp_child(i)
      if (p >= 0) then
        select case(parent % element(i) % adaptation % mark)
        case(100:108)
          ! frozen child from vertex refinement or cloning
          nc_part(2,p) = nc_part(2,p) + 1
        case(201:212)
          ! frozen children from edge refinement
          nc_part(2,p) = nc_part(2,p) + 2
        case(401:406)
          ! frozen children from face refinement
          nc_part(2,p) = nc_part(2,p) + 4
        case(800)
          ! frozen children from volume refinement
          nc_part(2,p) = nc_part(2,p) + 8
        case(1000)
          ! active child from cloning
          nc_part(1,p) = nc_part(1,p) + 1
        case(8000)
          ! active children from regular refinement
          nc_part(1,p) = nc_part(1,p) + 8
        end select
      end if
    end do

    ! IDs of local elements in their target partitions .........................

    call ComputeChildOffsets(parent%comm_parts, nc_part, oc_part)

    do i = 1, this % n_elem
      p = tp_child(i)
      if (p >= 0) then

        m = parent % element(i) % adaptation % mark

        if (m < 1000) then

          select case(m)

          ! clone
          case(100); id_child(:,:,:,i) = oc_part(2,p) + 1

          ! vertex
          case(101); id_child(1,1,1,i) = oc_part(2,p) + 1
          case(102); id_child(2,1,1,i) = oc_part(2,p) + 1
          case(103); id_child(1,2,1,i) = oc_part(2,p) + 1
          case(104); id_child(2,2,1,i) = oc_part(2,p) + 1
          case(105); id_child(1,1,2,i) = oc_part(2,p) + 1
          case(106); id_child(2,1,2,i) = oc_part(2,p) + 1
          case(107); id_child(1,2,2,i) = oc_part(2,p) + 1
          case(108); id_child(2,2,2,i) = oc_part(2,p) + 1

          ! edge
          case(201); id_child(:,1,1,i) = oc_part(2,p) + ic_edge
          case(202); id_child(:,2,1,i) = oc_part(2,p) + ic_edge
          case(203); id_child(:,1,2,i) = oc_part(2,p) + ic_edge
          case(204); id_child(:,2,2,i) = oc_part(2,p) + ic_edge
          case(205); id_child(1,:,1,i) = oc_part(2,p) + ic_edge
          case(206); id_child(2,:,1,i) = oc_part(2,p) + ic_edge
          case(207); id_child(1,:,2,i) = oc_part(2,p) + ic_edge
          case(208); id_child(2,:,2,i) = oc_part(2,p) + ic_edge
          case(209); id_child(1,1,:,i) = oc_part(2,p) + ic_edge
          case(210); id_child(2,1,:,i) = oc_part(2,p) + ic_edge
          case(211); id_child(1,2,:,i) = oc_part(2,p) + ic_edge
          case(212); id_child(2,2,:,i) = oc_part(2,p) + ic_edge

          ! face
          case(401); id_child(1,:,:,i) = oc_part(2,p) + ic_face
          case(402); id_child(2,:,:,i) = oc_part(2,p) + ic_face
          case(403); id_child(:,1,:,i) = oc_part(2,p) + ic_face
          case(404); id_child(:,2,:,i) = oc_part(2,p) + ic_face
          case(405); id_child(:,:,1,i) = oc_part(2,p) + ic_face
          case(406); id_child(:,:,2,i) = oc_part(2,p) + ic_face

          ! volume
          case(800); id_child (:,:,:,i) = oc_part(2,p) + ic_regular

          end select

          oc_part(2,p) = oc_part(2,p) + m / 100

        else

          select case(m)
          case(1000)
            id_child(:,:,:,i) = oc_part(1,p) + 1
            oc_part(1,p)      = oc_part(1,p) + 1
          case(8000)
            id_child(:,:,:,i) = oc_part(1,p) + ic_regular
            oc_part(1,p)      = oc_part(1,p) + 8
          end select

        end if
      end if
    end do

    ! IDs of the ghosts' masters in their target partition .....................

    if (this % n_ghost > 0) then
      mark_buf     = ElementTransferBuffer_3D(parent, mark)
      id_child_buf = ElementTransferBuffer_3D(parent, id_child)
      call mark_buf     % Transfer(parent, mark    , tag = 1001)
      call id_child_buf % Transfer(parent, id_child, tag = 1002)
      call id_child_buf % Merge(mark)
      call id_child_buf % Merge(id_child)
    end if

    allocate(this % mark(this%n_elem + this%n_ghost))
    do i = 1, size(this%mark)
      this % mark(i) = mark(1,1,1,i)
    end do

    call move_alloc(id_child, this % id_child)
    call move_alloc(nc_part , this % nc_part )

  end subroutine BuildMap

  !-----------------------------------------------------------------------------
  !> Computes the offsets for numbering the redistributed elements
  !>
  !> On input, `nc_part(m,n)` is the number of local children of class `m`
  !> contributed to the new child partition `n`. Currently `m=1` refers to
  !> active children emerging from regular refinement and `m=2` to frozen
  !> children which form a buffer around refinement zones.
  !> As elements are numbered sequentially class by class, an offset is needed
  !> for each class. For parent partition `p`, the offset of children of class
  !> `m` contributed to child partition `n` equals the sum of corresponding
  !> children contributed by parents `q < p`, i.e.
  !>
  !>       oc_part[p](1,n) = sum(q < p) nc_part[q](1,n)
  !>       oc_part[p](2,n) = sum(q < p) nc_part[q](2,n) + sum(q) nc_part[q](1,n)
  !>
  !> The partial sum is evaluated using one-sided communication based on MPI's
  !> window facility.

  subroutine ComputeChildOffsets(comm_parts, nc_part, oc_part)
    type(MPI_Comm), intent(in) :: comm_parts
    integer, asynchronous, intent(in)  :: nc_part(:,0:)
    integer, asynchronous, intent(out) :: oc_part(:,0:)

    type(MPI_Win) :: window
    integer(MPI_ADDRESS_KIND) :: integer_extent, lb
    integer(MPI_ADDRESS_KIND) :: buf_size
    integer(MPI_ADDRESS_KIND) :: target_disp = 0

    integer, allocatable :: na_loc(:), na_tot(:)

    integer :: disp_unit
    integer :: proc, n_proc
    integer :: i, n

    ! preliminaries ............................................................

    call MPI_Comm_rank(comm_parts, proc)
    call MPI_Comm_size(comm_parts, n_proc)

    n = size(nc_part)

    ! create MPI window ........................................................

    call MPI_Type_get_extent(MPI_INTEGER, lb, integer_extent)
    buf_size  = integer_extent * n
    disp_unit = int(integer_extent)

    call MPI_Win_create( base       =  oc_part        &
                       , size       =  buf_size       &
                       , disp_unit  =  disp_unit      &
                       , info       =  MPI_INFO_NULL  &
                       , comm       =  comm_parts     &
                       , win        =  window         )

    ! compute partial sums via accumulation ....................................

    oc_part = 0

    call MPI_Win_fence(MPI_MODE_NOSTORE + MPI_MODE_NOPRECEDE, window)

    ! add local nc_part to oc_part in processes of higher rank
    do i = proc + 1, n_proc - 1
      call MPI_Accumulate( origin_addr     = nc_part       &
                         , origin_count    = n             &
                         , origin_datatype = MPI_INTEGER   &
                         , target_rank     = i             &
                         , target_disp     = target_disp   &
                         , target_count    = n             &
                         , target_datatype = MPI_INTEGER   &
                         , op              = MPI_SUM       &
                         , win             = window        )
    end do

    call MPI_Win_fence( MPI_MODE_NOSTORE + MPI_MODE_NOPUT + MPI_MODE_NOSUCCEED &
                      , window )

    call MPI_Win_free(window)

    ! adjust offset of frozen elements .........................................

    allocate(na_loc(0:size(nc_part,2)-1), source = nc_part(1,:))
    allocate(na_tot, mold = na_loc)
    call XMPI_Allreduce(na_loc, na_tot, MPI_SUM, comm_parts)

    oc_part(2,:) = oc_part(2,:) + na_tot

  end subroutine ComputeChildOffsets

  !=============================================================================

end module Child_Distribution_Map__3D
