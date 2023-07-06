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
  !>  – `tp_child(1:n_elem + n_ghost)`
  !>     child target partitions of local elements and ghosts, or `-1` if none
  !>
  !>  – `id_child(2,2,2, 1:n_elem + n_ghost)`
  !>     child IDs in the target partition or `0` if none
  !>
  !>  – `nc_part(2, 0:n_parts -1)`
  !>     number of active (1) and frozen (2) children contributed to partitions

  type ChildDistributionMap_3D
    integer :: n_parts = 0                    !< num target partitions
    integer :: n_elem  = 0                    !< num elements
    integer :: n_ghost = 0                    !< num ghosts
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

    type(ElementTransferBuffer_3D), allocatable, asynchronous :: id_child_buf

    integer, parameter :: &
        ic_regular (2,2,2) = reshape( [1,2,3,4,5,6,7,8], [2,2,2] ), &
        ic_face    (2,2)   = reshape( [1,2,3,4],         [2,2]   ), &
        ic_edge    (2)     = reshape( [1,2],             [2]     )

    integer, allocatable :: id_child(:,:,:,:), nc_part(:,:)

    integer :: oc_part(2,0:n_parts-1)
    integer :: i, m, p

    ! basic initialization .....................................................

    this % n_parts  = n_parts
    this % n_elem   = parent % n_elem
    this % n_ghost  = parent % n_ghost
    this % tp_child = tp_child
    this % tp_child = max(this % tp_child, -1)

    allocate(id_child(2,2,2, 1 : this%n_elem + this%n_ghost), source = 0)
    allocate(nc_part (2    , 0 : n_parts - 1               ), source = 0)

    ! number of active and frozen children per target partition ................
    do i = 1, this % n_elem
      p = tp_child(i)
      if (p >= 0) then
        select case(parent % element(i) % adaptation % mark)
        case(100)
          ! active children from regular refinement
          nc_part(1,p) = nc_part(1,p) + 8
        case(50)
          ! frozen children from regular refinement
          nc_part(2,p) = nc_part(2,p) + 8
        case(1:6)
          ! frozen children from face refinement
          nc_part(2,p) = nc_part(2,p) + 4
        case(7:18)
          ! frozen children from edge refinement
          nc_part(2,p) = nc_part(2,p) + 2
        case(19:26)
          ! frozen child from vertex refinement
          nc_part(2,p) = nc_part(2,p) + 1
        end select
      end if
    end do

    ! IDs of local elements in their target partitions .........................

    call ComputeChildOffsets(parent%comm_parts, nc_part, oc_part)

    do i = 1, this % n_elem
      p = tp_child(i)
      if (p >= 0) then

        m = parent % element(i) % adaptation % mark

        if (m == 100) then

          ! active children from regular refinement
          id_child (:,:,:,i) = oc_part(1,p) + ic_regular
          oc_part  (1    ,p) = oc_part(1,p) + 8

        else if (m == 50) then

          ! frozen children from regular refinement
          id_child (:,:,:,i) = oc_part(2,p) + ic_regular
          oc_part  (2    ,p) = oc_part(2,p) + 8

        else if (m <= 6) then

          ! frozen children from refinement of face m
          select case(m)
          case(1); id_child(1,:,:,i) = oc_part(2,p) + ic_face
          case(2); id_child(2,:,:,i) = oc_part(2,p) + ic_face
          case(3); id_child(:,1,:,i) = oc_part(2,p) + ic_face
          case(4); id_child(:,2,:,i) = oc_part(2,p) + ic_face
          case(5); id_child(:,:,1,i) = oc_part(2,p) + ic_face
          case(6); id_child(:,:,2,i) = oc_part(2,p) + ic_face
          end select
          oc_part(2,p) = oc_part(2,p) + 4

        else if (m <= 18) then

          ! frozen children from edge refinement
          select case(m)
          case( 7); id_child(:,1,1,i) = oc_part(2,p) + ic_edge
          case( 8); id_child(:,2,1,i) = oc_part(2,p) + ic_edge
          case( 9); id_child(:,1,2,i) = oc_part(2,p) + ic_edge
          case(10); id_child(:,2,2,i) = oc_part(2,p) + ic_edge
          case(11); id_child(1,:,1,i) = oc_part(2,p) + ic_edge
          case(12); id_child(2,:,1,i) = oc_part(2,p) + ic_edge
          case(13); id_child(1,:,2,i) = oc_part(2,p) + ic_edge
          case(14); id_child(2,:,2,i) = oc_part(2,p) + ic_edge
          case(15); id_child(1,1,:,i) = oc_part(2,p) + ic_edge
          case(16); id_child(2,1,:,i) = oc_part(2,p) + ic_edge
          case(17); id_child(1,2,:,i) = oc_part(2,p) + ic_edge
          case(18); id_child(2,2,:,i) = oc_part(2,p) + ic_edge
          end select
          oc_part(2,p) = oc_part(2,p) + 2

        else if (m <= 26) then

          ! frozen child from vertex refinement at vertex m - 18
          select case(m)
          case(19); id_child(1,1,1,i) = oc_part(2,p) + 1
          case(20); id_child(2,1,1,i) = oc_part(2,p) + 1
          case(21); id_child(1,2,1,i) = oc_part(2,p) + 1
          case(22); id_child(2,2,1,i) = oc_part(2,p) + 1
          case(23); id_child(1,1,2,i) = oc_part(2,p) + 1
          case(24); id_child(2,1,2,i) = oc_part(2,p) + 1
          case(25); id_child(1,2,2,i) = oc_part(2,p) + 1
          case(26); id_child(2,2,2,i) = oc_part(2,p) + 1
          end select
          oc_part(2,p) = oc_part(2,p) + 1

        end if
      end if
    end do

    ! IDs of the ghosts' masters in their target partition .....................

    if (this % n_ghost > 0) then
      id_child_buf = ElementTransferBuffer_3D(parent, id_child)
      call id_child_buf % Transfer(parent, id_child, tag = 1000)
      call id_child_buf % Merge(id_child)
    end if

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
  !> This sum is evaluated using one-sided communication based on MPI's
  !> window facility.

  subroutine ComputeChildOffsets(comm_parts, nc_part, oc_part)
    type(MPI_Comm), intent(in) :: comm_parts
    integer, asynchronous, intent(in)  :: nc_part(:,0:)
    integer, asynchronous, intent(out) :: oc_part(:,0:)

    type(MPI_Win) :: window
    integer(MPI_ADDRESS_KIND) :: integer_extent, lb
    integer(MPI_ADDRESS_KIND) :: buf_size
    integer(MPI_ADDRESS_KIND) :: target_disp = 0
    integer :: disp_unit
    integer :: old_n_parts, old_part
    integer :: i, n

    ! preliminaries ............................................................

    call MPI_Comm_rank(comm_parts, old_part)
    call MPI_Comm_size(comm_parts, old_n_parts)

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

    ! compute offsets via accumulation .........................................

    oc_part = 0

    call MPI_Win_fence(MPI_MODE_NOSTORE + MPI_MODE_NOPRECEDE, window)

    do i = old_part + 1, old_n_parts - 1
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

    oc_part(2,:) = oc_part(2,:) + sum(nc_part(1,:))

  end subroutine ComputeChildOffsets

  !=============================================================================

end module Child_Distribution_Map__3D
