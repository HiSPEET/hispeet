module Element_Distribution_Map__3D
  use XMPI
  use Mesh__3D
  use Element_Transfer_Buffer__3D

  implicit none
  private

  public :: ElementDistributionMap_3D

  !-----------------------------------------------------------------------------
  !> Map defining the redistribution of local element data
  !>
  !> The map refers to the elements and ghosts of the `mesh` that was given at
  !> instantiation. Hence
  !>
  !>     n_elem  = mesh % n_elem
  !>     n_ghost = mesh % n_ghost
  !>
  !> whereas `n_parts` is the number of new partitions and generally different
  !> from `mesh%n_parts`. The dimensions and the meaning of the array components
  !> are as follows:
  !>
  !>  – `tp_elem(1:n_elem + n_ghost)` with
  !>       *  `tp_elem(1:n_elem) `    target partition of the local elements
  !>       *  `tp_elem(n_elem+1:)`    target partition of the ghosts' masters
  !>
  !>
  !>  – `id_elem(1:n_elem + n_ghost)` with
  !>       *  `tid_elem(1:n_elem) `   element IDs in target partition
  !>       *  `tid_elem(n_elem+1:)`   ghosts' master IDs in target partition
  !>
  !>  - `ne_part(0:n_parts-1)` number of elements contributed to the new
  !>     partitions, ghosts are not counted
  !>
  !> For elements or ghosts `e` with no target partition are set as follows
  !>
  !>     tp_elem(e) = -1
  !>     id_elem(e) =  0

  type ElementDistributionMap_3D
    integer :: n_parts = 0             !< num target partitions
    integer :: n_elem  = 0             !< num elements
    integer :: n_ghost = 0             !< num ghosts
    integer, allocatable :: tp_elem(:) !< element target partitions
    integer, allocatable :: id_elem(:) !< element IDs in target partitions
    integer, allocatable :: ne_part(:) !< num elements targeted to partition (0:)
  end type ElementDistributionMap_3D

  ! constructor
  interface ElementDistributionMap_3D
    procedure New_Map
  end interface

contains

  !-----------------------------------------------------------------------------
  !> New distribution map

  function New_Map(mesh, n_parts, tp_elem) result(this)

    type(ElementDistributionMap_3D) :: this

    class(Mesh_3D), intent(in) :: mesh
    integer, intent(in) :: n_parts
    integer, intent(in) :: tp_elem(:)

    call BuildMap(this, mesh, n_parts, tp_elem)

  end function New_Map

  !-----------------------------------------------------------------------------
  !> Build distribution map including ghost contributions

  subroutine BuildMap(this, mesh, n_parts, tp_elem)
    class(ElementDistributionMap_3D), target, intent(inout) :: this
    class(Mesh_3D), intent(in) :: mesh
    integer, intent(in) :: n_parts
    integer, intent(in) :: tp_elem(mesh%n_elem + mesh%n_ghost)

    ! internal variables .......................................................

    type(ElementTransferBuffer_3D), allocatable, asynchronous :: id_elem_buf
    integer, contiguous , pointer :: id_elem_val(:,:,:,:)
    integer :: id_part(0:n_parts-1)
    integer :: i, p

    ! basic initialization .....................................................

    this % n_parts = n_parts
    this % n_elem  = mesh % n_elem
    this % n_ghost = mesh % n_ghost
    this % tp_elem = tp_elem
    this % tp_elem = max(this % tp_elem, -1)

    if (allocated(this % id_elem)) deallocate(this % id_elem)
    if (allocated(this % ne_part)) deallocate(this % ne_part)

    allocate(this % id_elem(1 : this%n_elem + this%n_ghost), source = -1)
    allocate(this % ne_part(0 : n_parts-1), source = 0)

    ! number of local elements targeted to partition ...........................

    do i = 1, this % n_elem
      p = tp_elem(i)
      if (p >= 0) then
        this % ne_part(p) = this % ne_part(p) + 1
      end if
    end do

    ! IDs of local elements in their target partitions .........................

    call ComputeElementOffsets(mesh%comm_parts, this%ne_part, id_part)

    do i = 1, this % n_elem
      p = tp_elem(i)
      if (p >= 0) then
        id_part(p) = id_part(p) + 1    ! increment element counter
        this % id_elem(i) = id_part(p) ! element ID in new partition
      end if
    end do

    ! IDs of the ghosts' masters in their target partition .....................

    if (this % n_ghost > 0) then
      id_elem_val(1:1,1:1,1:1,1:this%n_elem+this%n_ghost) => this % id_elem
      id_elem_buf = ElementTransferBuffer_3D(mesh, id_elem_val)
      call id_elem_buf % Transfer(mesh, id_elem_val, tag = 1000)
      call id_elem_buf % Merge(id_elem_val)
    end if

  end subroutine BuildMap

  !-----------------------------------------------------------------------------
  !> Computes the offsets for numbering the redistributed elements

  subroutine ComputeElementOffsets(comm_parts, ne_part, id_part)
    type(MPI_Comm), intent(in) :: comm_parts
    integer, asynchronous, intent(in)  :: ne_part(0:)
    integer, asynchronous, intent(out) :: id_part(0:)

    type(MPI_Win) :: window
    integer(MPI_ADDRESS_KIND) :: integer_extent, lb
    integer(MPI_ADDRESS_KIND) :: buf_size
    integer(MPI_ADDRESS_KIND) :: target_disp = 0
    integer :: disp_unit
    integer :: old_n_parts, old_part
    integer :: i, np

    ! preliminaries ............................................................

    call MPI_Comm_rank(comm_parts, old_part)
    call MPI_Comm_size(comm_parts, old_n_parts)

    np = size(ne_part)

    ! create MPI window ........................................................

    call MPI_Type_get_extent(MPI_INTEGER, lb, integer_extent)
    buf_size  = integer_extent * np
    disp_unit = int(integer_extent)

    call MPI_Win_create( base       =  id_part        &
                       , size       =  buf_size       &
                       , disp_unit  =  disp_unit      &
                       , info       =  MPI_INFO_NULL  &
                       , comm       =  comm_parts &
                       , win        =  window         )

    ! compute offsets via accumulation .........................................

    id_part = 0

    call MPI_Win_fence(MPI_MODE_NOSTORE + MPI_MODE_NOPRECEDE, window)

    do i = old_part + 1, old_n_parts - 1
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

end module Element_Distribution_Map__3D
