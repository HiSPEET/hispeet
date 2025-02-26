module Restrict_Adaptation_Pattern__3D
  use Execution_Control
  use Mesh__3D
  use Data_Exchange__3D

  implicit none
  private

  public :: RestrictAdaptationPattern_3D

  logical :: ibset_is_safe = .false.

contains

  !-----------------------------------------------------------------------------
  !> Restrict adaptation pattern from child to parent level
  !>
  !> Upgrades the adaptation mark for all parent elements with active children
  !> to `max(mark, child_mark+1000)`, where `child_mark` is the highest mark of
  !> all children.
  !>
  !> Additionally, the value `2^(i-1)` is added to `mark` for each child `i`
  !> that exists and is retained. These indicators are later used to produce a
  !> consistent adaptation pattern.
  !>
  !> The current implementation uses the intrinsics `IBset` and `BTest` to set
  !> and query the indicators, which requires that `IBset(0,i) == 2^i` is
  !> fulfilled for `i < 8`. To verify this, a compatibility check is carried on
  !> the first call.

  subroutine RestrictAdaptationPattern_3D(child, parent)
    class(Mesh_3D), intent(in)    :: child
    class(Mesh_3D), intent(inout) :: parent

    type(DataExchangeMap_3D), allocatable :: send_map(:)
    type(DataExchangeMap_3D), allocatable :: recv_map(:)

    type(DataExchangeSendBuf_3D), allocatable :: send_buf(:)
    type(DataExchangeRecvBuf_3D), allocatable :: recv_buf(:)

    integer, allocatable :: child_mark(:,:)
    integer, allocatable :: send_mark(:)
    integer, allocatable :: recv_mark(:)
    integer :: n_recv, n_send, n_cpe
    integer :: c, i, o

    ! compatibility check ......................................................

    if (.not. ibset_is_safe) then
      do i = 0, 7
        if (IBset(0,i) /= 2**i) then
          call Error('RestrictAdaptationPattern_3D', 'IBset is incompatible')
        end if
      end do
      ibset_is_safe = .true.
    end if

    ! initialization ...........................................................

    ! number of children per element with regular refinement
    select case(parent % refinement)
    case('c')
      n_cpe = 1  ! cloning
    case('s')
      n_cpe = 8  ! subdividing
    case default
      call Error('RestrictAdaptationPattern_3D','parent%refinement not set')
    end select

    n_send = child % n_parent
    allocate(send_map(n_send), send_buf(n_send))
    do i = 1, n_send
      send_map(i) = DataExchangeMap_3D(child % map_parent(i))
    end do

    n_recv = parent % n_child
    allocate(recv_map(n_recv), recv_buf(n_recv))
    do i = 1, n_recv
      recv_map(i) = DataExchangeMap_3D(parent % map_child(i))
    end do

    allocate(child_mark(n_cpe, child  % n_cluster), source = -1)
    allocate( send_mark(       child  % n_cluster), source = -1)
    allocate( recv_mark(       parent % n_elem   ), source = -1)

    ! prepare child marks for sending ..........................................

    do i = 1, child % n_elem_active
      associate(element => child % element(i))
        c = element % cluster_id
        o = element % cluster_oct
        child_mark(o,c) = element % adaptation % mark
      end associate
    end do

    do c = 1, child % n_cluster
      send_mark(c) = maxval(child_mark(:,c))
      if (send_mark(c) <= 0) cycle
      if (n_cpe == 8) then
        ! subdividing: set bit 0:7 if corresponding child is generated
        do i = 1, 8
          if (child_mark(i,c) > 0) then
            send_mark(c) = send_mark(c) + IBset(0,i-1)
          end if
        end do
      else
        ! cloning: set bits 0:7
        do i = 1, 8
          send_mark(c) = send_mark(c) + IBset(0,i-1)
        end do
      end if
    end do

    ! pass child marks to parent  ..............................................

    do i = 1, n_send
      call send_buf(i) % Extract_Data(send_map(i), send_mark)
      call send_buf(i) % Send_Start()
    end do

    do i = 1, n_recv
      call recv_buf(i) % Init(recv_map(i), recv_mark)
      call recv_buf(i) % Recv_Start()
    end do

    do i = 1, n_send
      call send_buf(i) % Send_Finish()
    end do

    do i = 1, n_recv
      call recv_buf(i) % Recv_Finish()
      call recv_buf(i) % Assign_Data(recv_mark)
    end do

    ! upgrade parent marks .....................................................

    do i = 1, parent % n_elem
      associate(adaptation => parent % element(i) % adaptation)
        if (adaptation % refinement < 1000) then
          cycle ! skip elements with no active children
        end if
        adaptation % mark = max(adaptation % mark, recv_mark(i) + 1000)
      end associate
    end do

  end subroutine RestrictAdaptationPattern_3D

  !=============================================================================

end module Restrict_Adaptation_Pattern__3D
