module Restrict_Adaptation_Pattern
  use Mesh__3D
  use Data_Exchange__3D

  implicit none
  private

  public :: RestrictAdaptationPattern

contains

  !-----------------------------------------------------------------------------
  !> Restrict adaptation pattern from child to parent level
  !>
  !> Upgrades the adaptation mark for all parent elements with active children
  !> to `max(mark, child_mark+1)`, where `child_mark` is the highest mark of all
  !> children.

  subroutine RestrictAdaptationPattern(child, parent)
    class(Mesh_3D), intent(in)    :: child
    class(Mesh_3D), intent(inout) :: parent

    type(DataExchangeMap_3D), allocatable :: send_map(:)
    type(DataExchangeMap_3D), allocatable :: recv_map(:)

    type(DataExchangeSendBuf_3D), allocatable :: send_buf(:)
    type(DataExchangeRecvBuf_3D), allocatable :: recv_buf(:)

    integer, allocatable :: send_mark(:)
    integer, allocatable :: recv_mark(:)
    integer :: n_recv, n_send
    integer :: i


    ! initialization ...........................................................

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

    allocate(send_mark(child  % n_cluster), source = -1)
    allocate(recv_mark(parent % n_elem   ), source = -1)

    ! prepare child marks for sending ..........................................

    do i = 1, child % n_elem
      associate(element => child % element(i))
        send_mark(element%cluster_id) = max( send_mark(element%cluster_id) &
                                           , element % adaptation % mark   )
      end associate
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
        if (adaptation % refinement < 100) cycle
        adaptation % mark = max(adaptation % mark, recv_mark(i) + 1)
      end associate
    end do

  end subroutine RestrictAdaptationPattern

end module Restrict_Adaptation_Pattern
