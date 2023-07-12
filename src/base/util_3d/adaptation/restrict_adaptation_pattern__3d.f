module Restrict_Adaptation_Pattern__3D
  use Mesh__3D
  use Data_Exchange__3D

  implicit none
  private

  public :: RestrictAdaptationPattern_3D

contains

  !-----------------------------------------------------------------------------
  !> Restrict adaptation pattern from child to parent level
  !>
  !> Upgrades the adaptation mark for all parent elements with active children
  !> to `max(mark, child_mark+1)`, where `child_mark` is the highest mark of all
  !> children.

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
    integer :: n_recv, n_send
    integer :: c, i, o


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

    allocate(child_mark(8, child  % n_cluster), source = -1)
    allocate( send_mark(   child  % n_cluster), source = -1)
    allocate( recv_mark(   parent % n_elem   ), source = -1)

    ! prepare child marks for sending ..........................................

    do i = 1, child % n_elem
      associate(element => child % element(i))
        if (element % frozen) cycle
        c = element % cluster_id
        o = element % cluster_oct
        child_mark(o,c) = element % adaptation % mark
      end associate
    end do

    do c = 1, child % n_cluster
      send_mark(c) = maxval(child_mark(:,c))
      if (send_mark(c) <= 0) cycle
      do i = 0, 7
        if (child_mark(i+1,c) > 0) then
          send_mark(c) = send_mark(c) + IBset(0,i)
        end if
      end do
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
!### CHECK
if (i == 1) then
write(*,'(99(G0,1X))') 'parent % element(1) % adaptation % refinement =', &
                        parent % element(1) % adaptation % refinement
write(*,'(99(G0,1X))') 'parent % element(1) % adaptation % mark       =', &
                        parent % element(1) % adaptation % mark
write(*,'(99(G0,1X))') 'recv_mark       (1)                           =', &
                        recv_mark       (1)
end if
!### CHECK END
        if (adaptation % refinement < 100) cycle
        adaptation % mark = max(adaptation % mark, recv_mark(i) + 1000)
!### CHECK
if (i == 1) then
write(*,'(99(G0,1X))') 'parent % element(1) % adaptation % refinement =', &
                        parent % element(1) % adaptation % refinement
end if
!### CHECK END
      end associate
    end do

  end subroutine RestrictAdaptationPattern_3D

  !=============================================================================

end module Restrict_Adaptation_Pattern__3D
