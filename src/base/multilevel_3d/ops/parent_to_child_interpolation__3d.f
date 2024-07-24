module Parent_To_Child_Interpolation__3D
  use Kind_Parameters
  use Mesh__3D
  use Data_Exchange__3D
  implicit none
  private

contains

  subroutine ParentToChildInterpolation_3D(parent, child, iop, v_p, v_c)
    class(Mesh_3D), intent(in) :: parent
    class(Mesh_3D), intent(in) :: child
    class(CoarseToFineInterpolationOperator_1D) :: iop
    real(RNP), contiguous, intent(in)    :: v_p(0:,0:,0:,:,:)
    real(RNP), contiguous, intent(inout) :: v_c(0:,0:,0:,:,:)

    type(DataExchangeMap_3D), allocatable, save :: send_map(:)
    type(DataExchangeMap_3D), allocatable, save :: recv_map(:)

    type(DataExchangeSendBuf_3D), allocatable, save :: send_buf(:)
    type(DataExchangeRecvBuf_3D), allocatable, save :: recv_buf(:)

    real(RNP), allocatable, save :: v_r(:,:,:,:,:)

    integer, save :: n_comp, n_recv, n_send

    integer :: i

    ! initialization ...........................................................

    ! number of components
    n_comp = size(u_p,5)

    ! parent send map and send buffer
    n_send = parent % n_child
    allocate(send_map(n_send), send_buf(n_send))
    do i = 1, n_send
      send_map(i) = DataExchangeMap_3D(parent % map_child(i))
    end do

    ! child receive map and receive buffer
    n_recv = child % n_parent
    allocate(recv_map(n_recv), recv_buf(n_recv))
    do i = 1, n_recv
      recv_map(i) = DataExchangeMap_3D(child % map_parent(i))
    end do

    ! prepare parent data received by child
    allocate(v_r(0:iop%po_c,0:iop%po_c,0:iop%po_c,child%n_cluster,n_comp))

    ! pass parent data to children .............................................

    ! parent extracts data and starts sending
    do i = 1, n_send
      call send_buf(i) % Extract_Data(send_map(i), v = v_p)
      call send_buf(i) % Send_Start()
    end do

    ! child prepares buffer and starts receiving
    do i = 1, n_recv
      call recv_buf(i) % Init(recv_map(i), v = v_r)
      call recv_buf(i) % Recv_Start()
    end do

    ! parent finalizes sending
    do i = 1, n_send
      call send_buf(i) % Send_Finish()
    end do

    ! child finalizes receiving and assigns data to
    do i = 1, n_recv
      call recv_buf(i) % Recv_Finish()
      call recv_buf(i) % Assign_Data(v_r)
    end do

    ! interpolate parent data to child variable ................................

  end subroutine ParentToChildInterpolation_3D

end module Parent_To_Child_Interpolation__3D
