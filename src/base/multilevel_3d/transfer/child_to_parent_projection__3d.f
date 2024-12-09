!> summary:  Projection of mesh data from child to parent level
!> author:   Joerg Stiller
!> date:     2024/08/07
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Child_To_Parent_Projection__3D
  use Kind_Parameters
  use Array_Assignments
  use HP__Coarsening_Operator__1D
  use Mesh__3D
  use Data_Exchange__3D
  use TPO__AAA__3D
  use TPO__8To1__3D
  implicit none
  private

  public :: ChildToParentProjection_3D

  interface ChildToParentProjection_3D
    module procedure ChildToParentProjection_A
    module procedure ChildToParentProjection_S
  end interface

contains

  !-----------------------------------------------------------------------------
  !> 3D child to parent projection of array data

  subroutine ChildToParentProjection_A(child, parent, pop, v_c, v_p)
    class(Mesh_3D), intent(in) :: child
      !< child mesh partition
    class(Mesh_3D), intent(in) :: parent
      !< parent mesh partition
    class(HP_CoarseningOperator_1D), intent(in) :: pop
      !< projection operator
    real(RNP), contiguous, intent(in) :: v_c(0:,0:,0:,:,:)
      !< child data
    real(RNP), contiguous, intent(inout) :: v_p(0:,0:,0:,:,:)
      !< parent data

    call ChildToParentProjection_X(child, parent, pop, size(v_c,5), v_c, v_p)

  end subroutine ChildToParentProjection_A

  !-----------------------------------------------------------------------------
  !> 3D child to parent interpolation of scalar data

  subroutine ChildToParentProjection_S(child, parent, pop, v_c, v_p)
    class(Mesh_3D), intent(in) :: child
      !< child mesh partition
    class(Mesh_3D), intent(in) :: parent
      !< parent mesh partition
    class(HP_CoarseningOperator_1D), intent(in) :: pop
      !< projection operator
    real(RNP), contiguous, intent(in) :: v_c(0:,0:,0:,:)
      !< child data
    real(RNP), contiguous, intent(inout) :: v_p(0:,0:,0:,:)
      !< parent data

    call ChildToParentProjection_X(child, parent, pop, 1, v_c, v_p)

  end subroutine ChildToParentProjection_S

  !-----------------------------------------------------------------------------
  !> 3D child to parent interpolation using arrays of explicit shape

  subroutine ChildToParentProjection_X(child, parent, pop, n_comp, v_c, v_p)
    class(Mesh_3D), intent(in) :: child
      !< child mesh partition
    class(Mesh_3D), intent(in) :: parent
      !< parent mesh partition
    class(HP_CoarseningOperator_1D), intent(in) :: pop
      !< projection operator
    integer,intent(in) :: n_comp
      !< number of array components
    real(RNP), target, intent(in) :: &
      v_c(0:pop%po_f, 0:pop%po_f, 0:pop%po_f, child % n_elem, n_comp)
      !< child data
    real(RNP), intent(inout) :: &
      v_p(0:pop%po_c, 0:pop%po_c, 0:pop%po_c, parent% n_elem, n_comp)
      !< parent data

    type(DataExchangeMap_3D), allocatable, save :: send_map(:)
    type(DataExchangeMap_3D), allocatable, save :: recv_map(:)

    type(DataExchangeSendBuf_3D), allocatable, save :: send_buf(:)
    type(DataExchangeRecvBuf_3D), allocatable, save :: recv_buf(:)

    real(RNP), allocatable, save :: v_s(:,:,:,:,:)

    logical, save :: skip_frozen = .true.
    integer, save :: n_recv, n_send

    integer :: n_cluster_active = 0
    integer :: i, k

    ! initialization ...........................................................

    ! child send map and send buffer
    n_send = child % n_parent
    allocate(send_map(n_send), send_buf(n_send))
    do i = 1, n_send
      send_map(i) = DataExchangeMap_3D(child % map_parent(i), skip_frozen)
    end do

    ! parent receive map and receive buffer
    n_recv = parent % n_child
    allocate(recv_map(n_recv), recv_buf(n_recv))
    do i = 1, n_recv
      recv_map(i) = DataExchangeMap_3D(parent % map_child(i), skip_frozen)
    end do

    ! number of active clusters
    select case(parent%refinement)
    case('c')
      n_cluster_active = child%n_elem_active
    case('s')
      n_cluster_active = child%n_elem_active / 8
    end select

    ! prepare parent data send by child
    allocate(v_s(0:pop%po_c,0:pop%po_c,0:pop%po_c,n_cluster_active,n_comp))

    ! project child data to send variable ......................................

    select case (pop % mode)
    case(0)
      ! identity
      call SetArray( v_s(:,:,:,1:n_cluster_active,:) &
                   , v_c(:,:,:,1:n_cluster_active,:) &
                   , multi = .true.                  )
    case(1)
      ! p-refinement
      do k = 1, n_comp
        call TPO_AAA( pop%A(:,:,1)                    &
                    , v_c(:,:,:,1:n_cluster_active,k) &
                    , v_s(:,:,:,1:n_cluster_active,k) )
      end do
    case(2)
      ! h- or hp-coarsening
      do k = 1, n_comp
        call TPO_8To1( pop%A                             &
                     , v_c(:,:,:,1:n_cluster_active*8,k) &
                     , v_s(:,:,:,1:n_cluster_active  ,k) )
      end do
    end select

    ! pass child data to parents ...............................................

    ! child extracts data and starts sending
    do i = 1, n_send
      call send_buf(i) % Extract_Data(send_map(i), v = v_s)
      call send_buf(i) % Send_Start()
    end do

    ! parent prepares buffer and starts receiving
    do i = 1, n_recv
      call recv_buf(i) % Init(recv_map(i), v = v_p)
      call recv_buf(i) % Recv_Start()
    end do

    ! child finalizes sending
    do i = 1, n_send
      call send_buf(i) % Send_Finish()
    end do

    ! parent finalizes receiving and assigns data
    do i = 1, n_recv
      call recv_buf(i) % Recv_Finish()
      call recv_buf(i) % Assign_Data(v = v_p)
    end do

    ! finalization .............................................................

    deallocate(send_buf, send_map)
    deallocate(recv_buf, recv_map)
    deallocate(v_s)

  end subroutine ChildToParentProjection_X

  !=============================================================================

end module Child_To_Parent_Projection__3D
