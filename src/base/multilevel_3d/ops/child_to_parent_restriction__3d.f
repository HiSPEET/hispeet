!> summary:  Restriction of mesh data from child to parent level
!> author:   Joerg Stiller
!> date:     2024/09/03
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Child_To_Parent_Restriction__3D
  use Kind_Parameters
  use Constants
  use Array_Assignments
  use HP__Refinement_Operator__1D
  use Mesh__3D
  use Data_Exchange__3D
  use TPO__AAA__3D
  use TPO__1To8__3D
  implicit none
  private

  public :: ChildToParentRestriction_3D

  interface ChildToParentRestriction_3D
    module procedure ChildToParentRestriction_A
    module procedure ChildToParentRestriction_S
  end interface

contains

  !-----------------------------------------------------------------------------
  !> 3D child to parent restriction of array data

  subroutine ChildToParentRestriction_A(child, parent, iop, v_c, v_p)
    class(Mesh_3D), intent(in) :: child
      !< child mesh partition
    class(Mesh_3D), intent(in) :: parent
      !< parent mesh partition
    class(HP_RefinementOperator_1D), intent(in) :: iop
      !< parent to child interpolation operator
    real(RNP), contiguous, intent(in) :: v_c(0:,0:,0:,:,:)
      !< child data
    real(RNP), contiguous, intent(inout) :: v_p(0:,0:,0:,:,:)
      !< parent data

    call ChildToParentRestriction_X(child, parent, iop, size(v_p,5), v_c, v_p)

  end subroutine ChildToParentRestriction_A

  !-----------------------------------------------------------------------------
  !> 3D child to parent restriction of scalar data

  subroutine ChildToParentRestriction_S(child, parent, iop, v_c, v_p)
    class(Mesh_3D), intent(in) :: child
      !< child mesh partition
    class(Mesh_3D), intent(in) :: parent
      !< parent mesh partition
    class(HP_RefinementOperator_1D), intent(in) :: iop
      !< parent to child interpolation operator
    real(RNP), contiguous, intent(in) :: v_c(0:,0:,0:,:)
      !< child data
    real(RNP), contiguous, intent(inout) :: v_p(0:,0:,0:,:)
      !< parent data

    call ChildToParentRestriction_X(child, parent, iop, 1, v_c, v_p)

  end subroutine ChildToParentRestriction_S

  !-----------------------------------------------------------------------------
  !> 3D child to parent restriction using arrays of explicit shape

  subroutine ChildToParentRestriction_X(child, parent, iop, n_comp, v_c, v_p)
    class(Mesh_3D), intent(in) :: child
      !< child mesh partition
    class(Mesh_3D), intent(in) :: parent
      !< parent mesh partition
    class(HP_RefinementOperator_1D), intent(in) :: iop
      !< parent to child interpolation operator
    integer, intent(in) :: n_comp
      !< number of array components
    real(RNP), target, intent(in) :: &
      v_c(0:iop%po_f, 0:iop%po_f, 0:iop%po_f, child % n_elem, n_comp)
      !< child data
    real(RNP), intent(inout) :: &
      v_p(0:iop%po_c, 0:iop%po_c, 0:iop%po_c, parent% n_elem, n_comp)
      !< parent data

    type(DataExchangeMap_3D), allocatable, save :: send_map(:)
    type(DataExchangeMap_3D), allocatable, save :: recv_map(:)

    type(DataExchangeSendBuf_3D), allocatable, save :: send_buf(:)
    type(DataExchangeRecvBuf_3D), allocatable, save :: recv_buf(:)

    real(RNP), allocatable, save :: v_s(:,:,:,:,:)

    logical, save :: skip_frozen = .true.
    integer, save :: n_recv, n_send

    real(RNP), allocatable :: A(:,:,:), B(:,:)
    integer :: n_cluster_active = 0
    integer :: smooth = 0
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
    allocate(v_s(0:iop%po_c,0:iop%po_c,0:iop%po_c,n_cluster_active,n_comp))

    ! transposed interpolation operators
    select case(iop % mode)
    case(1)
      allocate(A(0:iop%po_c, 0:iop%po_f, 1))
      A(0:,0:,1) = transpose(iop % A(:,:,1))
    case(2)
      allocate(A(0:iop%po_c, 0:iop%po_f, 2))
      A(0:,0:,1) = transpose(iop % A(:,:,1))
      A(0:,0:,2) = transpose(iop % A(:,:,2))
      allocate(B(0:iop%po_f,2), source = ZERO)
    end select

    ! restrict child data to send variable .....................................

    select case (iop % mode)
    case(0)
      ! identity
      call SetArray( v_s(:,:,:,1:n_cluster_active,:) &
                   , v_c(:,:,:,1:n_cluster_active,:) &
                   , multi = .true.                  )
    case(1)
      ! p-refinement
      do k = 1, n_comp
        call TPO_AAA( A(:,:,1)                        &
                    , v_c(:,:,:,1:n_cluster_active,k) &
                    , v_s(:,:,:,1:n_cluster_active,k) )
      end do
    case(2)
      ! h- or hp-coarsening
      do k = 1, n_comp
        call TPO_8To1( A, B, smooth                      &
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

  end subroutine ChildToParentRestriction_X

  !=============================================================================

end module Child_To_Parent_Restriction__3D
