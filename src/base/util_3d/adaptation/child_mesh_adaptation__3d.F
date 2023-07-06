module Child_Mesh_Adaptation__3D
  use Kind_Parameters
  use Constants
  use Execution_Control
  use XMPI
  use Parent_To_Child_Interpolation__1D
  use Mesh__3D
  use Mesh_Element__3D
  use Mesh_Element_Indexing__3D
  use Element_Transfer_Buffer__3D
  use Element_Distribution_Map__3D
  use Element_Distribution_Data__3D
  use Child_Distribution_Map__3D
  use Partitioner_Interface__3D
  use Data_Exchange__3D

  implicit none
  private

  public :: PartitioningOptions_3D   ! from Partitioner_Interface__3D
  public :: ChildMeshAdaptation_3D

  !=============================================================================
  ! module procedures defined in submodules

  interface

    !---------------------------------------------------------------------------
    !> Computes the target family partitions and first child element IDs

    module subroutine BuildChildMap(opt, parent, map)
      class(PartitioningOptions_3D), intent(in)  :: opt
      class(Mesh_3D),                intent(in)  :: parent
      type(ChildDistributionMap_3D), intent(out) :: map
    end subroutine BuildChildMap

    !---------------------------------------------------------------------------
    !> Build child distribution data for sending from parent to child partitions

    module subroutine BuildChildData(parent, child, map, child_data)
      class(Mesh_3D), intent(in) :: parent
      class(Mesh_3D), intent(in) :: child
      class(ChildDistributionMap_3D),    intent(in)  :: map
      class(ElementDistributionData_3D), intent(out) :: child_data(0:)
    end subroutine BuildChildData

    !---------------------------------------------------------------------------
    !> Build connections between parent, new child, old child and grandchildren

    module subroutine BuildConnections( parent, new_child_map, new_child &
                                      , old_child, grandchild, exch_plan )

      class(Mesh_3D),                       intent(inout) :: parent
      class(ChildDistributionMap_3D),       intent(in)    :: new_child_map
      class(Mesh_3D),                       intent(in)    :: new_child
      class(Mesh_3D),             optional, intent(in)    :: old_child
      class(Mesh_3D),             optional, intent(inout) :: grandchild
      class(DataExchangePlan_3D), optional, intent(out)   :: exch_plan

    end subroutine BuildConnections

  end interface

contains

  !=============================================================================
  ! Adaptation routines

  !-----------------------------------------------------------------------------
  !> Generation of an adapted child mesh

  subroutine ChildMeshAdaptation_3D( opt, parent, new_child, old_child &
                                   , grandchild, exch_plan             )

    ! arguments ................................................................

    class(PartitioningOptions_3D),        intent(in)    :: opt
    class(Mesh_3D),                       intent(inout) :: parent
    class(Mesh_3D),                       intent(out)   :: new_child
    class(Mesh_3D),             optional, intent(in)    :: old_child
    class(Mesh_3D),             optional, intent(inout) :: grandchild
    class(DataExchangePlan_3D), optional, intent(out)   :: exch_plan

    ! internal variables .......................................................

    type(MeshAttributes_3D)         , allocatable, save :: attrib
    type(ChildDistributionMap_3D)   , allocatable, save :: new_child_map
    type(ElementDistributionData_3D), allocatable, save :: send_data(:)
    type(ElementDistributionData_3D), allocatable, save :: recv_data(:)

    integer :: i

    ! sanity check .............................................................

    if (present(old_child)) then
      if (.not. present(exch_plan)) then
        call Error( 'ChildMeshAdaptation_3D' &
                  , 'exch_plan is required, when old_child is present' &
                  , 'Child_Mesh_Adaptation__3D' )
      end if
    end if

    ! prepare new partitions ...................................................

    !$omp master

    attrib = MeshAttributes_3D(parent)
    attrib % structured = .false.
    attrib % is_root = .false.
    attrib % dx = HALF * attrib % dx

    call new_child % Init_Mesh_3D(attrib, parent % comm_world)

    new_child % n_parts = opt % n_parts

    ! use identity for mapping parts to processes, so far
    if (new_child % proc < new_child % n_parts) then
      new_child % part = new_child % proc
    end if

    call new_child % BuildCommunicator()

    ! prepare send data
    if (parent % part >= 0) then
      allocate(send_data(0:new_child%n_parts-1))
    end if

    ! prepare receive data
    if (new_child % part >= 0) then
      allocate(recv_data(0:parent%n_parts-1))
      do i = 0, parent%n_parts-1
        recv_data(i) % proc = parent % proc_part(i)
        recv_data(i) % comm = parent % comm_world
      end do
    end if

    !$omp end master

    ! build new child element map ..............................................

    !$omp master
    if (parent % part >= 0) then
      allocate(new_child_map)
      call BuildChildMap(opt, parent, new_child_map)
    end if
    !$omp end master

    ! create child element data ................................................

    if (parent % part >= 0) then
      call BuildChildData(parent, new_child, new_child_map, send_data)
    end if

    ! transfer child element data ................................................

    !$omp master

    ! send dimensions (skip local part)
    if (parent % part >= 0) then
      do i = 0, new_child%n_parts-1
        if (send_data(i) % proc == parent % proc) cycle
        call send_data(i) % Send_Dimensions()
      end do
    end if

    ! receive dimensions (skip local part)
    if (new_child % part >= 0) then
      do i = 0, parent%n_parts-1
        if (recv_data(i) % proc == parent % proc) cycle
        call recv_data(i) % Recv_Dimensions()
      end do
    end if

    ! send data (skip local part)
    if (parent % part >= 0) then
      do i = 0, new_child%n_parts-1
        if (send_data(i) % proc == parent % proc) cycle
        call send_data(i) % Send_Data()
      end do
    end if

    ! receive data (skip local part)
    if (new_child % part >= 0) then
      do i = 0, parent%n_parts-1
        if (recv_data(i) % proc == parent % proc) cycle
        call recv_data(i) % Recv_Data()
      end do
    end if

    ! finalize transfer
    if (parent % part >= 0) then
      do i = 0, new_child%n_parts-1
        if (send_data(i) % proc == parent % proc) cycle
        call send_data(i) % Send_Finish()
      end do
    end if
    if (new_child % part >= 0) then
      do i = 0, parent%n_parts-1
        if (recv_data(i) % proc == parent % proc) cycle
        call recv_data(i) % Recv_Finish()
      end do
    end if

    !$omp end master

    ! assign child element data ................................................

    !$omp master

    ! number of elements
    if (parent % part >= 0) then
      new_child % n_elem = sum(recv_data % n_elem) &
                         + send_data(parent % part) % n_elem
    else
      new_child % n_elem = sum(recv_data % n_elem)
    end if

    allocate(new_child % element(new_child % n_elem))

    !$omp end master

    ! assign received data
    if (new_child % part >= 0) then
      do i = 0, parent%n_parts-1
        if (recv_data(i) % proc /= parent % proc) then
          call recv_data(i) % Assign_Data(new_child)
        else
          call send_data(parent%part) % Assign_Data(new_child)
        end if
      end do
    end if
!### CHECK
write(*,'(99(G0,1X))') 'child mesh adaptation: new_child element(65) parent_id   =', &
  new_child%element(65)%adaptation%parent_id
write(*,'(99(G0,1X))') 'child mesh adaptation: new_child element(65) cluster_oct =', &
  new_child%element(65)%cluster_oct
write(*,'(99(G0,1X))') 'child mesh adaptation: send_data n_elem =', &
  send_data(0) % n_elem
write(*,'(99(G0,1X))') 'child mesh adaptation: send_data element(64) id          =', &
  send_data(0) % element(64)%id
write(*,'(99(G0,1X))') 'child mesh adaptation: send_data element(64) parent_id   =', &
  send_data(0) % element(64)%adaptation%parent_id
write(*,'(99(G0,1X))') 'child mesh adaptation: send_data element(64) cluster_oct =', &
  send_data(0) % element(64)%cluster_oct
write(*,'(99(G0,1X))') 'child mesh adaptation: send_data element(65) id          =', &
  send_data(0) % element(65)%id
write(*,'(99(G0,1X))') 'child mesh adaptation: send_data element(65) parent_id   =', &
  send_data(0) % element(65)%adaptation%parent_id
write(*,'(99(G0,1X))') 'child mesh adaptation: send_data element(65) cluster_oct =', &
  send_data(0) % element(65)%cluster_oct
!### CHECK END

    ! complete mesh.............................................................

    if (new_child % part >= 0) then

      !$omp master
      new_child % n_elem_frozen = count(new_child % element % frozen)
      new_child % n_elem_active = new_child % n_elem - new_child % n_elem_frozen
      !$omp end master

      call new_child % IdentifyVertices()
!### CHECK
!! do i = 1, new_child % n_elem
!!   write(*,'(99(G0,1X))') 'elem:',i,', vid =',new_child%element(i)%vertex%id
!!   if (any(new_child%element(i)%vertex%id < 1)) stop "SOMETHING IS WRONG"
!! end do
!### CHECK END
      call new_child % IdentifyEdges()
!### CHECK
!! do i = 1, new_child % n_elem
!!   write(*,'(99(G0,1X))') 'elem:',i,', eid =',new_child%element(i)%edge%id
!!   if (any(new_child%element(i)%edge%id < 1)) stop "SOMETHING IS WRONG"
!! end do
!### CHECK END
      call new_child % BuildFaces()
      call new_child % BuildBoundaryFaces()
      call new_child % BuildLinks()
      call new_child % BuildGhosts()
      call new_child % BuildCuboids()
      call new_child % IdentifyRanks()
      call new_child % IdentifyClusters()
      call new_child % BuildMapToParent()
!### CHECK
write(*,'(99(G0,1X))') '# CMA  I'
!### CHECK END

    end if

    ! complete connections and save retained data ..............................
    !
    ! - mapping from parent to new children
    !     + `parent % n_child`
    !     + `parent % map_child`
    !
    ! - mapping from grandchildren to new children
    !     + `grandchild % n_parent`
    !     + `grandchild % map_parent`
    !     + `grandchild % element % adaptation % parent_proc`
    !     + `grandchild % element % adaptation % parent_id`
    !
    !  - mapping between old and new children
    !     + `exch_plan` for redistributing retained data

    call BuildConnections( parent, new_child_map, new_child &
                         , old_child, grandchild, exch_plan )

    ! copy retained element attributes .........................................
    !
    ! - `new_child % element % adaptation % refinement`
    ! - `new_child % element % adaptation % child_proc`
    ! - `new_child % element % adaptation % mark`
!### CHECK
write(*,'(99(G0,1X))') '# CMA  J'
!### CHECK END

    if (present(old_child)) then
      call CopyRetainedData(old_child, new_child, exch_plan)
    end if
!### CHECK
write(*,'(99(G0,1X))') '# CMA  K'
!### CHECK END

    ! connect new child with grandchild ........................................
    !
    ! - `new_child % n_child`
    ! - `new_child % map_child`

    if (new_child % part >= 0 .and. present(grandchild)) then
      call new_child % BuildMapToChild()
    end if
!### CHECK
write(*,'(99(G0,1X))') '# CMA  L'
!### CHECK END

    ! finalization .............................................................

    !$omp master

    parent    % is_top  = .false.
    new_child % is_root = .false.

    if (allocated(attrib))         deallocate(attrib)
    if (allocated(new_child_map))  deallocate(new_child_map)
    if (allocated(send_data))      deallocate(send_data)
    if (allocated(recv_data))      deallocate(recv_data)

    !$omp end master
!### CHECK
write(*,'(99(G0,1X))') '# CMA  X'
!### CHECK END

  end subroutine ChildMeshAdaptation_3D

  !-----------------------------------------------------------------------------
  !> Copy retained element attributes

  subroutine CopyRetainedData(old_child, new_child, exch_plan)
    class(Mesh_3D),             intent(in)    :: old_child
    class(Mesh_3D),             intent(inout) :: new_child
    class(DataExchangePlan_3D), intent(in)    :: exch_plan

    type(DataExchangeSendBuf_3D), allocatable :: send_buf(:)
    type(DataExchangeRecvBuf_3D), allocatable :: recv_buf(:)

    integer, allocatable :: send_attrib(:,:)
    integer, allocatable :: recv_attrib(:,:)
    integer :: n_recv, n_send
    integer :: i

    associate( send_map => exch_plan % send_map &
             , recv_map => exch_plan % recv_map )

      ! initialization .........................................................

      n_send = size(send_map)
      n_recv = size(recv_map)

      allocate(send_buf(n_send), send_attrib(3,old_child%n_elem_active))
      allocate(recv_buf(n_recv), recv_attrib(3,new_child%n_elem_active))

      ! prepare send data ......................................................

      do i = 1, old_child % n_elem_active
        send_attrib(1,i) = old_child % element(i) % adaptation % refinement
        send_attrib(2,i) = old_child % element(i) % adaptation % child_proc
        send_attrib(3,i) = old_child % element(i) % adaptation % mark
      end do

      ! pass attributes to new children ........................................

      do i = 1, n_send
        call send_buf(i) % Extract_Data(send_map(i), send_attrib)
        call send_buf(i) % Send_Start()
      end do

      do i = 1, n_recv
        call recv_buf(i) % Init(recv_map(i), recv_attrib)
        call recv_buf(i) % Recv_Start()
      end do

      do i = 1, n_send
        call send_buf(i) % Send_Finish()
      end do

      do i = 1, n_recv
        call recv_buf(i) % Recv_Finish()
        call recv_buf(i) % Assign_Data(recv_attrib)
      end do

      ! extract received attributes ............................................

      do i = 1, new_child % n_elem_active
        associate(element => new_child % element(i))
          element % adaptation % refinement = recv_attrib(1,i)
          element % adaptation % child_proc = recv_attrib(2,i)
          element % adaptation % mark       = recv_attrib(3,i)
        end associate
      end do

    end associate

  end subroutine CopyRetainedData

  !=============================================================================

end module Child_Mesh_Adaptation__3D

