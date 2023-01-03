module Child_Mesh_Adaptation__3D
  use Constants
  use XMPI
  use Mesh__3D
  use Element_Transfer_Buffer__3D
  use Element_Distribution_Map__3D
  use Element_Distribution_Data__3D
  use Child_Distribution_Map__3D
  use Partitioner_Interface__3D

  implicit none
  private

  public :: PartitioningOptions_3D   ! from Partitioner_Interface__3D
  public :: ChildMeshAdaptation_3D

contains

  !=============================================================================
  ! Adaptation routines

  !-----------------------------------------------------------------------------
  !> Generation of an adapted child mesh

  subroutine ChildMeshAdaptation_3D(opt, parent, child)
    class(PartitioningOptions_3D), intent(in) :: opt
    class(Mesh_3D), intent(inout) :: parent
    class(Mesh_3D), intent(out)   :: child

    ! internal variables .......................................................

    type(MeshAttributes_3D)       :: attrib
    type(ChildDistributionMap_3D) :: child_map
    integer :: i

    !---------------------------------------------------------------------------
    !$omp master

    ! prepare new partitions ...................................................

    attrib = MeshAttributes_3D(parent)
    attrib % structured = .false.
    attrib % is_root = .false.
    attrib % dx = HALF * attrib % dx

    call child % Init_Mesh_3D(attrib, parent % comm_world)

    child % n_parts = opt % n_parts

    ! use identity for mapping parts to processes, so far
    if (child % proc < child % n_parts) then
      child % part = child % proc
    end if

    call child % BuildCommunicator()

    ! build child element map ..................................................

    if (parent % part >= 0) then
      call BuildChildMap(opt, parent, child_map)
    end if

    !---------------------------------------------------------------------------
    !$omp end master

  end subroutine ChildMeshAdaptation_3D

  !-----------------------------------------------------------------------------
  !> Computes the target family partitions and first child element IDs

  subroutine BuildChildMap(opt, parent, map)
    class(PartitioningOptions_3D), intent(in)  :: opt
    class(Mesh_3D),                intent(in)  :: parent
    type(ChildDistributionMap_3D), intent(out) :: map

    integer, allocatable, target :: tp_child(:)
    integer, contiguous, pointer :: tp_child_val(:,:,:,:)
    type(ElementTransferBuffer_3D), allocatable, asynchronous :: tp_child_buf

    allocate(tp_child(parent%n_elem + parent%n_ghost), source = -1)

    ! graph-based partitioning
    call ParMETIS_Partitioner(opt, parent, tp_child)

    ! transfer target partition IDs to ghosts
    if (parent % n_ghost > 0) then
      tp_child_val(1:1,1:1,1:1,1:parent%n_elem+parent%n_ghost) => tp_child
      tp_child_buf = ElementTransferBuffer_3D(parent, tp_elem_val)
      call tp_child_buf % Transfer(parent, tp_child_val, tag = 1000)
      call tp_child_buf % Merge(tp_child_val)
    end if

    map = ChildDistributionMap_3D(parent, opt%n_parts, tp_child)

  end subroutine BuildChildMap

  !=============================================================================
  ! Generation of child elements

  !-----------------------------------------------------------------------------
  !> Build child distribution data for sending from parent to child partitions

  subroutine BuildDistributionData(parent, child, map, dist_data)

    class(Mesh_3D), intent(in) :: parent
    class(Mesh_3D), intent(in) :: child
    class(ChildDistributionMap_3D),    intent(in)  :: map
    class(ElementDistributionData_3D), intent(out) :: dist_data(0:)

    ! internal variables .......................................................

    integer :: ce_part(0:map%n_parts-1)  ! element counters
    integer :: cn_part(0:map%n_parts-1)  ! neighbor counters
    integer :: cp_part(0:map%n_parts-1)  ! element point counters

    integer :: nc, nn, nn_max, np, tp
    integer :: i, j, k, l

    if (parent % part < 0) return

    ! MPI and dimensions .......................................................

    ! basic initialization
    do tp = 0, map%n_parts-1
      dist_data(tp) % proc       = child % proc_part(tp)
      dist_data(tp) % comm       = child % comm_world
      dist_data(tp) % n_element  = 0
      dist_data(tp) % n_neighbor = 0
      dist_data(tp) % n_point    = 0
    end do

    ! number of elements and upper bounds for number of neighbors and points
    ! per partition
    ce_part = 0
    cn_part = 0
    cp_part = 0
    nn_max  = 0
    do l = 1, parent % n_elem
      associate(element => parent % element(l))
        tp = map % tp_child(l)
        if (tp >= 0) then

          ! number of child elements
          select case(element % adaptation % mark)
          case(1:6)
            nc = 4
          case(7:18)
            nc = 2
          case(19:26)
            nc = 1
          case default
            nc = 8
          end select

          ! max number of neighbors at faces, edges and vertices
          i = max( element % face   % n_neighbor, 1)
          j = max( element % edge   % n_neighbor, 1)
          k = max( element % vertex % n_neighbor, 1)

          ! max number of neighbors per child
          nn = 6 * i + 12 * j + 8 * k
          nn_max = max(nn, nn_max)

          ! number of element points per child
          np = (element % geometry % po + 1) ** 3

          ! add to counters
          ce_part(tp) = ce_part(tp) + nc
          cn_part(tp) = cn_part(tp) + nc * nn
          cp_part(tp) = cp_part(tp) + nc * np

        end if
      end associate
    end do

    ! build components .........................................................

    ! allocate components
    do tp = 0, map%n_parts-1
      allocate( dist_data(tp) % element              ( ce_part(tp)    ) )
      allocate( dist_data(tp) % neighbor_id          ( cn_part(tp)    ) )
      allocate( dist_data(tp) % neighbor_part        ( cn_part(tp)    ) )
      allocate( dist_data(tp) % neighbor_component   ( cn_part(tp)    ) )
      allocate( dist_data(tp) % neighbor_orientation ( cn_part(tp)    ) )
      allocate( dist_data(tp) % x_e                  ( cp_part(tp), 3 ) )
    end do

  end subroutine BuildDistributionData

  !=============================================================================

end module Child_Mesh_Adaptation__3D
