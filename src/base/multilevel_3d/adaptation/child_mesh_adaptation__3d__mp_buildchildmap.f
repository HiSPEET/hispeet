submodule(Child_Mesh_Adaptation__3D) MP_BuildChildMap
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Computes the target family partitions and first child element IDs

  module subroutine BuildChildMap(opt, parent, map)
    class(MeshPartitionerOptions_3D), intent(in)  :: opt
    class(Mesh_3D),                   intent(in)  :: parent
    type(ChildDistributionMap_3D),    intent(out) :: map

    integer, allocatable, target, save :: tp_child(:)
    integer, contiguous, pointer, save :: tp_child_val(:,:,:,:) => null()
    type(ElementTransferBuffer_3D), allocatable, asynchronous, save :: tp_child_buf

    character(len=:), allocatable :: prefix
    logical :: logging
    integer :: n_parts

    logging = parent%proc == 0 .and. log_level > 0 .or. &
              parent%proc  > 0 .and. log_level > 1
    prefix  = LoggingPrefix('BuildChildMap', parent%proc)

    if (logging) then
      print '(2A)', prefix, 'start'
    end if

    allocate(tp_child(parent%n_elem + parent%n_ghost), source = -1)

    ! graph-based partitioning
    call MeshPartitioner_3D(opt, parent, tp_child, n_parts)

    ! transfer target partition IDs to ghosts
    if (parent % n_ghost > 0) then
      tp_child_val(1:1,1:1,1:1,1:parent%n_elem+parent%n_ghost) => tp_child
      tp_child_buf = ElementTransferBuffer_3D(parent, tp_child_val)
      call tp_child_buf % Transfer(parent, tp_child_val, tag = 1000)
      call tp_child_buf % Merge(tp_child_val)
    end if

    map = ChildDistributionMap_3D(parent, n_parts, tp_child)

    if (allocated(tp_child))     deallocate(tp_child)
    if (allocated(tp_child_buf)) deallocate(tp_child_buf)
    tp_child_val => null()

    if (logging) then
      print '(2A)', prefix, 'exit'
    end if

  end subroutine BuildChildMap

  !=============================================================================

end submodule MP_BuildChildMap
