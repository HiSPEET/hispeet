submodule(Child_Mesh_Adaptation__3D) MP_BuildChildMap
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Computes the target family partitions and first child element IDs

  module subroutine BuildChildMap(opt, parent, map)
    class(PartitioningOptions_3D), intent(in)  :: opt
    class(Mesh_3D),                intent(in)  :: parent
    type(ChildDistributionMap_3D), intent(out) :: map

    integer, allocatable, target :: tp_child(:)
    integer, contiguous, pointer :: tp_child_val(:,:,:,:)
    type(ElementTransferBuffer_3D), allocatable, asynchronous :: tp_child_buf

!### CHECK
print '(99(G0,1X))', 'BCM 0, proc',parent%proc
!### CHECK END
    allocate(tp_child(parent%n_elem + parent%n_ghost), source = -1)
!### CHECK
print '(99(G0,1X))', 'BCM 1, proc',parent%proc,', size(tp_child) =',size(tp_child)
!### CHECK END

    ! graph-based partitioning
    call ParMETIS_Partitioner_3D(opt, parent, tp_child)
!### CHECK
print '(99(G0,1X))', 'BCM 2, proc',parent%proc
!### CHECK END

    ! transfer target partition IDs to ghosts
    if (parent % n_ghost > 0) then
      tp_child_val(1:1,1:1,1:1,1:parent%n_elem+parent%n_ghost) => tp_child
      tp_child_buf = ElementTransferBuffer_3D(parent, tp_child_val)
      call tp_child_buf % Transfer(parent, tp_child_val, tag = 1000)
      call tp_child_buf % Merge(tp_child_val)
    end if
!### CHECK
print '(99(G0,1X))', 'BCM 3, proc',parent%proc
!### CHECK END

    map = ChildDistributionMap_3D(parent, opt%n_parts, tp_child)
!### CHECK
print '(99(G0,1X))', 'BCM X, proc',parent%proc
!### CHECK END

  end subroutine BuildChildMap

  !=============================================================================

end submodule MP_BuildChildMap
