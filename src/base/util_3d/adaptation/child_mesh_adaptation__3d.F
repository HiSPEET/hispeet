module Child_Mesh_Adaptation__3D
  use XMPI
  use Mesh__3D
  use Element_Transfer_Buffer__3D
  use Element_Distribution_Map__3D
  use Partitioner_Interface__3D

  implicit none
  private

  public :: PartitioningOptions_3D   ! from Partitioner_Interface__3D
  public :: ChildMeshAdaptation_3D

contains

  !=============================================================================
  ! Adaptation routines

  !-----------------------------------------------------------------------------
  !> Computes the target family partitions and first child element IDs

  subroutine BuildFamilyDistributionMap(opt, parent_mesh, family_map)
    class(ChildMeshAdaptationOptions_3D), intent(in) :: opt
    class(Mesh_3D), intent(in) :: parent_mesh
    type(ElementDistributionMap_3D), intent(out) :: family_map

    integer, allocatable, target :: tp_elem(:)
    integer, contiguous, pointer :: tp_elem_val(:,:,:,:)
    type(ElementTransferBuffer_3D), allocatable, asynchronous :: tp_elem_buf

    allocate(tp_elem(parent_mesh%n_elem + parent_mesh%n_ghost), source = -1)

    ! graph-based partitioning
    call ParMETIS_Partitioner(opt, parent_mesh, tp_elem)

  end subroutine BuildFamilyDistributionMap

  !=============================================================================

end module Child_Mesh_Adaptation__3D
