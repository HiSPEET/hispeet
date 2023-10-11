module Mesh_Map_To_Child__3D
  use XMPI
  implicit none
  private

  public :: MeshMapToChild_3D

  !-----------------------------------------------------------------------------
  !> Map supporting the transfer of data to and from child elements
  !>
  !> The entries in `id_elem` are sorted according 1) the element ID, 2) the
  !> activity of the corresponding child elements. Elements with active children
  !> precede those with inactive ones.
  !>
  !> @note
  !> The communicator `comm` is identical to `comm_world` in the embedding data
  !> structure `mesh`, but is included to allow communication without access to
  !> the latter.

  type MeshMapToChild_3D
    type(MPI_Comm) :: comm !< MPI communicator, usually comm_world
    integer :: proc        !< MPI rank of process keeping the child partition
    integer :: n_elem      !< num elements contributing to child partition
    integer :: n_active    !< num elements with active children
    integer :: n_frozen    !< num elements with frozen children
    integer, allocatable :: id_elem(:) !< IDs of contributing elements
  end type MeshMapToChild_3D

  !=============================================================================

end module Mesh_Map_To_Child__3D
