module Mesh_Map_To_Parent__3D
  use XMPI
  implicit none
  private

  public :: MeshMapToParent_3D

  !-----------------------------------------------------------------------------
  !> Map supporting the transfer of data to and from parent elements
  !>
  !> The IDs of sibling element clusters whose parent is kept by process `proc`
  !> are stored in `id_cluster`. The entries are sorted according to 1) the
  !> parent element ID and 2) the element activity. Active clusters precede
  !> inactive ones.
  !>
  !> @note
  !> The communicator `comm` is identical to `comm_world` in the embedding data
  !> structure `mesh`, but is included to allow communication without access to
  !> the latter.

  type MeshMapToParent_3D
    type(MPI_Comm) :: comm !< MPI communicator, usually comm_world
    integer :: proc        !< MPI rank of process keeping the parent partition
    integer :: n_cluster   !< number of clusters originating from `proc`
    integer :: n_active    !< number of active clusters
    integer :: n_frozen    !< number of frozen clusters
    integer, allocatable :: id_cluster(:) !< IDs of related child clusters
  end type MeshMapToParent_3D

  !=============================================================================

end module Mesh_Map_To_Parent__3D
