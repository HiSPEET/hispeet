!> summary:  Generation of an intracommunicator between active mesh partitions
!> author:   Joerg Stiller
!> date:     2022/06/12
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(Mesh__3D) MP_BuildCommunicator
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Build communicator between active partitions

  module subroutine BuildCommunicator(mesh)
    class(Mesh_3D), intent(inout) :: mesh !< mesh partition

    type(MPI_Comm)  :: comm_active
    type(MPI_Group) :: group_world
    type(MPI_Group) :: group_parts
    integer, allocatable :: parts(:)
    integer :: active, n_parts
    integer :: i

    ! distinguish between active and inactive processes
    if (mesh % part >= 0) then
      active = 1
    else
      active = 0
    end if

    ! split world communicator into active and inactive processes
    call MPI_Comm_split(mesh%comm_world, active, mesh%proc, comm_active)

    if (active > 0) then

      mesh % comm_parts = comm_active

      ! build map from partition to world communicator
      n_parts = mesh % n_parts
      allocate(mesh % proc_part(0:n_parts-1), parts(0:n_parts-1))
      do i = 0, n_parts-1
        parts(i) = i
      end do
      call MPI_Comm_group(mesh % comm_world, group_world)
      call MPI_Comm_group(mesh % comm_parts, group_parts)
      call MPI_Group_translate_ranks( group_parts, n_parts, parts,  &
                                      group_world, mesh % proc_part )
      call MPI_Group_free(group_world)
      call MPI_Group_free(group_parts)

    else

      mesh % comm_parts = MPI_COMM_NULL
      call MPI_Comm_free(comm_active)

    end if

  end subroutine BuildCommunicator

  !=============================================================================

end submodule MP_BuildCommunicator
