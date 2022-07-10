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
    class(Mesh_3D), target, intent(inout) :: mesh !< mesh partition

    type(MPI_Comm)  :: comm_split
    type(MPI_Group) :: group_world, group_split, group_active
    integer, allocatable :: ranks(:)
    integer :: active
    integer :: i

    associate( comm_world => mesh % comm_world &
             , n_parts    => mesh % n_parts    )

      ! preparations ...........................................................

      ! prepare list of active processes
      allocate(mesh % proc_active(0:n_parts-1))

      ! distinguish between active and inactive processes
      if (mesh % part >= 0) then
        active = 1
      else
        active = 0
      end if

      ! array for selection of group ranks
      allocate(ranks(0:n_parts-1))
      do i = 0, n_parts-1
        ranks(i) = i
      end do

      ! split world communicator into active and inactive processes
      call MPI_Comm_split(comm_world, active, mesh%part, comm_split)

      ! identify associated groups
      call MPI_Comm_group(comm_world, group_world) ! identical for all
      call MPI_Comm_group(comm_split, group_split) ! differs between in/active

      if (active > 0) then

        ! communicator and maps for active processes ...........................

        ! intra-partition communicator
        mesh % comm_parts = comm_split

        ! list of active processes within comm_world
        call MPI_Group_translate_ranks( group_split, n_parts, ranks     &
                                      , group_world, mesh % proc_active )

        mesh % proc_part(0:) => mesh % proc_active

        call MPI_Group_free(group_world)
        call MPI_Group_free(group_split)

      else

        ! communicator and maps for inactive process ...........................

        ! deactivate intra-partition
        mesh % comm_parts = MPI_COMM_NULL

        ! list of active processes within comm_world
        call MPI_Group_difference(group_world, group_split, group_active )
        call MPI_Group_translate_ranks( group_active, n_parts, ranks  &
                                      , group_world, mesh % proc_active )

        mesh % proc_part => null()

        call MPI_Comm_free(comm_split)
        call MPI_Group_free(group_world)
        call MPI_Group_free(group_split)
        call MPI_Group_free(group_active)

      end if

    end associate

  end subroutine BuildCommunicator

  !=============================================================================

end submodule MP_BuildCommunicator
