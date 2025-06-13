!> summary:  Adapt a multilevel mesh
!> author:   Joerg Stiller
!> date:     2024/12/30
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(ML__Mesh__3D) MP_Adapt
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Adapt multilevel mesh

  module subroutine Adapt(this, part_opt, x_plan)
    class(ML_Mesh_3D),                      intent(inout) :: this
    class(MeshPartitionerOptions_3D),       intent(in)    :: part_opt(:)
    type(DataExchangePlan_3D), allocatable, intent(out)   :: x_plan(:)

    type(Mesh_3D), allocatable, save :: old_mesh(:)
    integer, save :: l_top_loc, l_top_new, l_top_old

    type(MeshPartitionerOptions_3D) :: part_opt_next
    integer :: l

    ! initialization ...........................................................

    l_top_old = size(this % mesh)
    do l = l_top_old, 1, -1
      if (this % mesh(l) % n_elem > 0) then
        if (any(this % mesh(l) % element % adaptation % mark > 0)) exit
      end if
    end do
    l_top_loc = l + 1
    call XMPI_Allreduce(l_top_loc, l_top_new, MPI_MAX, this%mesh(1)%comm_world)

    call move_alloc(this%mesh, old_mesh)
    allocate(this%mesh(l_top_new))
    allocate(x_plan(min(l_top_old, l_top_new)))

    ! make adaptation pattern consistent .....................................

    do l = l_top_old, 2, -1
      call GlobalizeAdaptationPattern_3D(old_mesh(l))
      call RestrictAdaptationPattern_3D(old_mesh(l), old_mesh(l-1))
    end do
    call GlobalizeAdaptationPattern_3D(old_mesh(1))

    ! root level partitioning ..................................................

    call ProcessAdaptationPattern_3D(old_mesh(1))

    part_opt_next = part_opt(1)
    part_opt_next % n_con_root = min(4, l_top_new, part_opt_next % n_con_root)

    if (max(old_mesh(1)%n_parts, part_opt(1)%n_parts) == 1) then
      this % mesh(1) = old_mesh(1)
      x_plan(1) % identity = .true.
      allocate(x_plan(1) % send_map(0))
      allocate(x_plan(1) % recv_map(0))
    else if (old_mesh(1) % is_top) then
      call RootMeshPartitioning_3D( opt      = part_opt_next &
                                  , old_mesh = old_mesh(1)   &
                                  , new_mesh = this%mesh(1)  &
                                  , x_plan   = x_plan(1)     )
    else
      call RootMeshPartitioning_3D( opt      = part_opt_next &
                                  , old_mesh = old_mesh(1)   &
                                  , new_mesh = this%mesh(1)  &
                                  , child    = old_mesh(2)   &
                                  , x_plan   = x_plan(1)     )
    end if

    ! adaptation ...............................................................

    do l = 1, l_top_new-1

      this%mesh(l)%refinement = old_mesh(l) % refinement

      if (l > 1) then
        call ProcessAdaptationPattern_3D(this%mesh(l))
      end if

      part_opt_next = part_opt(l+1)
      part_opt_next % n_con_child = min( 3                           &
                                       , l_top_new - l               &
                                       , part_opt_next % n_con_child )

      select case(l_top_old - l)
      case(0)
        call ChildMeshAdaptation_3D( opt        = part_opt_next  &
                                   , parent     = this%mesh(l)   &
                                   , new_child  = this%mesh(l+1) )
      case(1)
        call ChildMeshAdaptation_3D( opt        = part_opt_next  &
                                   , parent     = this%mesh(l)   &
                                   , new_child  = this%mesh(l+1) &
                                   , old_child  = old_mesh(l+1)  &
                                   , x_plan     = x_plan(l+1)    )
      case(2:)
        call ChildMeshAdaptation_3D( opt        = part_opt_next  &
                                   , parent     = this%mesh(l)   &
                                   , new_child  = this%mesh(l+1) &
                                   , old_child  = old_mesh(l+1)  &
                                   , grandchild = old_mesh(l+2)  &
                                   , x_plan     = x_plan(l+1)    )
      end select

    end do

    ! finalization .............................................................

    deallocate(old_mesh)

  end subroutine Adapt

  !=============================================================================

end submodule MP_Adapt
