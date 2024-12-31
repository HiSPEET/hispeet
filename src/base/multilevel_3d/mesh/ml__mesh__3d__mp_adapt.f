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

  module subroutine Adapt(this, partition, exch_plan)
    class(ML_Mesh_3D),                      intent(inout) :: this
    class(PartitioningOptions_3D),          intent(in)    :: partition(:)
    type(DataExchangePlan_3D), allocatable, intent(out)   :: exch_plan(:)

    type(Mesh_3D), allocatable, save :: old_mesh(:)
    integer, save :: l_top_loc, l_top_new, l_top_old

    integer :: e, l

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
    allocate(exch_plan(min(l_top_old, l_top_new)))

    ! make adaptation pattern consistent .....................................

    do l = l_top_new-1, 2, -1
      call GlobalizeAdaptationPattern_3D(old_mesh(l))
      call RestrictAdaptationPattern_3D(old_mesh(l), old_mesh(l-1))
    end do
    call GlobalizeAdaptationPattern_3D(old_mesh(1))

    ! root level partitioning ..................................................

    call ProcessAdaptationPattern_3D(old_mesh(1))

    if (max(old_mesh(1)%n_parts, partition(1)%n_parts) == 1) then
      this % mesh(1) = old_mesh(1)
    else if (old_mesh(1) % is_top) then
      call RootMeshPartitioning_3D( opt       = partition(1) &
                                  , old_mesh  = old_mesh(1)  &
                                  , new_mesh  = this%mesh(1) &
                                  , exch_plan = exch_plan(1) )
    else
      call RootMeshPartitioning_3D( opt       = partition(1) &
                                  , old_mesh  = old_mesh(1)  &
                                  , new_mesh  = this%mesh(1) &
                                  , child     = old_mesh(2)  &
                                  , exch_plan = exch_plan(1) )
    end if

    ! adaptation ...............................................................

    do l = 1, l_top_new-1

      this%mesh(l)%refinement = old_mesh(l) % refinement

      if (l > 1) then
        call ProcessAdaptationPattern_3D(this%mesh(l))
      end if

      select case(l_top_old - l)
      case(0)
        call ChildMeshAdaptation_3D( opt        = partition(l+1) &
                                   , parent     = this%mesh(l)   &
                                   , new_child  = this%mesh(l+1) )
      case(1)
        call ChildMeshAdaptation_3D( opt        = partition(l+1) &
                                   , parent     = this%mesh(l)   &
                                   , new_child  = this%mesh(l+1) &
                                   , old_child  = old_mesh(l+1)  &
                                   , exch_plan  = exch_plan(l+1) )
      case(2:)
        call ChildMeshAdaptation_3D( opt        = partition(l+1) &
                                   , parent     = this%mesh(l)   &
                                   , new_child  = this%mesh(l+1) &
                                   , old_child  = old_mesh(l+1)  &
                                   , grandchild = old_mesh(l+2)  &
                                   , exch_plan  = exch_plan(l+1) )
      end select

    end do

    ! finalization .............................................................

    deallocate(old_mesh)

  end subroutine Adapt

  !=============================================================================

end submodule MP_Adapt
