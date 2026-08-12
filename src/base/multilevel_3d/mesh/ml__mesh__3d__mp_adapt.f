!------------------------------------------------------------------------------!
! This file is part of HiSPEET: High-order Spectral Element Techniques         !
!                                                                              !
! Copyright (C) 2026 by the HiSPEET authors and the                            !
! Chair of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany                 !
!                                                                              !
! HiSPEET is free software: you can redistribute it and/or modify              !
! it under the terms of the GNU General Public License as published by         !
! the Free Software Foundation, either version 3 of the License, or            !
! (at your option) any later version.                                          !
!                                                                              !
! HiSPEET is distributed in the hope that it will be useful,                   !
! but WITHOUT ANY WARRANTY; without even the implied warranty of               !
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                         !
! See the GNU General Public License for more details.                         !
!                                                                              !
! You should have received a copy of the GNU General Public License            !
! along with HiSPEET. If not, see <http://www.gnu.org/licenses/>.              !
!------------------------------------------------------------------------------!

!> summary:  Adapt a multilevel mesh
!> author:   Joerg Stiller
!> date:     2024/12/30
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
    character(len=:), allocatable :: prefix
    integer :: l, m

    if (log_level > 0) then
      prefix = LoggingPrefix('ML__Mesh__3D::Adapt', this%mesh(1)%proc)
    end if

    ! initialization ...........................................................

    l_top_old = size(this % mesh)
    do l = l_top_old, 1, -1
      if (this % mesh(l) % n_elem > 0) then
        m = maxval(this % mesh(l) % element % adaptation % mark)
        ! m > 0: at least one element is refined
        ! m = 0: no element is refined, but at least one is retained
        if (m >= 0) exit
      else
        m = -1
      end if
    end do

    l_top_loc = l + min(1,m)
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

    if (max(old_mesh(1)%n_parts, part_opt_next%n_parts) == 1) then
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

      if (log_level > 0) then
        ! check SFC integrity
        if (this%mesh(l) % has_sfc .and. this%mesh(l+1) % n_elem > 0) then
          block
            integer :: max_rk, min_rk
            max_rk = maxval(this % mesh(l+1) % element % sfc_rank)
            min_rk = minval(this % mesh(l+1) % element % sfc_rank)
            if (max_rk - min_rk + 1 /= this % mesh(l+1) % n_elem) then
              print '(9G0)', prefix,'level  = ', l+1
              print '(9G0)', prefix,'n_elem = ', this % mesh(l+1) % n_elem
              print '(9G0)', prefix,'max_rk = ', max_rk
              print '(9G0)', prefix,'min_rk = ', min_rk
              call Error('Adapt','SFC broken','ML__Mesh__3D')
            end if
          end block
        end if
      end if

    end do

    ! finalization .............................................................

    deallocate(old_mesh)

  end subroutine Adapt

  !=============================================================================

end submodule MP_Adapt
