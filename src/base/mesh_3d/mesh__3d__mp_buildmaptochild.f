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

!> summary:  Generation of map to child elements
!> author:   Joerg Stiller
!> date:     2023/04/11
!===============================================================================

submodule(Mesh__3D) MP_BuildMapToChild
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Generation of map to child elements
  !>
  !> The routine provides
  !>
  !>   - `mesh % n_child`
  !>   - `mesh % map_child(:)`
  !>
  !> and requires `mesh % element`
  !>
  !>   -  `mesh % element % adaptation % refinement`
  !>   -  `mesh % element % adaptation % child_proc`

  module subroutine BuildMapToChild(mesh)
    class(Mesh_3D), intent(inout) :: mesh  !< mesh parition

    integer, allocatable :: n_elem(:)    ! total num elements  per proc
    integer, allocatable :: n_active(:)  ! num active elements per proc
    integer, allocatable :: n_frozen(:)  ! num frozen elements per proc
    integer, allocatable :: map_proc(:)  ! map entry per proc
    integer :: n_proc
    integer :: c, e, p, p_min, p_max

    !$omp master

    if (allocated(mesh % map_child)) then
      deallocate(mesh % map_child)
    end if

    if (mesh % n_elem < 1) then

      ! top or empty mesh ......................................................

      mesh % n_child = 0
      allocate(mesh % map_child(0))

    else

      ! initialization .........................................................

      call MPI_Comm_size(mesh % comm_world, n_proc)

      allocate( n_elem  (0:n_proc-1), source = 0)
      allocate( n_active(0:n_proc-1), source = 0)
      allocate( n_frozen(0:n_proc-1), source = 0)
      allocate( map_proc(0:n_proc-1), source = 0)

      p_min = n_proc-1
      p_max = 0

      ! counts .................................................................

      do e = 1, mesh % n_elem
        p = mesh % element(e) % adaptation % child_proc
        if (p >= 0) then
          p_min = min(p_min, p)
          p_max = max(p_max, p)
          select case(mesh % element(e) % adaptation % refinement)
          case(1000,8000)
            n_active(p) = n_active(p) + 1
          case(100:800)
            n_frozen(p) = n_frozen(p) + 1
          end select
        end if
      end do

      ! setup of maps ..........................................................

      c = 0
      do p = p_min, p_max
        n_elem(p) = n_active(p) + n_frozen(p)
        if (n_elem(p) > 0) then
          c = c + 1
          map_proc(p) = c
        end if
      end do
      mesh % n_child = c
      allocate(mesh % map_child(c))

      c = 0
      do p = p_min, p_max
        if (n_elem(p) > 0) then
          c = c + 1
          mesh % map_child(c) % comm     = mesh % comm_world
          mesh % map_child(c) % proc     = p
          mesh % map_child(c) % n_elem   = n_elem(p)
          mesh % map_child(c) % n_active = n_active(p)
          mesh % map_child(c) % n_frozen = n_frozen(p)
        end if
      end do

      do c = 1, mesh % n_child
        allocate(mesh % map_child(c) % id_elem( mesh % map_child(c) % n_elem ))
      end do

      ! build element lists ....................................................

      n_active = 0
      n_frozen = 0

      associate(map => mesh % map_child)
        do e = 1, mesh % n_elem
          p = mesh % element(e) % adaptation % child_proc
          if (p >= 0) then
            c = map_proc(p)
            select case(mesh % element(e) % adaptation % refinement)
            case(1000,8000)
              n_active(p) = n_active(p) + 1
              map(c) % id_elem(n_active(p)) = e
            case(100:800)
              n_frozen(p) = n_frozen(p) + 1
              map(c) % id_elem(map(c) % n_active + n_frozen(p)) = e
            end select
          end if
        end do
      end associate

    end if

    !$omp end master
    !$omp barrier

  end subroutine BuildMapToChild

  !=============================================================================

end submodule MP_BuildMapToChild
