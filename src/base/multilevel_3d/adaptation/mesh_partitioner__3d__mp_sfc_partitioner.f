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

!> summary:  Mesh partitioning using space filling curves
!> author:   Joerg Stiller
!> date:     2026/01/27
!===============================================================================

submodule(Mesh_Partitioner__3D) MP_SFC_Partitioner
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Partitioner based on space filling curves

  module subroutine SFC_Partitioner(opt, mesh, tp_elem, n_parts)
    class(MeshPartitionerOptions_3D), intent(in) :: opt
    class(Mesh_3D), intent(in) :: mesh !< current mesh partition
    integer, intent(out) :: tp_elem(:) !< new/child element target partitions
    integer, intent(out) :: n_parts    !< actual number of partitions

    ! internal variables .......................................................

    integer, allocatable, save :: map(:)      ! element ID per SFC rank
    integer, allocatable, save :: wgt(:)      ! accumulated weight per SFC rank
    integer, allocatable, save :: wgt_part(:) ! accumulated partition weights
    integer, allocatable, save :: wgt_dist(:) ! weight distribution

    integer, save :: nn, nn_loc

    character(len=:), allocatable :: prefix
    logical   :: logging
    integer   :: c_active, c_frozen
    integer   :: rk_min, rk_max, sfc_rank_0
    integer   :: wgt_0, wgt_sum
    real(RNP) :: wgt_avg
    integer   :: e, k, m, p

    logging = mesh%proc == 0 .and. log_level > 0 .or. &
              mesh%proc  > 0 .and. log_level > 1
    prefix  = LoggingPrefix('SFC_Partitioner', mesh%proc)

    associate(ne => mesh % n_elem)

      ! initialization .........................................................

      if (logging) then
        print '(2A)', prefix, 'initialization'
      end if

      c_active = max(opt%c_active, 1)
      c_frozen = max(opt%c_frozen, 0)

      allocate(map(ne), wgt(ne), source = 0)
      allocate(wgt_part(0:mesh%n_parts-1), source = 0)

      ! map local SFC rank to element ..........................................

      if (logging) then
        print '(2A)', prefix, 'map local SFC rank to element'
      end if

      if (ne > 0) then
        rk_min = minval(mesh % element % sfc_rank)
        rk_max = maxval(mesh % element % sfc_rank)
        if (rk_max - rk_min + 1 /= ne) then
          if (log_level > 0) then
            print '(9(G0,X))', prefix, 'proc   =',mesh%proc
            print '(9(G0,X))', prefix, 'n_elem =',mesh%n_elem
            print '(9(G0,X))', prefix, 'rk_min =',rk_min
            print '(9(G0,X))', prefix, 'rk_max =',rk_max
          end if
          call Error('SFC_Partitioner','SFC fragmented','Mesh_Partitioner__3D')
        end if
      else
        rk_min = 0
        rk_max = 0
      end if

      ! offset for converting global into local SFC ranks
      sfc_rank_0 = 1 - rk_min

      do e = 1, ne
        map(mesh%element(e)%sfc_rank + sfc_rank_0) = e
      end do

      if (opt % child) then

        ! weights for child partitioning .......................................

        if (logging) then
          print '(2A)', prefix, 'weights for child partitioning'
        end if

        ! individual weights and number of local contributions
        nn_loc = 0
        do e = 1, ne
          k = mesh % element(e) % sfc_rank + sfc_rank_0 ! local rank
          m = mesh % element(e) % adaptation % mark     ! adaptation mark
          if (m >= 1000) then
            wgt(k) = m/1000 * c_active
            nn_loc = nn_loc + 1
          else if (m >= 100 .and. c_frozen > 0) then
            wgt(k) = m/100 * c_frozen
            nn_loc = nn_loc + 1
          end if
        end do

        ! local accumulation
        do k = 2, ne
          wgt(k) = wgt(k-1) + wgt(k)
        end do

      else

        ! weights for mesh repartitioning ......................................

        if (logging) then
          print '(2A)', prefix, 'weights for mesh repartitioning'
        end if

        ! local accumulation using unit weights
        do k = 1, ne
          wgt(k) = k
        end do

        ! number of local contributions
        nn_loc = ne

      end if

      ! globalization ...........................................................

      if (logging) then
        print '(2A)', prefix, 'globalization'
      end if

      ! accumulated partition weights
      call MPI_Allgather( sendbuf   = wgt(ne)           &
                        , sendcount = 1                 &
                        , sendtype  = MPI_INTEGER       &
                        , recvbuf   = wgt_part          &
                        , recvcount = 1                 &
                        , recvtype  = MPI_INTEGER       &
                        , comm      = mesh % comm_parts )

      ! globalize accumulated weigths
      wgt_0 = sum(wgt_part(0:mesh%part-1))
      if (mesh % part > 0) then
        do k = 1, ne
          wgt(k) = wgt(k) + wgt_0
        end do
      end if
      wgt_sum = sum(wgt_part)

      ! total number of contributions
      call XMPI_Allreduce(nn_loc, nn, MPI_SUM, mesh%comm_parts)

      ! number of parts, avoiding empty partitions
      if (.not. opt%child) then
        n_parts = min( opt%n_parts, nn )
      else if (mesh%refinement == 'c') then
        n_parts = min( opt%n_parts, nn &
                     , nint(real(wgt_sum,RNP) / (c_active + c_frozen)) )
      else
        n_parts = min( opt%n_parts, nn &
                     , nint(real(wgt_sum,RNP) / (8 * (c_active + c_frozen))) )
      end if

      ! distribution of accumulated weights
      allocate(wgt_dist(0:n_parts-1))
      wgt_avg = wgt_sum / real(n_parts,RNP)
      do p = 0, n_parts-2
        wgt_dist(p) = nint(wgt_avg * (p+1))
      end do
      wgt_dist(n_parts-1) = wgt_sum

      ! target partitions ......................................................

      if (logging) then
        print '(2A)', prefix, 'target partitions'
      end if

      ! intialization
      tp_elem = -1

      ! start partition
      do p = n_parts-1, 1, -1
        if (wgt_0 > wgt_dist(p-1)) exit
      end do

      ! assign element target partitions
      do k = 1, ne

        if (opt % child) then
          ! skip elements with no children
          if (mesh % element(map(k)) % adaptation % mark < 100) cycle
        end if

        if (wgt(k) > wgt_dist(p)) then
          ! switch to next target partition
          p = p + 1
        end if

        tp_elem(map(k)) = p

      end do

      ! finalization ...........................................................

      deallocate(map, wgt, wgt_part, wgt_dist)

      if (logging) then
        print '(2A)', prefix, 'exit'
      end if

    end associate

  end subroutine SFC_Partitioner

  !=============================================================================

end submodule MP_SFC_Partitioner
