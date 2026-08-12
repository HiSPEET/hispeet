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

!> summary:  Control of logging levels
!> author:   Joerg Stiller
!> date:     2023/02/18
!===============================================================================

module Logging_Levels
  use XMPI
  implicit none

  integer :: log_level = 0
  integer :: log_level_inner_iteration = 0
  integer :: log_level_outer_iteration = 0
  integer :: log_level_multigrid_cycle = 0

contains

  !-----------------------------------------------------------------------------
  !> MPI broadcasting of logging levels

  subroutine XMPI_Bcast_LoggingLevels(root, comm)
    integer,        intent(in) :: root !< rank of broadcast root
    type(MPI_Comm), intent(in) :: comm !< MPI communicator

    type(MPI_Request) :: request(4)

    call XMPI_Ibcast( log_level                , root, comm, request(1) )
    call XMPI_Ibcast( log_level_inner_iteration, root, comm, request(2) )
    call XMPI_Ibcast( log_level_outer_iteration, root, comm, request(3) )
    call XMPI_Ibcast( log_level_multigrid_cycle, root, comm, request(4) )

    call MPI_Waitall( size(request), request, MPI_STATUSES_IGNORE )

  end subroutine XMPI_Bcast_LoggingLevels

  !-----------------------------------------------------------------------------
  !> Unique prefix for logging output

  function LoggingPrefix(name, rank, n_proc) result(prefix)
    character(len=*),  intent(in) :: name   !< program of procedure name
    integer, optional, intent(in) :: rank   !< process rank
    integer, optional, intent(in) :: n_proc !< number of processes
    character(len=:), allocatable :: prefix

    integer :: l

    if (present(rank)) then

      if (present(n_proc)) then
        l = int(log10(dble(n_proc))) + 1
      else
        l = max(int(log10(dble(max(rank,1)))), 2) + 1
      end if
      allocate(character(len=l+2) :: prefix)
      write(prefix,'(I0,A)') rank, ':'
      prefix = '> ' // trim(name) // ' #' // prefix

    else

      prefix = '> ' // trim(name) // ': '

    end if

  end function LoggingPrefix

  !=============================================================================

end module Logging_Levels
