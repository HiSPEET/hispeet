!> summary:  Control of logging levels
!> author:   Joerg Stiller
!> date:     2023/02/18
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
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

  !=============================================================================

end module Logging_Levels
