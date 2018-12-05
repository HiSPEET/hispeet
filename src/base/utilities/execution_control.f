!> summary:  Exception handling and execution control
!> author:   Joerg Stiller
!> date:     2013/05/12
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Exception handling and execution control accounting for MPI and OpenMP
!>
!> Note that with OpenMP, all instructions are executed by the master thread.
!===============================================================================

module Execution_Control
  use, intrinsic :: ISO_Fortran_Env, only: ERR => ERROR_UNIT, OU => OUTPUT_UNIT
  use MPI_F08
  implicit none
  private

  public :: Error
  public :: Warning
  public :: TracingMessage


  !-----------------------------------------------------------------------------
  ! Private entities

  ! control parameters
  logical :: initialized = .false.
  logical :: parallel    = .false.
  integer :: pid         = -1

  ! switch-on trace level
  integer :: min_trace_level = 0

contains

!-------------------------------------------------------------------------------
!> Initialization of the module

subroutine Init_Execution_Control

  if (initialized) return
  call MPI_Initialized(parallel)
  if (parallel) then
    call MPI_Comm_rank(MPI_COMM_WORLD, pid)
  end if

  initialized = .true.

end subroutine Init_Execution_Control

!-------------------------------------------------------------------------------
!> Error handling in a parallel excecution environment.

subroutine Error(procedure, message, module)
  use, intrinsic :: ISO_Fortran_Env, only: ERR => ERROR_UNIT
  character(LEN=*),           intent(in) :: procedure !< procedure name
  character(LEN=*),           intent(in) :: message   !< error message
  character(LEN=*), optional, intent(in) :: module    !< module name

  !$omp master

  if (.not. initialized) then
    call Init_Execution_Control
  end if

  if (present(module)) then
    write(ERR,'(/,5A)', advance='NO') 'ERROR (', module,'::', procedure, ') '
  else
    write(ERR,'(/,3A)', advance='NO') 'ERROR (', procedure, ') '
  end if
  if (parallel) then
    write(ERR,'(A,I0)', advance='NO') 'in process ', pid
  end if
  write(ERR,'(2A,/)') ': ', message

  flush(ERR)

  if (parallel) then
    call MPI_Abort(MPI_COMM_WORLD, 1)
  end if
  stop

  !$omp end master

end subroutine Error

!-------------------------------------------------------------------------------
!> Generates a warning in a parallel excecution environment.

subroutine Warning(procedure, message, module)
  character(LEN=*),           intent(in) :: procedure !< procedure name
  character(LEN=*),           intent(in) :: message   !< warning message
  character(LEN=*), optional, intent(in) :: module    !< module name

  !$omp master

  if (.not. initialized) then
    call Init_Execution_Control
  end if

  if (present(module)) then
    write(ERR,'(/,5A)', advance='NO') 'WARNING (', module,'::', procedure, ') '
  else
    write(ERR,'(/,3A)', advance='NO') 'WARNING (', procedure, ') '
  end if
  if (parallel) then
    write(ERR,'(A,I0)', advance='NO') 'in process ', pid
  end if
  write(ERR,'(2A,/)') ': ', message

  flush(ERR)

  !$omp end master

end subroutine Warning

!-------------------------------------------------------------------------------
!> Print a tracing message
!>
!> The message will be printed if its level is greater than or equal to
!> min_trace_level.

subroutine TracingMessage(level, label, pos, message)
  integer,                    intent(in) :: level   !< message level
  character(len=*),           intent(in) :: label   !< identification string
  integer,          optional, intent(in) :: pos     !< identification number
  character(len=*), optional, intent(in) :: message !< additional information

  character(len=80) :: label_, pid_, pos_, message_

  if (level < min_trace_level) return

  !$omp master

  if (.not. initialized) then
    call Init_Execution_Control
  end if

  label_   = ''
  pid_     = ''
  pos_     = ''
  message_ = ''

  write(label_,'(A)') label
  if (parallel) then
    write(pid_,'(A,I0)') ' in process ', pid
  end if
  if (present(pos)) then
    write(pos_,'(A,I0)') ', pos ', pos
  end if
  if (present(message)) then
    write(message_,'(2A)') ': ', trim(message)
  end if
  write(OU,'(4A)') trim(label_), trim(pid_), trim(pos_), trim(message_)

  flush(OU)

  !$omp end master

end subroutine TracingMessage

!===============================================================================

end module Execution_Control