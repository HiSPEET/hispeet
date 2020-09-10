!> summary:  Wrappers for selected procedures defined in the OpenMP module
!> author:   Immo Huismann, Joerg Stiller
!> date:     2014/08/28
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> This module provides the OpenMP functionalities to determine the number of
!> threads and the thread ID.
!==============================================================================

module OpenMP_Binding
#ifdef _OPENMP
  use OMP_Lib, only: OMP_Get_Num_Threads,  &
                     OMP_Get_Thread_Num
#endif

  implicit none
  private

  public :: OMP_Num_Threads
  public :: OMP_Thread_ID

contains

  !-----------------------------------------------------------------------------
  !> Returns the current number of threads.

  integer function OMP_Num_Threads() result(num_threads)

#ifdef _OPENMP
    num_threads = OMP_Get_Num_Threads()
#else
    num_threads = 1
#endif

  end function OMP_Num_Threads

  !-----------------------------------------------------------------------------
  !> Returns the current thread ID

  integer function OMP_Thread_ID() result(thread_id)

#ifdef _OPENMP
    thread_id = OMP_Get_Thread_Num()
#else
    thread_id = 0
#endif

  end function OMP_Thread_ID

  !=============================================================================

end module OpenMP_Binding
