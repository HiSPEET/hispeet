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

!> summary:  Wrappers for selected procedures defined in the OpenMP module
!> author:   Immo Huismann, Joerg Stiller
!> date:     2014/08/28
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
