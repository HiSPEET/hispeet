!> summary:  A small module to enable simple cache blocking
!> author:   Immo Huismann, Joerg Stiller
!> date:     2016/12/03
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### A small module to enable simple cache blocking
!===============================================================================

module Cache_Blocking
  implicit none
  private

  public :: DeviseCacheBlocks

  integer :: double_size   = storage_size(1D0)/8  !< byte size of doubles
  integer :: cache_size(2) = [256, 2048]*1024     !< cache sizes in bytes
  integer :: l_block_min   = 5                    !< minimum block length
  real    :: safety_factor = 2                    !< safety factor

contains

!-------------------------------------------------------------------------------
!> Computes the number of tasks per cache block and the number of blocks

subroutine DeviseCacheBlocks(s_task, n_task, l_block, n_block)
  integer, intent(in)  :: s_task   !< size of one task in terms of used doubles
  integer, intent(in)  :: n_task   !< total number of tasks
  integer, intent(out) :: l_block  !< length of one block
  integer, intent(out) :: n_block  !< number of blocks

  integer :: tasks_per_cache(2)

  tasks_per_cache = int( cache_size/(s_task * double_size) * safety_factor )

  if (tasks_per_cache(1) >= l_block_min) then
    l_block = tasks_per_cache(1)
  else if (tasks_per_cache(2) >= l_block_min) then
    l_block = tasks_per_cache(2)
  else
    l_block = n_task
  end if

  n_block = 1 + (n_task - 1) / l_block

end subroutine DeviseCacheBlocks

!===============================================================================

end module Cache_Blocking