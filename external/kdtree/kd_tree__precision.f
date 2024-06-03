!
!(c) Matthew Kennel, Institute for Nonlinear Science (2004)
!
! Licensed under the Academic Free License version 1.1 found in file LICENSE
! with additional provisions found in that same file.
!
! Adapted to HiSPEET by Joerg Stiller (2024)

module KD_Tree__Precision
  implicit none
  public

! integer, parameter :: KDTree_RP = kind(1E0) ! single precision
  integer, parameter :: KDTree_RP = kind(1D0) ! double precision

end module KD_Tree__Precision
