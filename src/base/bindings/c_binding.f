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

!> summary:  Collection of procedures promoting the access to C libraries
!> author:   Joerg Stiller
!> date:     2014/05/15
!>
!>### Collection of procedures promoting the access to C libraries.
!===============================================================================

module C_Binding
  use, intrinsic ::  ISO_C_Binding
  implicit none
  public

contains

!-------------------------------------------------------------------------------
!> Returns a C char array matching the given Fortran character variable

function CString(fs) result(cs)
  character(len=*), intent(in) :: fs  !< character variable
  character(C_CHAR) :: cs(len_trim(fs)+1)

  integer :: i, l

  l = size(cs)
  do i = 1, l-1
    cs(i) = fs(i:i)
  end do
  cs(l) = C_NULL_CHAR

end function CString

!-------------------------------------------------------------------------------
!> Returns a Fortran character variable matching the given C char array

function FString(cs) result(fs)
  character(C_CHAR), intent(in) :: cs(:) ! C string
  character(len=size(cs)-1) :: fs

  integer ::  i, l

  l = len(fs)
  fs = ''
  do i = 1, l
    if (cs(i) == C_NULL_CHAR) exit
    fs(i:i) = cs(i)
  end do

end function FString

!===============================================================================

end module C_Binding
