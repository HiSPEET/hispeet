!> summary:  Collection of procedures promoting the access to C libraries
!> author:   Joerg Stiller
!> date:     2014/05/15
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
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