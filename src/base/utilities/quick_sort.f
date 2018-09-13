!> summary:  Quick sort algorithm for vectors, pairs and triplets
!> author:   Joerg Stiller
!> date:     2014/01/10
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Quick sort algorithm for vectors, pairs and triplets
!>
!> @todo
!>
!>   *  real version of Sort
!>   *  generic implementation once Fortran is mature enough
!===============================================================================

module Quick_Sort
  use Kind_Parameters,   only: RNP
  implicit none
  private

  public :: Sort
  public :: SortIndex
  public :: SortPairs
  public :: SortTriplets

  !-----------------------------------------------------------------------------
  !> Sorts a given one-dimensional array
  !>
  !> Returns the given array x such that
  !>
  !>     x(j) <= x(k)  for all j < k

  interface Sort
    module procedure  IntegerSort
  end interface Sort

  !-----------------------------------------------------------------------------
  !> Determines a permutation ordering a given one-dimensional array.
  !>
  !> Given a 1D array x the procedure returns a permutation vector i such that
  !>
  !>     x(i(j)) <= x(i(k))  for all j < k

  interface SortIndex
    module procedure  IntegerSortIndex
    module procedure  RealSortIndex
  end interface

  !-----------------------------------------------------------------------------
  !> Sorts a given two-dimensional array according to columns i1:i1+1.
  !>
  !> Sorts the given 2D array x such that
  !>
  !>         x(i1, i) <  x(i1, k)
  !>     or
  !>         x(i1  , i) == x(i1  , k)  and
  !>         x(i1+1, i) <= x(i1+1, k)
  !>
  !> The optional argument i1 defaults to 1.

  interface SortPairs
    module procedure  IntegerSortPairs
  end interface SortPairs

  !-----------------------------------------------------------------------------
  !> Sorts a given two-dimensional according to first three columns
  !>
  !> Sorts the given 2D array x such that for all  i < k
  !>
  !>         x(1, i) <  x(1, k)
  !>     or
  !>         x(1, i) == x(1, k)  and
  !>         x(2, i) <  x(2, k)
  !>     or
  !>         x(1, i) == x(1, k)  and
  !>         x(2, i) == x(2, k)  and
  !>         x(3, i) <= x(3, k)

  interface SortTriplets
    module procedure  IntegerSortTriplets
  end interface SortTriplets

contains

!-------------------------------------------------------------------------------
!> Returns x such that x(j) <= x(k) for all j < k

pure recursive subroutine IntegerSort(x, i1, i2)
  integer,           intent(inout) :: x(:) !< array, sorted on output
  integer, optional, intent(in)    :: i1   !< start index
  integer, optional, intent(in)    :: i2   !< terminal index

  integer :: j1, j2, j10, j20

  if (size(x) < 2) return

  select case (present(i1) .and. present(i2))
  case (.false.) ! assuming first call
    j10 = lbound(x,1)
    j20 = ubound(x,1)
  case (.true.)  ! assuming secondary call
    if (i1 >= i2) return
    j10 = i1
    j20 = i2
  end select
  j1 = j10
  j2 = j20

  search: do
    do
      if ( x(j1) > x(j2) ) then
        x([j1,j2]) = x([j2,j1])
        exit
      end if
      j2 = j2-1
      if (j1 == j2) exit search
    end do
    do
      if ( x(j1) > x(j2) ) then
        x([j1,j2]) = x([j2,j1])
        exit
      end if
      j1 = j1+1
      if (j1 == j2) exit search
    end do
  end do search

  call IntegerSort(x, j10 , j1-1)
  call IntegerSort(x, j2+1, j20 )

end subroutine IntegerSort

!===============================================================================

!-------------------------------------------------------------------------------
!> Returns i such that x(i(j)) <= x(i(k)) for all j < k

pure recursive subroutine IntegerSortIndex(x, i, i1, i2)
  integer,           intent(in)    :: x(:)       !< array to be sorted
  integer,           intent(inout) :: i(size(x)) !< ordered indices of x
  integer, optional, intent(in)    :: i1         !< start index
  integer, optional, intent(in)    :: i2         !< terminal index

  integer :: j, j1, j2, j10, j20

  if (size(x) < 2) then
    i = 1
    return
  end if

  select case (present(i1) .and. present(i2))
  case (.false.) ! assuming first call
    j10 = lbound(x,1)
    j20 = ubound(x,1)
    i = [ (j, j=j10,j20) ]
  case (.true.)  ! assuming secondary call
    if (i1 >= i2) return
    j10 = i1
    j20 = i2
  end select
  j1 = j10
  j2 = j20

  search: do
    do
      if ( x(i(j1)) > x(i(j2)) ) then
        i([j1,j2]) = i([j2,j1])
        exit
      end if
      j2 = j2-1
      if (j1 == j2) exit search
    end do
    do
      if ( x(i(j1)) > x(i(j2)) ) then
        i([j1,j2]) = i([j2,j1])
        exit
      end if
      j1 = j1+1
      if (j1 == j2) exit search
    end do
  end do search

  call IntegerSortIndex(x, i, j10 , j1-1)
  call IntegerSortIndex(x, i, j2+1, j20 )
end subroutine IntegerSortIndex

!-------------------------------------------------------------------------------
!> Returns i such that x(i(j)) <= x(i(k)) for all j < k

pure recursive subroutine RealSortIndex(x, i, i1, i2)
  real(RNP),         intent(in)    :: x(:)       !< array to be sorted
  integer,           intent(inout) :: i(size(x)) !< ordered indices of x
  integer, optional, intent(in)    :: i1         !< start index
  integer, optional, intent(in)    :: i2         !< terminal index

  integer :: j, j1, j2, j10, j20

  if (size(x) < 2) then
    i = 1
    return
  end if

  select case (present(i1) .and. present(i2))
  case (.false.) ! assuming first call
    j10 = lbound(x,1)
    j20 = ubound(x,1)
    i = [ (j, j=j10,j20) ]
  case (.true.)  ! assuming secondary call
    if (i1 >= i2) return
    j10 = i1
    j20 = i2
  end select
  j1 = j10
  j2 = j20

  search: do
    do
      if ( x(i(j1)) > x(i(j2)) ) then
        i([j1,j2]) = i([j2,j1])
        exit
      end if
      j2 = j2-1
      if (j1 == j2) exit search
    end do
    do
      if ( x(i(j1)) > x(i(j2)) ) then
        i([j1,j2]) = i([j2,j1])
        exit
      end if
      j1 = j1+1
      if (j1 == j2) exit search
    end do
  end do search

  call RealSortIndex(x, i, j10 , j1-1)
  call RealSortIndex(x, i, j2+1, j20 )
end subroutine RealSortIndex

!===============================================================================

!-------------------------------------------------------------------------------
!> Sorts x(:,:) according to columns i1:i1+1

pure subroutine IntegerSortPairs(x, i1)
  integer, intent(inout)        :: x(:,:) !< array being sorted, size(x,1) > i1
  integer, intent(in), optional :: i1     !< first column to sort, 1 if absent

  integer, dimension(size(x,2)) :: p
  integer :: i, j, k, s

  if (size(p) < 2) return

  ! define index shift
  if (present(i1)) then
    s = max(0, i1-1)
  else
    s = 0
  end if

  ! sort x with respect to first index
  call SortIndex(x(1+s,:), p)
  x(:,:) = x(:,p)

  ! sort sections of x with same first index
  k = 1+s
  i = 1
  do
    do  j = i+1, size(p)
      if (x(k,j) /= x(k,i)) exit
    end do
    j = j - 1
    if (j - i > 0) then
      call SortIndex(x(k+1,i:j), p(i:j))
      p(i:j)   = p(i:j) + i - 1
      x(:,i:j) = x(:,p(i:j))
    end if
    i = j + 1
     if (i >= size(p)) exit
  end do

end subroutine IntegerSortPairs

!===============================================================================

!-------------------------------------------------------------------------------
!> Sorts x(:,:) according to first three columns

pure subroutine IntegerSortTriplets(x)
  integer, intent(inout) :: x(:,:) !< array of triplets being sorted

  integer, dimension(size(x,2)) :: p  ! permutation vector
  integer :: i, j

  if (size(x,1) < 3) return
  if (size(x,2) < 2) return

  ! sort x with respect to first index
  call SortIndex(x(1,:), p)
  x(:,:) = x(:,p)

  ! sort sections of x with same first index
  i = 1
  do
    do  j = i+1, size(p)
      if (x(1,j) /= x(1,i)) exit
    end do
    j = j - 1
    if (j - i > 0) then
      call IntegerSortPairs(x(:,i:j), i1=2)
    end if
    i = j + 1
    if (i >= size(p)) exit
  end do

end subroutine IntegerSortTriplets

!=============================================================================

end module Quick_Sort