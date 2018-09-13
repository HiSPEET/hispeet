!> summary:   Accelerated implementation of array assignments
!> author:    Joerg Stiller
!> date:      2016/12/09
!> license:   Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Accelerated implementation of array assignments
!>
!> @note
!> With OpenACC, all arguments must be present on the accelerator device.
!> @endnote
!===============================================================================

module Array_Assignments
  use Kind_Parameters, only: IXL, RNP
  use MPI_Binding
  implicit none
  private

  public :: AssignScalar
  public :: AssignArray
  public :: MergeArrays
  public :: ScaleArray
  public :: CalibrateArray

  interface AssignScalar
    module procedure AssignScalar_1
    module procedure AssignScalar_2
    module procedure AssignScalar_3
    module procedure AssignScalar_4
    module procedure AssignScalar_5
    module procedure AssignScalar_6
  end interface

  interface AssignArray
    module procedure AssignArray_3
    module procedure AssignArray_4
    module procedure AssignArray_5
  end interface

  interface MergeArrays
    module procedure MergeArrays_3
    module procedure MergeArrays_4
    module procedure MergeArrays_5
  end interface

  interface ScaleArray
    module procedure ScaleArray_3
    module procedure ScaleArray_4
    module procedure ScaleArray_5
  end interface

  interface CalibrateArray
    module procedure CalibrateArray_3
    module procedure CalibrateArray_4
    module procedure CalibrateArray_5
  end interface

contains

!===============================================================================
! AssignScalar

!-------------------------------------------------------------------------------
!> Assigment of a scalar s to 1D array a

subroutine AssignScalar_1(a, s)
  real(RNP),         intent(out) :: a(:)  !< target array
  real(RNP),         intent(in)  :: s     !< assigned scalar

  call AssignScalar_X(size(a), 1, a, s)

end subroutine AssignScalar_1

!-------------------------------------------------------------------------------
!> Assigment of a scalar s to 2D array a

subroutine AssignScalar_2(a, s, multi)
  real(RNP),         intent(out) :: a(:,:) !< target array
  real(RNP),         intent(in)  :: s      !< assigned scalar
  logical, optional, intent(in)  :: multi  !< switch to multiple components

  integer :: nc

  nc = 1
  if (present(multi)) then
    if (multi) nc = size(a,2)
  end if

  call AssignScalar_X(size(a)/nc, nc, a, s)

end subroutine AssignScalar_2

!-------------------------------------------------------------------------------
!> Assigment of a scalar s to 3D array a

subroutine AssignScalar_3(a, s, multi)
  real(RNP),         intent(out) :: a(:,:,:) !< target array
  real(RNP),         intent(in)  :: s        !< assigned scalar
  logical, optional, intent(in)  :: multi    !< switch to multiple components

  integer :: nc

  nc = 1
  if (present(multi)) then
    if (multi) nc = size(a,3)
  end if

  call AssignScalar_X(size(a)/nc, nc, a, s)

end subroutine AssignScalar_3

!-------------------------------------------------------------------------------
!> Assigment of a scalar s to 4D array a

subroutine AssignScalar_4(a, s, multi)
  real(RNP),         intent(out) :: a(:,:,:,:) !< target array
  real(RNP),         intent(in)  :: s          !< assigned scalar
  logical, optional, intent(in)  :: multi      !< switch to multiple components

  integer :: nc

  nc = 1
  if (present(multi)) then
    if (multi) nc = size(a,4)
  end if

  call AssignScalar_X(size(a)/nc, nc, a, s)

end subroutine AssignScalar_4

!-------------------------------------------------------------------------------
!> Assigment of a scalar s to 5D array a

subroutine AssignScalar_5(a, s, multi)
  real(RNP),         intent(out) :: a(:,:,:,:,:) !< target array
  real(RNP),         intent(in)  :: s            !< assigned scalar
  logical, optional, intent(in)  :: multi        !< switch to multiple comps.

  integer :: nc

  nc = 1
  if (present(multi)) then
    if (multi) nc = size(a,5)
  end if

  call AssignScalar_X(size(a)/nc, nc, a, s)

end subroutine AssignScalar_5

!-------------------------------------------------------------------------------
!> Assigment of a scalar s to 6D array a

subroutine AssignScalar_6(a, s, multi)
  real(RNP),         intent(out) :: a(:,:,:,:,:,:) !< target array
  real(RNP),         intent(in)  :: s              !< assigned scalar
  logical, optional, intent(in)  :: multi          !< switch to multiple comps.

  integer :: nc

  nc = 1
  if (present(multi)) then
    if (multi) nc = size(a,6)
  end if

  call AssignScalar_X(size(a)/nc, nc, a, s)

end subroutine AssignScalar_6

!-------------------------------------------------------------------------------
!> Assigment of a scalar s to multi-component array a (eXplicit)

subroutine AssignScalar_X(ne, nc, a, s)
  integer,   intent(in)  :: ne       !< number of entries per component
  integer,   intent(in)  :: nc       !< number of components
  real(RNP), intent(out) :: a(ne,nc) !< target array
  real(RNP), intent(in)  :: s        !< assigned scalar

  integer :: i, j

  !$acc parallel loop collapse(2) present(a)
  do j = 1, nc
    !$omp do
    do i = 1, ne
      a(i,j) = s
    end do
    !$omp end do nowait
  end do
  !$omp barrier

end subroutine AssignScalar_X

!===============================================================================
! AssignArray

!-------------------------------------------------------------------------------
!> Assigment of array b to array a (3D)

subroutine AssignArray_3(a, b, multi)
  real(RNP),         intent(out) :: a(:,:,:) !< target array
  real(RNP),         intent(in)  :: b(:,:,:) !< assigned array
  logical, optional, intent(in)  :: multi    !< switch to multiple components

  integer :: nc

  nc = 1
  if (present(multi)) then
    if (multi) nc = size(a,3)
  end if

  call AssignArray_X(size(a)/nc, nc, a, b)

end subroutine AssignArray_3

!-------------------------------------------------------------------------------
!> Assigment of array b to array a (4D)

subroutine AssignArray_4(a, b, multi)
  real(RNP),         intent(out) :: a(:,:,:,:) !< target array
  real(RNP),         intent(in)  :: b(:,:,:,:) !< assigned array
  logical, optional, intent(in)  :: multi      !< switch to multiple components

  integer :: nc

  nc = 1
  if (present(multi)) then
    if (multi) nc = size(a,4)
  end if

  call AssignArray_X(size(a)/nc, nc, a, b)

end subroutine AssignArray_4

!-------------------------------------------------------------------------------
!> Assigment of array b to array a (5D)

subroutine AssignArray_5(a, b, multi)
  real(RNP),         intent(out) :: a(:,:,:,:,:) !< target array
  real(RNP),         intent(in)  :: b(:,:,:,:,:) !< assigned array
  logical, optional, intent(in)  :: multi        !< switch to multiple comps.

  integer :: nc

  nc = 1
  if (present(multi)) then
    if (multi) nc = size(a,5)
  end if

  call AssignArray_X(size(a)/nc, nc, a, b)

end subroutine AssignArray_5

!-------------------------------------------------------------------------------
!> Assigment of array b to array a (multi-component eXplicit)

subroutine AssignArray_X(ne, nc, a, b)
  integer,   intent(in)  :: ne       !< number of entries per component
  integer,   intent(in)  :: nc       !< number of components
  real(RNP), intent(out) :: a(ne,nc) !< target array
  real(RNP), intent(in)  :: b(ne,nc) !< assigned array

  integer :: i, j

  !$acc parallel loop collapse(2) present(a,b)
  do j = 1, nc
    !$omp do
    do i = 1, ne
      a(i,j) = b(i,j)
    end do
    !$omp end do nowait
  end do
  !$omp barrier

end subroutine AssignArray_X

!===============================================================================
! MergeArrays

!-------------------------------------------------------------------------------
!> Performs  a = alpha * a + beta * b + s  with vectors a,b and scalar s (3D)

subroutine MergeArrays_3(alpha, a, beta, b, s, multi)
  real(RNP),           intent(in)    :: alpha    !< coefficient to a
  real(RNP),           intent(in)    :: beta     !< coefficient to b
  real(RNP),           intent(inout) :: a(:,:,:) !< target array
  real(RNP),           intent(in)    :: b(:,:,:) !< merged array
  real(RNP), optional, intent(in)    :: s        !< scalar [0]
  logical,   optional, intent(in)    :: multi    !< switch to multiple comps.

  integer :: nc

  nc = 1
  if (present(multi)) then
    if (multi) nc = size(a,3)
  end if

  call MergeArrays_X(size(a)/nc, nc, alpha, a, beta, b, s)

end subroutine MergeArrays_3

!-------------------------------------------------------------------------------
!> Performs  a = alpha * a + beta * b + s  with vectors a,b and scalar s (4D)

subroutine MergeArrays_4(alpha, a, beta, b, s, multi)
  real(RNP),           intent(in)    :: alpha      !< coefficient to a
  real(RNP),           intent(in)    :: beta       !< coefficient to b
  real(RNP),           intent(inout) :: a(:,:,:,:) !< target array
  real(RNP),           intent(in)    :: b(:,:,:,:) !< merged array
  real(RNP), optional, intent(in)    :: s          !< scalar [0]
  logical,   optional, intent(in)    :: multi      !< switch to multiple comps.

  integer :: nc

  nc = 1
  if (present(multi)) then
    if (multi) nc = size(a,4)
  end if

  call MergeArrays_X(size(a)/nc, nc, alpha, a, beta, b, s)

end subroutine MergeArrays_4

!-------------------------------------------------------------------------------
!> Performs  a = alpha * a + beta * b + s  with vectors a,b and scalar s (5D)

subroutine MergeArrays_5(alpha, a, beta, b, s, multi)
  real(RNP),           intent(in)    :: alpha        !< coefficient to a
  real(RNP),           intent(in)    :: beta         !< coefficient to b
  real(RNP),           intent(inout) :: a(:,:,:,:,:) !< target array
  real(RNP),           intent(in)    :: b(:,:,:,:,:) !< merged array
  real(RNP), optional, intent(in)    :: s            !< scalar [0]
  logical,   optional, intent(in)    :: multi        !< switch to multiple comps.

  integer :: nc

  nc = 1
  if (present(multi)) then
    if (multi) nc = size(a,5)
  end if

  call MergeArrays_X(size(a)/nc, nc, alpha, a, beta, b, s)

end subroutine MergeArrays_5

!-------------------------------------------------------------------------------
!> Performs  a = alpha * a + beta * b + s  with vectors a,b and scalar s
!> with multi-component arrays a and b

subroutine MergeArrays_X(ne, nc, alpha, a, beta, b, s)
  integer,             intent(in)    :: ne       !< num entries per component
  integer,             intent(in)    :: nc       !< num components
  real(RNP),           intent(in)    :: alpha    !< coefficient to a
  real(RNP),           intent(in)    :: beta     !< coefficient to b
  real(RNP),           intent(inout) :: a(ne,nc) !< target array
  real(RNP),           intent(in)    :: b(ne,nc) !< merged array
  real(RNP), optional, intent(in)    :: s        !< scalar [0]

  integer :: i, j

  real(RNP) :: c

  if (present(s)) then
    c = s
  else
    c = 0
  end if

  if (abs(alpha) < tiny(alpha)) then

    !$acc parallel loop collapse(2) present(a,b)
    do j = 1, nc
      !$omp do
      do i = 1, ne
        a(i,j) = beta * b(i,j) + c
      end do
      !$omp end do nowait
    end do
    !$omp barrier

  else

    !$acc parallel loop collapse(2) present(a,b)
    do j = 1, nc
      !$omp do
      do i = 1, ne
        a(i,j) = alpha * a(i,j) + beta * b(i,j) + c
      end do
      !$omp end do nowait
    end do
    !$omp barrier

  end if

end subroutine MergeArrays_X

!===============================================================================
! ScaleArray

!-------------------------------------------------------------------------------
!> Scale single-component array a by a factor of s (3D)

subroutine ScaleArray_3(a, s, multi)
  real(RNP),         intent(inout) :: a(:,:,:) !< target array
  real(RNP),         intent(in)    :: s        !< scaling factor
  logical, optional, intent(in)    :: multi    !< switch to multiple comps.

  integer :: nc

  nc = 1
  if (present(multi)) then
    if (multi) nc = size(a,3)
  end if

  call ScaleArray_X(size(a)/nc, nc, a, s)

end subroutine ScaleArray_3

!-------------------------------------------------------------------------------
!> Scale single-component array a by a factor of s (4D)

subroutine ScaleArray_4(a, s, multi)
  real(RNP),         intent(inout) :: a(:,:,:,:) !< target array
  real(RNP),         intent(in)    :: s          !< scaling factor
  logical, optional, intent(in)    :: multi      !< switch to multiple comps.

  integer :: nc

  nc = 1
  if (present(multi)) then
    if (multi) nc = size(a,4)
  end if

  call ScaleArray_X(size(a)/nc, nc, a, s)

end subroutine ScaleArray_4

!-------------------------------------------------------------------------------
!> Scale single-component array a by a factor of s (5D)

subroutine ScaleArray_5(a, s, multi)
  real(RNP),         intent(inout) :: a(:,:,:,:,:) !< target array
  real(RNP),         intent(in)    :: s            !< scaling factor
  logical, optional, intent(in)    :: multi        !< switch to multiple comps.

  integer :: nc

  nc = 1
  if (present(multi)) then
    if (multi) nc = size(a,5)
  end if

  call ScaleArray_X(size(a)/nc, nc, a, s)

end subroutine ScaleArray_5

!-------------------------------------------------------------------------------
!> Scale multi-component array a by a factor of s (eXplicit)

subroutine ScaleArray_X(ne, nc, a, s)
  integer,   intent(in)    :: ne       !< number of entries per component
  integer,   intent(in)    :: nc       !< number of entries per component
  real(RNP), intent(inout) :: a(ne,nc) !< target array
  real(RNP), intent(in)    :: s        !< scaling factor

  integer :: i, j

  !$acc parallel loop collapse(2) present(a)
  do j = 1, nc
    !$omp do
    do i = 1, ne
      a(i,j) = s * a(i,j)
    end do
    !$omp end do nowait
  end do
  !$omp barrier

end subroutine ScaleArray_X

!===============================================================================
! CalibrateArray

!-------------------------------------------------------------------------------
!> Remove non-zero arithmetic mean value (3D)

subroutine CalibrateArray_3(a, comm)
  real(RNP),                intent(inout) :: a(:,:,:) !< array
  type(MPI_Comm), optional, intent(in)    :: comm     !< MPI communicator

  call CalibrateArray_X(size(a), a, comm)

end subroutine CalibrateArray_3

!-------------------------------------------------------------------------------
!> Remove non-zero arithmetic mean value (4D)

subroutine CalibrateArray_4(a, comm)
  real(RNP),                intent(inout) :: a(:,:,:,:) !< array
  type(MPI_Comm), optional, intent(in)    :: comm       !< MPI communicator

  call CalibrateArray_X(size(a), a, comm)

end subroutine CalibrateArray_4

!-------------------------------------------------------------------------------
!> Remove non-zero arithmetic mean value (4D)

subroutine CalibrateArray_5(a, comm)
  real(RNP),                intent(inout) :: a(:,:,:,:,:) !< array
  type(MPI_Comm), optional, intent(in)    :: comm         !< MPI communicator

  call CalibrateArray_X(size(a), a, comm)

end subroutine CalibrateArray_5

!-------------------------------------------------------------------------------
!> Remove non-zero arithmetic mean value (eXplicit)

subroutine CalibrateArray_X(ne, a, comm)
  integer,                  intent(in)    :: ne    !< number of entries
  real(RNP),                intent(inout) :: a(ne) !< array
  type(MPI_Comm), optional, intent(in)    :: comm  !< MPI communicator

  ! local variables, some declared save to become shared with OpenMP
  real(RNP), save :: s_loc, s_glob, a_mean
  integer(IXL) :: n_loc, n_glob
  integer :: i

  ! mean value .................................................................

  s_loc = 0

  !$omp do reduction(+:s_loc)
  !$acc data present(a)
  !$acc parallel loop reduction(+:s_loc)
  do i = 1, ne
    s_loc = s_loc + a(i)
  end do
  !$acc end data

  ! OpenACC:
  ! Note that the reduction variable is automatically updated on the host!

  if (present(comm)) then

    !$omp master
    n_loc = ne
    call MPI_Allreduce(n_loc, n_glob, 1, MPI_INTEGER_IXL, MPI_SUM, comm)
    call MPI_Allreduce(s_loc, s_glob, 1, MPI_REAL_RNP   , MPI_SUM, comm)
    a_mean = s_glob / n_glob
    !$omp end master
    !$omp barrier

  else

    a_mean = s_loc / max(ne, 1)

  end if

  ! remove mean value ..........................................................

  !$omp do
  !$acc parallel loop present(a)
  do i = 1, ne
    a(i) = a(i) - a_mean
  end do

end subroutine CalibrateArray_X

!===============================================================================

end module Array_Assignments
