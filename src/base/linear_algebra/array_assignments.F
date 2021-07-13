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
  use Kind_Parameters, only: RNP
  use MPI_Binding
  use Execution_Control
  implicit none
  private

  public :: SetArray
  public :: MergeArrays
  public :: ScaleArray
  public :: CalibrateArray

  interface SetArray
    ! set array by assigning a scalar
    module procedure SetArray_1S
    module procedure SetArray_2S
    module procedure SetArray_3S
    module procedure SetArray_4S
    module procedure SetArray_5S
    module procedure SetArray_6S
    ! set array by assigning a matching array
    module procedure SetArray_3A
    module procedure SetArray_4A
    module procedure SetArray_5A
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
! SetArray

!-------------------------------------------------------------------------------
!> Assigment of a scalar s to 1D array a

subroutine SetArray_1S(a, s)
  real(RNP),         intent(out) :: a(:)  !< target array
  real(RNP),         intent(in)  :: s     !< assigned scalar

  call SetArray_XS(size(a), 1, a, s)

end subroutine SetArray_1S

!-------------------------------------------------------------------------------
!> Assigment of a scalar s to 2D array a

subroutine SetArray_2S(a, s, multi)
  real(RNP),         intent(out) :: a(:,:) !< target array
  real(RNP),         intent(in)  :: s      !< assigned scalar
  logical, optional, intent(in)  :: multi  !< switch to multiple components

  integer :: nc

  nc = 1
  if (present(multi)) then
    if (multi) nc = size(a,2)
  end if

  call SetArray_XS(size(a)/nc, nc, a, s)

end subroutine SetArray_2S

!-------------------------------------------------------------------------------
!> Assigment of a scalar s to 3D array a

subroutine SetArray_3S(a, s, multi)
  real(RNP),         intent(out) :: a(:,:,:) !< target array
  real(RNP),         intent(in)  :: s        !< assigned scalar
  logical, optional, intent(in)  :: multi    !< switch to multiple components

  integer :: nc

  nc = 1
  if (present(multi)) then
    if (multi) nc = size(a,3)
  end if

  call SetArray_XS(size(a)/nc, nc, a, s)

end subroutine SetArray_3S

!-------------------------------------------------------------------------------
!> Assigment of a scalar s to 4D array a

subroutine SetArray_4S(a, s, multi)
  real(RNP),         intent(out) :: a(:,:,:,:) !< target array
  real(RNP),         intent(in)  :: s          !< assigned scalar
  logical, optional, intent(in)  :: multi      !< switch to multiple components

  integer :: nc

  nc = 1
  if (present(multi)) then
    if (multi) nc = size(a,4)
  end if

  call SetArray_XS(size(a)/nc, nc, a, s)

end subroutine SetArray_4S

!-------------------------------------------------------------------------------
!> Assigment of a scalar s to 5D array a

subroutine SetArray_5S(a, s, multi)
  real(RNP),         intent(out) :: a(:,:,:,:,:) !< target array
  real(RNP),         intent(in)  :: s            !< assigned scalar
  logical, optional, intent(in)  :: multi        !< switch to multiple comps.

  integer :: nc

  nc = 1
  if (present(multi)) then
    if (multi) nc = size(a,5)
  end if

  call SetArray_XS(size(a)/nc, nc, a, s)

end subroutine SetArray_5S

!-------------------------------------------------------------------------------
!> Assigment of a scalar s to 6D array a

subroutine SetArray_6S(a, s, multi)
  real(RNP),         intent(out) :: a(:,:,:,:,:,:) !< target array
  real(RNP),         intent(in)  :: s              !< assigned scalar
  logical, optional, intent(in)  :: multi          !< switch to multiple comps.

  integer :: nc

  nc = 1
  if (present(multi)) then
    if (multi) nc = size(a,6)
  end if

  call SetArray_XS(size(a)/nc, nc, a, s)

end subroutine SetArray_6S

!-------------------------------------------------------------------------------
!> Assigment of a scalar s to multi-component array a (eXplicit)

subroutine SetArray_XS(ne, nc, a, s)
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

end subroutine SetArray_XS

!-------------------------------------------------------------------------------
!> Assigment of array b to array a (3D)

subroutine SetArray_3A(a, b, multi)
  real(RNP),         intent(out) :: a(:,:,:) !< target array
  real(RNP),         intent(in)  :: b(:,:,:) !< assigned array
  logical, optional, intent(in)  :: multi    !< switch to multiple components

  integer :: nc

#ifdef DEBUG
  if (any(shape(a) /= shape(b))) then
    call Error('SetArray_3A', 'shape mismatch', 'Array_Assignments')
  end if
#endif

  nc = 1
  if (present(multi)) then
    if (multi) nc = size(a,3)
  end if

  call SetArray_XA(size(a)/nc, nc, a, b)

end subroutine SetArray_3A

!-------------------------------------------------------------------------------
!> Assigment of array b to array a (4D)

subroutine SetArray_4A(a, b, multi)
  real(RNP),         intent(out) :: a(:,:,:,:) !< target array
  real(RNP),         intent(in)  :: b(:,:,:,:) !< assigned array
  logical, optional, intent(in)  :: multi      !< switch to multiple components

  integer :: nc

#ifdef DEBUG
  if (any(shape(a) /= shape(b))) then
    call Error('SetArray_4A', 'shape mismatch', 'Array_Assignments')
  end if
#endif

  nc = 1
  if (present(multi)) then
    if (multi) nc = size(a,4)
  end if

  call SetArray_XA(size(a)/nc, nc, a, b)

end subroutine SetArray_4A

!-------------------------------------------------------------------------------
!> Assigment of array b to array a (5D)

subroutine SetArray_5A(a, b, multi)
  real(RNP),         intent(out) :: a(:,:,:,:,:) !< target array
  real(RNP),         intent(in)  :: b(:,:,:,:,:) !< assigned array
  logical, optional, intent(in)  :: multi        !< switch to multiple comps.

  integer :: nc

#ifdef DEBUG
  if (any(shape(a) /= shape(b))) then
    call Error('SetArray_5A', 'shape mismatch', 'Array_Assignments')
  end if
#endif

  nc = 1
  if (present(multi)) then
    if (multi) nc = size(a,5)
  end if

  call SetArray_XA(size(a)/nc, nc, a, b)

end subroutine SetArray_5A

!-------------------------------------------------------------------------------
!> Assigment of array b to array a (multi-component eXplicit)

subroutine SetArray_XA(ne, nc, a, b)
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

end subroutine SetArray_XA

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

#ifdef DEBUG
  if (any(shape(a) /= shape(b))) then
    call Error('MergeArrays_3', 'shape mismatch', 'Array_Assignments')
  end if
#endif

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

#ifdef DEBUG
  if (any(shape(a) /= shape(b))) then
    call Error('MergeArrays_4', 'shape mismatch', 'Array_Assignments')
  end if
#endif

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

#ifdef DEBUG
  if (any(shape(a) /= shape(b))) then
    call Error('MergeArrays_5', 'shape mismatch', 'Array_Assignments')
  end if
#endif

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
  real(RNP), save :: s_loc = 0, s_glob = 0, a_mean
  integer :: n_loc, n_glob
  integer :: i

  ! mean value .................................................................

  !$omp do reduction(+:s_loc)
  !$acc data present(a)
  !$acc parallel loop reduction(+:s_loc)
  do i = 1, ne
    s_loc = s_loc + a(i)
  end do
  !$acc end data

  ! OpenACC:
  ! Note that the reduction variable is automatically updated on the host!

  !$omp master
  if (present(comm)) then
    n_loc = ne
    call MPI_Allreduce(n_loc, n_glob, 1, MPI_INTEGER , MPI_SUM, comm)
    call MPI_Allreduce(s_loc, s_glob, 1, MPI_REAL_RNP, MPI_SUM, comm)
    a_mean = s_glob / n_glob
  else
    a_mean = s_loc / max(ne, 1)
  end if
  s_loc = 0
  !$omp end master
  !$omp barrier

  ! remove mean value ..........................................................

  !$acc parallel loop present(a)
  !$omp do
  do i = 1, ne
    a(i) = a(i) - a_mean
  end do

end subroutine CalibrateArray_X

!===============================================================================

end module Array_Assignments
