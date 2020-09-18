!> summary:  Parallel accelerated array reductions
!> author:   Joerg Stiller
!> date:     2017/06/28
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Parallel accelerated array reductions
!===============================================================================

module Array_Reductions
  use Kind_Parameters, only: RNP
  use MPI_Binding
  implicit none
  private

  public :: SumOfElements
  public :: ScalarProduct
  public :: WeightedScalarProduct

  interface SumOfElements
    module procedure SumOfElements_3
    module procedure SumOfElements_4
    module procedure SumOfElements_5
  end interface

  interface ScalarProduct
    module procedure ScalarProduct_3
    module procedure ScalarProduct_4
    module procedure ScalarProduct_5
  end interface

  interface WeightedScalarProduct
    module procedure WeightedScalarProduct_3
    module procedure WeightedScalarProduct_4
    module procedure WeightedScalarProduct_5
  end interface

contains

!===============================================================================
! SumOfElements

!-------------------------------------------------------------------------------
!> Sum over all array elements (3D)
!>
!> If the MPI communicator comm is passed, the sum will be evaluated over all
!> processes. Otherwise, only the local contribution is returned.

real(RNP) function SumOfElements_3(a, comm) result(r)
  real(RNP),                intent(in) :: a(:,:,:) !< array
  type(MPI_Comm), optional, intent(in) :: comm     !< MPI communicator

  r = SumOfElements_X(size(a), a, comm)

end function SumOfElements_3

!-------------------------------------------------------------------------------
!> Sum over all array elements (4D)
!>
!> If the MPI communicator comm is passed, the sum will be evaluated over all
!> processes. Otherwise, only the local contribution is returned.

real(RNP) function SumOfElements_4(a, comm) result(r)
  real(RNP),                intent(in) :: a(:,:,:,:) !< array
  type(MPI_Comm), optional, intent(in) :: comm       !< MPI communicator

  r = SumOfElements_X(size(a), a, comm)

end function SumOfElements_4

!-------------------------------------------------------------------------------
!> Sum over all array elements (5D)
!>
!> If the MPI communicator comm is passed, the sum will be evaluated over all
!> processes. Otherwise, only the local contribution is returned.

real(RNP) function SumOfElements_5(a, comm) result(r)
  real(RNP),                intent(in) :: a(:,:,:,:,:) !< array
  type(MPI_Comm), optional, intent(in) :: comm         !< MPI communicator

  r = SumOfElements_X(size(a), a, comm)

end function SumOfElements_5

!-------------------------------------------------------------------------------
!> Sum over all array elements (eXplicit)
!>
!> If the MPI communicator comm is passed, the sum will be evaluated over all
!> processes. Otherwise, only the local contribution is returned.

real(RNP) function SumOfElements_X(n, a, comm) result(r)
  integer,                  intent(in) :: n    !< array size
  real(RNP),                intent(in) :: a(n) !< array
  type(MPI_Comm), optional, intent(in) :: comm !< MPI communicator

  ! local variables declared save to become shared with OpenMP
  real(RNP), save :: r_loc, r_glob
  integer  :: i

  ! local contribution .........................................................

  r_loc = 0

  !$omp do reduction(+:r_loc)
  !$acc data present(a)
  !$acc parallel loop reduction(+:r_loc)
  do i = 1, n
    r_loc = r_loc + a(i)
  end do
  !$acc end data

  ! OpenACC:
  ! Note that the reduction variable is automatically updated on the host!

  ! global result ..............................................................

  if (present(comm)) then

    !$omp master
    call MPI_Allreduce(r_loc, r_glob, 1, MPI_REAL_RNP, MPI_SUM, comm)
    !$omp end master
    !$omp barrier

    r = r_glob

  else

    r = r_loc

  end if

end function SumOfElements_X

!===============================================================================
! ScalarProduct

!-------------------------------------------------------------------------------
!> Scalar product of vectors a and b (3D)
!>
!> If the MPI communicator comm is passed, the scalar product will be evaluated
!> over all processes. Otherwise, only the local contribution is returned.

real(RNP) function ScalarProduct_3(a, b, comm) result(r)
  real(RNP),                intent(in) :: a(:,:,:) !< left operand
  real(RNP),                intent(in) :: b(:,:,:) !< right operand
  type(MPI_Comm), optional, intent(in) :: comm     !< MPI communicator

  r = ScalarProduct_X(size(a), a, b, comm)

end function ScalarProduct_3

!-------------------------------------------------------------------------------
!> Scalar product of vectors a and b (4D)
!>
!> If the MPI communicator comm is passed, the scalar product will be evaluated
!> over all processes. Otherwise, only the local contribution is returned.

real(RNP) function ScalarProduct_4(a, b, comm) result(r)
  real(RNP),                intent(in) :: a(:,:,:,:) !< left operand
  real(RNP),                intent(in) :: b(:,:,:,:) !< right operand
  type(MPI_Comm), optional, intent(in) :: comm       !< MPI communicator

  r = ScalarProduct_X(size(a), a, b, comm)

end function ScalarProduct_4

!-------------------------------------------------------------------------------
!> Scalar product of vectors a and b (5D)
!>
!> If the MPI communicator comm is passed, the scalar product will be evaluated
!> over all processes. Otherwise, only the local contribution is returned.

real(RNP) function ScalarProduct_5(a, b, comm) result(r)
  real(RNP),                intent(in) :: a(:,:,:,:,:) !< left operand
  real(RNP),                intent(in) :: b(:,:,:,:,:) !< right operand
  type(MPI_Comm), optional, intent(in) :: comm         !< MPI communicator

  r = ScalarProduct_X(size(a), a, b, comm)

end function ScalarProduct_5

!-------------------------------------------------------------------------------
!> Scalar product of vectors a and b (eXplicit)
!>
!> If the MPI communicator comm is passed, the scalar product will be evaluated
!> over all processes. Otherwise, only the local contribution is returned.

real(RNP) function ScalarProduct_X(n, a, b, comm) result(r)
  integer,                  intent(in) :: n    !< array size
  real(RNP),                intent(in) :: a(n) !< left operand
  real(RNP),                intent(in) :: b(n) !< right operand
  type(MPI_Comm), optional, intent(in) :: comm !< MPI communicator

  ! local variables declared save to become shared with OpenMP
  real(RNP), save :: r_loc, r_glob
  real(RNP) :: r_priv
  integer :: i

  ! local contribution .........................................................

  !$omp single
  r_loc = 0
  !$omp end single

  r_priv = 0

  !$omp do
  do i = 1, n
    r_priv = r_priv + a(i) * b(i)
  end do
  !$omp end do nowait

  !$omp critical
  r_loc = r_loc + r_priv
  !$omp end critical
  !$omp barrier

  ! global result ..............................................................

  if (present(comm)) then

    !$omp master
    call MPI_Allreduce(r_loc, r_glob, 1, MPI_REAL_RNP, MPI_SUM, comm)
    !$omp end master
    !$omp barrier

    r = r_glob

  else

    r = r_loc

  end if

end function ScalarProduct_X

!===============================================================================
! WeightedScalarProduct

!-------------------------------------------------------------------------------
!> Weighted scalar productof vectors a and b (3D)
!>
!> If the MPI communicator comm is passed, the scalar product will be evaluated
!> over all processes. Otherwise, only the local contribution is returned.

real(RNP) function WeightedScalarProduct_3(w, a, b, comm) result(r)
  real(RNP),                intent(in) :: w(:,:,:) !< weights
  real(RNP),                intent(in) :: a(:,:,:) !< left operand
  real(RNP),                intent(in) :: b(:,:,:) !< right operand
  type(MPI_Comm), optional, intent(in) :: comm     !< MPI communicator

  r = WeightedScalarProduct_X(size(w), w, a, b, comm)

end function WeightedScalarProduct_3

!-------------------------------------------------------------------------------
!> Weighted scalar productof vectors a and b (4D)
!>
!> If the MPI communicator comm is passed, the scalar product will be evaluated
!> over all processes. Otherwise, only the local contribution is returned.

real(RNP) function WeightedScalarProduct_4(w, a, b, comm) result(r)
  real(RNP),                intent(in) :: w(:,:,:,:) !< weights
  real(RNP),                intent(in) :: a(:,:,:,:) !< left operand
  real(RNP),                intent(in) :: b(:,:,:,:) !< right operand
  type(MPI_Comm), optional, intent(in) :: comm       !< MPI communicator

  r = WeightedScalarProduct_X(size(w), w, a, b, comm)

end function WeightedScalarProduct_4

!-------------------------------------------------------------------------------
!> Weighted scalar productof vectors a and b (5D)
!>
!> If the MPI communicator comm is passed, the scalar product will be evaluated
!> over all processes. Otherwise, only the local contribution is returned.

real(RNP) function WeightedScalarProduct_5(w, a, b, comm) result(r)
  real(RNP),                intent(in) :: w(:,:,:,:,:) !< weights
  real(RNP),                intent(in) :: a(:,:,:,:,:) !< left operand
  real(RNP),                intent(in) :: b(:,:,:,:,:) !< right operand
  type(MPI_Comm), optional, intent(in) :: comm         !< MPI communicator

  r = WeightedScalarProduct_X(size(w), w, a, b, comm)

end function WeightedScalarProduct_5

!-------------------------------------------------------------------------------
!> Weighted scalar productof vectors a and b (eXplicit)
!>
!> If the MPI communicator comm is passed, the scalar product will be evaluated
!> over all processes. Otherwise, only the local contribution is returned.

real(RNP) function WeightedScalarProduct_X(n, w, a, b, comm) result(r)
  integer,                  intent(in) :: n    !< array size
  real(RNP),                intent(in) :: w(n) !< weights
  real(RNP),                intent(in) :: a(n) !< left operand
  real(RNP),                intent(in) :: b(n) !< right operand
  type(MPI_Comm), optional, intent(in) :: comm !< MPI communicator

  ! local variables declared save to become shared with OpenMP
  real(RNP), save :: r_loc, r_glob
  real(RNP) :: r_priv
  integer  :: i

  ! local contribution .........................................................

  !$omp single
  r_loc = 0
  !$omp end single

  r_priv = 0

  !$omp do
  do i = 1, n
    r_priv = r_priv + w(i) * a(i) * b(i)
  end do
  !$omp end do nowait

  !$omp critical
  r_loc = r_loc + r_priv
  !$omp end critical
  !$omp barrier

  ! global result ..............................................................

  if (present(comm)) then

    !$omp master
    call MPI_Allreduce(r_loc, r_glob, 1, MPI_REAL_RNP, MPI_SUM, comm)
    !$omp end master
    !$omp barrier

    r = r_glob

  else

    r = r_loc

  end if

end function WeightedScalarProduct_X

!===============================================================================

end module Array_Reductions
