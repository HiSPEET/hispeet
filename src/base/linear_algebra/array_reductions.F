!> summary:  Parallel accelerated array reductions
!> author:   Joerg Stiller
!> date:     2017/06/28
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
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
    module procedure SumOfElements_2
    module procedure SumOfElements_3
    module procedure SumOfElements_4
    module procedure SumOfElements_5
  end interface

  interface ScalarProduct
    module procedure ScalarProduct_2
    module procedure ScalarProduct_3
    module procedure ScalarProduct_4
    module procedure ScalarProduct_5
  end interface

  interface WeightedScalarProduct
    module procedure WeightedScalarProduct_2
    module procedure WeightedScalarProduct_3
    module procedure WeightedScalarProduct_4
    module procedure WeightedScalarProduct_5
  end interface

contains

  !=============================================================================
  ! SumOfElements

  !-----------------------------------------------------------------------------
  !> Sum over all array elements (2D)
  !>
  !> If the MPI communicator comm is passed, the sum will be evaluated over all
  !> processes. Otherwise, only the local contribution is returned.

  real(RNP) function SumOfElements_2(a, comm) result(r)
    real(RNP),                intent(in) :: a(:,:) !< array
    type(MPI_Comm), optional, intent(in) :: comm   !< MPI communicator

    r = SumOfElements_X(size(a), a, comm)

  end function SumOfElements_2

  !-----------------------------------------------------------------------------
  !> Sum over all array elements (3D)
  !>
  !> If the MPI communicator comm is passed, the sum will be evaluated over all
  !> processes. Otherwise, only the local contribution is returned.

  real(RNP) function SumOfElements_3(a, comm) result(r)
    real(RNP),                intent(in) :: a(:,:,:) !< array
    type(MPI_Comm), optional, intent(in) :: comm     !< MPI communicator

    r = SumOfElements_X(size(a), a, comm)

  end function SumOfElements_3

  !-----------------------------------------------------------------------------
  !> Sum over all array elements (4D)
  !>
  !> If the MPI communicator comm is passed, the sum will be evaluated over all
  !> processes. Otherwise, only the local contribution is returned.

  real(RNP) function SumOfElements_4(a, comm) result(r)
    real(RNP),                intent(in) :: a(:,:,:,:) !< array
    type(MPI_Comm), optional, intent(in) :: comm       !< MPI communicator

    r = SumOfElements_X(size(a), a, comm)

  end function SumOfElements_4

  !-----------------------------------------------------------------------------
  !> Sum over all array elements (5D)
  !>
  !> If the MPI communicator comm is passed, the sum will be evaluated over all
  !> processes. Otherwise, only the local contribution is returned.

  real(RNP) function SumOfElements_5(a, comm) result(r)
    real(RNP),                intent(in) :: a(:,:,:,:,:) !< array
    type(MPI_Comm), optional, intent(in) :: comm         !< MPI communicator

    r = SumOfElements_X(size(a), a, comm)

  end function SumOfElements_5

  !-----------------------------------------------------------------------------
  !> Sum over all array elements (eXplicit)
  !>
  !> If the MPI communicator comm is passed, the sum will be evaluated over all
  !> processes. Otherwise, only the local contribution is returned.

  real(RNP) function SumOfElements_X(n, a, comm) result(r)
    integer,                  intent(in) :: n    !< array size
    real(RNP),                intent(in) :: a(n) !< array
    type(MPI_Comm), optional, intent(in) :: comm !< MPI communicator

    ! local variables declared save to become shared with OpenMP
    real(RNP), save :: r_loc = 0, r_glob = 0
    integer  :: i

    ! local contribution .......................................................

    !$acc data present(a)
    !$acc parallel loop reduction(+:r_loc)
    !$omp do reduction(+:r_loc)
    do i = 1, n
      r_loc = r_loc + a(i)
    end do
    !$acc end data

    ! OpenACC:
    ! Note that the reduction variable is automatically updated on the host!

    ! global result ............................................................

    !$omp master
    if (present(comm)) then
      call MPI_Allreduce(r_loc, r_glob, 1, MPI_REAL_RNP, MPI_SUM, comm)
    else
      r_glob = r_loc
    end if
    r_loc = 0
    !$omp end master
    !$omp barrier

    r = r_glob

  end function SumOfElements_X

  !=============================================================================
  ! ScalarProduct

  !-----------------------------------------------------------------------------
  !> Scalar product of vectors a and b (2D)
  !>
  !> With the MPI communicator comm passed, the scalar product is be evaluated
  !> over all processes. Otherwise, only the local contribution is returned.

  real(RNP) function ScalarProduct_2(a, b, comm) result(r)
    real(RNP),                intent(in) :: a(:,:) !< left operand
    real(RNP),                intent(in) :: b(:,:) !< right operand
    type(MPI_Comm), optional, intent(in) :: comm   !< MPI communicator

    r = ScalarProduct_X(size(a), a, b, comm)

  end function ScalarProduct_2

  !-----------------------------------------------------------------------------
  !> Scalar product of vectors a and b (3D)
  !>
  !> With the MPI communicator comm passed, the scalar product is be evaluated
  !> over all processes. Otherwise, only the local contribution is returned.

  real(RNP) function ScalarProduct_3(a, b, comm) result(r)
    real(RNP),                intent(in) :: a(:,:,:) !< left operand
    real(RNP),                intent(in) :: b(:,:,:) !< right operand
    type(MPI_Comm), optional, intent(in) :: comm     !< MPI communicator

    r = ScalarProduct_X(size(a), a, b, comm)

  end function ScalarProduct_3

  !-----------------------------------------------------------------------------
  !> Scalar product of vectors a and b (4D)
  !>
  !> With the MPI communicator comm passed, the scalar product is be evaluated
  !> over all processes. Otherwise, only the local contribution is returned.

  real(RNP) function ScalarProduct_4(a, b, comm) result(r)
    real(RNP),                intent(in) :: a(:,:,:,:) !< left operand
    real(RNP),                intent(in) :: b(:,:,:,:) !< right operand
    type(MPI_Comm), optional, intent(in) :: comm       !< MPI communicator

    r = ScalarProduct_X(size(a), a, b, comm)

  end function ScalarProduct_4

  !-----------------------------------------------------------------------------
  !> Scalar product of vectors a and b (5D)
  !>
  !> With the MPI communicator comm passed, the scalar product is be evaluated
  !> over all processes. Otherwise, only the local contribution is returned.

  real(RNP) function ScalarProduct_5(a, b, comm) result(r)
    real(RNP),                intent(in) :: a(:,:,:,:,:) !< left operand
    real(RNP),                intent(in) :: b(:,:,:,:,:) !< right operand
    type(MPI_Comm), optional, intent(in) :: comm         !< MPI communicator

    r = ScalarProduct_X(size(a), a, b, comm)

  end function ScalarProduct_5

  !-----------------------------------------------------------------------------
  !> Scalar product of vectors a and b (eXplicit)
  !>
  !> With the MPI communicator comm passed, the scalar product is be evaluated
  !> over all processes. Otherwise, only the local contribution is returned.

#if defined(OMP_REDUCTION_SERIAL)

  ! omp serial version
  real(RNP) function ScalarProduct_X(n, a, b, comm) result(r)
    integer,                  intent(in) :: n    !< array size
    real(RNP),                intent(in) :: a(n) !< left operand
    real(RNP),                intent(in) :: b(n) !< right operand
    type(MPI_Comm), optional, intent(in) :: comm !< MPI communicator

    real(RNP), save :: r_glob
    real(RNP) :: r_loc
    integer :: i

    !$omp master

    r_loc = 0
    do i = 1, n
      r_loc = r_loc + a(i) * b(i)
    end do

    if (present(comm)) then
      call MPI_Allreduce(r_loc, r_glob, 1, MPI_REAL_RNP, MPI_SUM, comm)
    else
      r_glob = r_loc
    end if

    !$omp end master
    !$omp barrier

    r = r_glob

  end function ScalarProduct_X

#elif defined(OMP_REDUCTION_ATOMIC)

  ! omp atomic version
  real(RNP) function ScalarProduct_X(n, a, b, comm) result(r)
    integer,                  intent(in) :: n    !< array size
    real(RNP),                intent(in) :: a(n) !< left operand
    real(RNP),                intent(in) :: b(n) !< right operand
    type(MPI_Comm), optional, intent(in) :: comm !< MPI communicator

    ! local variables declared save to become shared with OpenMP
    real(RNP), save :: r_loc, r_glob
    real(RNP) :: r_priv
    integer :: i

    ! local contribution .......................................................

    ! skipping initialization with _OPENMP defined failed in some cases
    !$omp single
    r_loc = 0
    !$omp end single

    r_priv = 0

    !$omp do private(i)
    do i = 1, n
      r_priv = r_priv + a(i) * b(i)
    end do
    !$omp end do nowait

    !$omp atomic
    r_loc = r_loc + r_priv
    !$omp barrier

    ! global result ............................................................

    !$omp master
    if (present(comm)) then
      call MPI_Allreduce(r_loc, r_glob, 1, MPI_REAL_RNP, MPI_SUM, comm)
    else
      r_glob = r_loc
    end if
    !$omp end master
    !$omp barrier

    r = r_glob

  end function ScalarProduct_X

#else

  ! default version using omp reduction
  real(RNP) function ScalarProduct_X(n, a, b, comm) result(r)
    integer,                  intent(in) :: n    !< array size
    real(RNP),                intent(in) :: a(n) !< left operand
    real(RNP),                intent(in) :: b(n) !< right operand
    type(MPI_Comm), optional, intent(in) :: comm !< MPI communicator

    ! local variables declared save to become shared with OpenMP
    real(RNP), save :: r_loc = 0, r_glob = 0
    integer :: i

    ! local contribution .......................................................

    !$omp do reduction(+:r_loc) private(i)
    do i = 1, n
      r_loc = r_loc + a(i) * b(i)
    end do

    ! global result ............................................................

    !$omp master
    if (present(comm)) then
      call MPI_Allreduce(r_loc, r_glob, 1, MPI_REAL_RNP, MPI_SUM, comm)
    else
      r_glob = r_loc
    end if
    r_loc = 0
    !$omp end master
    !$omp barrier

    r = r_glob

  end function ScalarProduct_X

#endif

  !=============================================================================
  ! WeightedScalarProduct

  !-----------------------------------------------------------------------------
  !> Weighted scalar product of vectors a and b (2D)
  !>
  !> With the MPI communicator comm passed, the scalar product is be evaluated
  !> over all processes. Otherwise, only the local contribution is returned.

  real(RNP) function WeightedScalarProduct_2(w, a, b, comm) result(r)
    real(RNP),                intent(in) :: w(:,:) !< weights
    real(RNP),                intent(in) :: a(:,:) !< left operand
    real(RNP),                intent(in) :: b(:,:) !< right operand
    type(MPI_Comm), optional, intent(in) :: comm   !< MPI communicator

    r = WeightedScalarProduct_X(size(w), w, a, b, comm)

  end function WeightedScalarProduct_2

  !-----------------------------------------------------------------------------
  !> Weighted scalar product of vectors a and b (3D)
  !>
  !> With the MPI communicator comm passed, the scalar product is be evaluated
  !> over all processes. Otherwise, only the local contribution is returned.

  real(RNP) function WeightedScalarProduct_3(w, a, b, comm) result(r)
    real(RNP),                intent(in) :: w(:,:,:) !< weights
    real(RNP),                intent(in) :: a(:,:,:) !< left operand
    real(RNP),                intent(in) :: b(:,:,:) !< right operand
    type(MPI_Comm), optional, intent(in) :: comm     !< MPI communicator

    r = WeightedScalarProduct_X(size(w), w, a, b, comm)

  end function WeightedScalarProduct_3

  !-----------------------------------------------------------------------------
  !> Weighted scalar product of vectors a and b (4D)
  !>
  !> With the MPI communicator comm passed, the scalar product is be evaluated
  !> over all processes. Otherwise, only the local contribution is returned.

  real(RNP) function WeightedScalarProduct_4(w, a, b, comm) result(r)
    real(RNP),                intent(in) :: w(:,:,:,:) !< weights
    real(RNP),                intent(in) :: a(:,:,:,:) !< left operand
    real(RNP),                intent(in) :: b(:,:,:,:) !< right operand
    type(MPI_Comm), optional, intent(in) :: comm       !< MPI communicator

    r = WeightedScalarProduct_X(size(w), w, a, b, comm)

  end function WeightedScalarProduct_4

  !-----------------------------------------------------------------------------
  !> Weighted scalar productof vectors a and b (5D)
  !>
  !> With the MPI communicator comm passed, the scalar product is be evaluated
  !> over all processes. Otherwise, only the local contribution is returned.

  real(RNP) function WeightedScalarProduct_5(w, a, b, comm) result(r)
    real(RNP),                intent(in) :: w(:,:,:,:,:) !< weights
    real(RNP),                intent(in) :: a(:,:,:,:,:) !< left operand
    real(RNP),                intent(in) :: b(:,:,:,:,:) !< right operand
    type(MPI_Comm), optional, intent(in) :: comm         !< MPI communicator

    r = WeightedScalarProduct_X(size(w), w, a, b, comm)

  end function WeightedScalarProduct_5

  !-----------------------------------------------------------------------------
  !> Weighted scalar productof vectors a and b (eXplicit)
  !>
  !> With the MPI communicator comm passed, the scalar product is be evaluated
  !> over all processes. Otherwise, only the local contribution is returned.

#if defined(OMP_REDUCTION_SERIAL)

  ! omp serial version

  real(RNP) function WeightedScalarProduct_X(n, w, a, b, comm) result(r)
    integer,                  intent(in) :: n    !< array size
    real(RNP),                intent(in) :: w(n) !< weights
    real(RNP),                intent(in) :: a(n) !< left operand
    real(RNP),                intent(in) :: b(n) !< right operand
    type(MPI_Comm), optional, intent(in) :: comm !< MPI communicator

    ! local variables declared save to become shared with OpenMP
    real(RNP), save :: r_glob
    real(RNP) :: r_loc
    integer  :: i

    !$omp master

    r_loc = 0
    do i = 1, n
      r_loc = r_loc + w(i) * a(i) * b(i)
    end do

    if (present(comm)) then
      call MPI_Allreduce(r_loc, r_glob, 1, MPI_REAL_RNP, MPI_SUM, comm)
    else
      r_glob = r_loc
    end if

    !$omp end master
    !$omp barrier

    r = r_glob

  end function WeightedScalarProduct_X

#elif defined(OMP_REDUCTION_ATOMIC)

  ! omp atomic version
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

    ! local contribution .......................................................

    ! skipping initialization with _OPENMP defined failed in some cases
    !$omp single
    r_loc = 0
    !$omp end single

    r_priv = 0

    !$omp do private(i)
    do i = 1, n
      r_priv = r_priv + w(i) * a(i) * b(i)
    end do
    !$omp end do nowait

    !$omp atomic
    r_loc = r_loc + r_priv
    !$omp barrier

    ! global result ............................................................

    !$omp master
    if (present(comm)) then
      call MPI_Allreduce(r_loc, r_glob, 1, MPI_REAL_RNP, MPI_SUM, comm)
    else
      r_glob = r_loc
    end if
    !$omp end master
    !$omp barrier

    r = r_glob

  end function WeightedScalarProduct_X

#else

  ! default version using omp reduction
  real(RNP) function WeightedScalarProduct_X(n, w, a, b, comm) result(r)
    integer,                  intent(in) :: n    !< array size
    real(RNP),                intent(in) :: w(n) !< weights
    real(RNP),                intent(in) :: a(n) !< left operand
    real(RNP),                intent(in) :: b(n) !< right operand
    type(MPI_Comm), optional, intent(in) :: comm !< MPI communicator

    ! local variables declared save to become shared with OpenMP
    real(RNP), save :: r_loc = 0, r_glob = 0
    integer  :: i

    ! local contribution .......................................................

    !$omp do reduction(+:r_loc) private(i)
    do i = 1, n
      r_loc = r_loc + w(i) * a(i) * b(i)
    end do

    ! global result ............................................................

    !$omp master
    if (present(comm)) then
      call MPI_Allreduce(r_loc, r_glob, 1, MPI_REAL_RNP, MPI_SUM, comm)
    else
      r_glob = r_loc
    end if
    r_loc = 0
    !$omp end master
    !$omp barrier

    r = r_glob

  end function WeightedScalarProduct_X

#endif

!===============================================================================

end module Array_Reductions
