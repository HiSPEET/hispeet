program Test_Scalar_Product
  !$ use OMP_Lib                ! only compiled with OpenMP enabled
  use Kind_Parameters
  use MPI_Binding
  use Array_Reductions
  implicit none

  type(MPI_Comm) :: comm
  real(RNP), allocatable, save :: a(:), b(:)
  real(RNP) :: ab, ab_omp, ab_ref
  real(RNP) :: t_omp, t_ref, t0
  integer   :: i_proc, n_proc, n_thread

  integer :: n_vect = 120000000   ! default vector length
  integer :: n_test = 100         ! default number of test repetitions
  integer :: ios, k, n
  character(len=64) :: arg

  ! initialization .............................................................

  call Init_MPI_Binding()

  comm = MPI_COMM_WORLD
  call MPI_Comm_rank(comm, i_proc)
  call MPI_Comm_size(comm, n_proc)

  ! command line arguments (optional) .........................................
  !   1st: vector length n_vect, adjusted to the nearest multiple of n_proc
  !   2nd: number of test repetitions n_test

  if (i_proc == 0) then

    if (command_argument_count() >= 1) then
      call get_command_argument(1, arg)
      read(arg, *, iostat=ios) n_vect
      if (ios /= 0 .or. n_vect < 1) then
        write(*,'(3A)') 'Error: invalid vector length "', trim(arg), '"'
        call MPI_Abort(comm, 1)
      end if
    end if

    if (command_argument_count() >= 2) then
      call get_command_argument(2, arg)
      read(arg, *, iostat=ios) n_test
      if (ios /= 0 .or. n_test < 1) then
        write(*,'(3A)') 'Error: invalid number of tests "', trim(arg), '"'
        call MPI_Abort(comm, 1)
      end if
    end if

    n = n_vect
    n_vect = max(1, nint(n / real(n_proc, RNP))) * n_proc
    if (n_vect /= n) then
      write(*,'(/,4(A,G0))') 'Note: n adjusted from ', n, ' to ', n_vect, &
                             ' (multiple of n_proc = ', n_proc, ')'
    end if

  end if

  call MPI_Bcast(n_vect, 1, MPI_INTEGER, 0, comm)
  call MPI_Bcast(n_test, 1, MPI_INTEGER, 0, comm)

  n_thread = 1
  !$omp parallel
  !$omp master
  !$ n_thread = omp_get_num_threads()
  !$omp end master
  !$omp end parallel

  if (i_proc == 0) then
    write(*,'(/,A)')  'Testing scalar product with hybrid parallelization'
    write(*,'(A,G0)') '  n_proc   = ', n_proc
    write(*,'(A,G0)') '  n_thread = ', n_thread
    write(*,'(A,G0)') '  n_vector = ', n_vect
    write(*,'(A,G0)') '  n_test   = ', n_test
  end if

  n = n_vect / n_proc
  allocate(a(n), b(n))
  call random_number(a)
  call random_number(b)

  ! reference: mean wall time of n_test repetitions ...........................

  call MPI_Barrier(comm)
  t0 = MPI_Wtime()
  do k = 1, n_test
    ab = sum(a*b)
    call MPI_Allreduce(ab, ab_ref, 1, MPI_REAL_RNP, MPI_SUM, comm)
  end do
  t_ref = (MPI_Wtime() - t0) / n_test

  ! function using OpenMP and MPI: mean wall time of n_test repetitions ........

  !$omp parallel private(ab, k)

  !$omp master
  call MPI_Barrier(comm)
  t0 = MPI_Wtime()
  !$omp end master
  !$omp barrier

  do k = 1, n_test
    ab = ScalarProduct(a, b, comm)
  end do

  !$omp master
  ab_omp = ab
  t_omp = (MPI_Wtime() - t0) / n_test
  !$omp end master

  !$omp end parallel

  ! results ....................................................................

  if (i_proc == 0) then
    write(*,'(/,A)')  'Test results'
    write(*,'(2(A,ES12.5))') '  reference:  ab = ', ab_ref,', t_wall = ', t_ref
    write(*,'(2(A,ES12.5))') '  hybrid:     ab = ', ab_omp,', t_wall = ', t_omp
  end if

  call MPI_Finalize()

end program Test_Scalar_Product

