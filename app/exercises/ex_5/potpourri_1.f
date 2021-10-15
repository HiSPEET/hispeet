program Potpourri_1
  use Kind_Parameters,       only: RDP
  use Potpourri_Print_MLUPS, only: Print_MLUPS
  use MPI_F08
  implicit none

  integer, parameter :: n            = 2000
  integer, parameter :: n_dof        = n**2
  integer, parameter :: n_repetition = 100

  real(RDP), allocatable :: u(:,:), v(:,:)
  real(RDP) :: start_time, stop_time
  integer   :: i_repetition

  !.............................................................................
  ! Allocation of data

  allocate(u(n,n), v(n,n))

  !.............................................................................
  ! Initialization

  u = 7
  v = 0
  call MPI_Init()

  !.............................................................................
  ! Test with baseline variant

  start_time = MPI_Wtime()

  do i_repetition = 1, n_repetition
    call Baseline(n, u, v)
  end do

  stop_time = MPI_Wtime()

  call Print_MLUPS(n_dof, n_repetition, start_time, stop_time)

  !.............................................................................
  ! Test with improved variant

  start_time = MPI_Wtime()

  do i_repetition = 1, n_repetition
    call Improved(n, u, v)
  end do

  stop_time = MPI_Wtime()

  call Print_MLUPS(n_dof, n_repetition, start_time, stop_time)

  !.............................................................................
  ! Cleanup

  call MPI_Finalize()

contains

!-------------------------------------------------------------------------------
!> Copy operator -- Baseline

subroutine Baseline(n_point,u,v)
  integer,   intent(in)  :: n_point
  real(RDP), intent(in)  :: u(n_point,n_point)
  real(RDP), intent(out) :: v(n_point,n_point)

  integer :: i,j

  do i = 1, n_point
  do j = 1, n_point
  
    v(i,j) = u(i,j)
    
  end do
  end do

end subroutine Baseline

!-------------------------------------------------------------------------------
!> Copy operator -- Improved

subroutine Improved(n_point,u,v)
  integer,   intent(in)  :: n_point
  real(RDP), intent(in)  :: u(n_point,n_point)
  real(RDP), intent(out) :: v(n_point,n_point)
  
  ! TODO

end subroutine Improved

!===============================================================================

end program Potpourri_1
