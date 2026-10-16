!> summary:   3D interpolation operator for Spectral-Element / DG methods
!> author:    Immo Huismann
!> date:      2016/11/22
!> license:   Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

program Potpourri_3
  use Kind_Parameters, only: RDP
  use Constants,       only: ONE, ZERO
  use Potpourri_Print_MLUPS, only: Print_MLUPS
  use MPI_F08
  implicit none

  integer, parameter :: n            = 100000
  integer, parameter :: n_dof        = n
  integer, parameter :: n_repetition = 100

  real(RDP), allocatable :: u(:), v(:), w(:)
  real(RDP) :: start_time, stop_time
  integer   :: i_repetition
  real(RDP) :: h

  !.............................................................................
  ! Allocation of data

  allocate(u(n), v(n), w(n))

  !.............................................................................
  ! Initialization

  u = 7
  v = [((-1)**i_repetition, i_repetition = 1, n)]
  h = ONE / (n - 1)

  call MPI_Init()

  !.............................................................................
  ! Test with baseline variant

  start_time = MPI_Wtime()

  do i_repetition = 1, n_repetition
    call Baseline(n, h, u, v, w)
  end do

  stop_time = MPI_Wtime()

  print *
  print '(A)', 'Baseline variant'
  call Print_MLUPS(n_dof, n_repetition, start_time, stop_time)

  !.............................................................................
  ! Test with improved variant

  start_time = MPI_Wtime()

  do i_repetition = 1, n_repetition
    call Improved(n, h, u, v, w)
  end do

  stop_time = MPI_Wtime()

  print *
  print '(A)', 'Improved variant'
  call Print_MLUPS(n_dof, n_repetition, start_time, stop_time)

  !.............................................................................
  ! Cleanup

  call MPI_Finalize()

contains

!-------------------------------------------------------------------------------
!> 1D upwind scheme for convection -- Baseline

subroutine Baseline(n, h, u, v, w)
  integer,   intent(in)  :: n    !< points per direction
  real(RDP), intent(in)  :: h    !< mesh width
  real(RDP), intent(in)  :: u(n) !< values to convect
  real(RDP), intent(in)  :: v(n) !< convection velocity
  real(RDP), intent(out) :: w(n) !< convection term

  integer :: i

  do i = 2, n-1

    if(v(i) > 0) then
      w(i) = (u(i  ) - u(i-1)) / h
    else
      w(i) = (u(i+1) - u(i  )) / h
    end if

  end do

end subroutine Baseline

!-------------------------------------------------------------------------------
!> 1D upwind scheme for convection -- Improved

subroutine Improved(n, h, u, v, w)
  integer,   intent(in)  :: n    !< points per direction
  real(RDP), intent(in)  :: h    !< mesh width
  real(RDP), intent(in)  :: u(n) !< values to convect
  real(RDP), intent(in)  :: v(n) !< convection velocity
  real(RDP), intent(out) :: w(n) !< convection term

  integer :: i

  ! TODO

end subroutine Improved

!===============================================================================

end program Potpourri_3
