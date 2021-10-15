!> summary:   3D interpolation operator for Spectral-Element / DG methods
!> author:    Immo Huismann
!> date:      2016/11/22
!> license: Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany

!===============================================================================

module Potpourri__Derivatives
  use Kind_Parameters, only: RDP
  implicit none
  private

  public :: FirstDerivative
  public :: SecondDerivative

contains

pure function FirstDerivative(h,u_left,u_right) result(du)
  real(RDP), intent(in) :: h
  real(RDP), intent(in) :: u_left
  real(RDP), intent(in) :: u_right
  real(RDP)             :: du

  du = (u_right - u_left) / (2 * h)

end function FirstDerivative

pure function SecondDerivative(h,u_left,u_middle,u_right) result(du)
  real(RDP), intent(in) :: h
  real(RDP), intent(in) :: u_left
  real(RDP), intent(in) :: u_middle
  real(RDP), intent(in) :: u_right
  real(RDP)             :: du

  du = (u_right - 2 * u_middle + u_left) / (h**2)

end function SecondDerivative

!===============================================================================

end module Potpourri__Derivatives

!-------------------------------------------------------------------------------

program Potpourri_2
  use Kind_Parameters,          only: RDP
  use Constants,                only: ONE
  use Potpourri_Print_MLUPS,    only: Print_MLUPS
  use Potpourri__Derivatives,   only: FirstDerivative, SecondDerivative
  use MPI_F08
  implicit none

  integer,   parameter :: n            = 10000
  integer,   parameter :: n_dof        = n
  integer,   parameter :: n_repetition = 100
  real(RDP), parameter :: h = ONE / (n - 1)

  real(RDP), allocatable :: u(:), v(:), w(:)
  real(RDP) :: start_time, stop_time
  integer   :: i_repetition

  !.............................................................................
  ! Allocation of data

  allocate(u(n), v(n), w(n))

  !.............................................................................
  ! Initialization

  u = 7
  v = 0

  call MPI_Init()

  !.............................................................................
  ! Test with baseline variant

  start_time = MPI_Wtime()

  do i_repetition = 1, n_repetition
    call Baseline(n, h, u, v, w)
  end do

  stop_time = MPI_Wtime()

  call Print_MLUPS(n_dof, n_repetition, start_time, stop_time)

  !.............................................................................
  ! Test with improved variant

  start_time = MPI_Wtime()

  do i_repetition = 1, n_repetition
    call Improved(n, h, u, v, w)
  end do

  stop_time = MPI_Wtime()

  call Print_MLUPS(n_dof, n_repetition, start_time, stop_time)

  !.............................................................................
  ! Cleanup

  call MPI_Finalize()

contains

!-------------------------------------------------------------------------------
!> Convection diffusion term for diffusivity lambda = 1 -- Baseline

subroutine Baseline(n, h, u, v, w)
  integer,   intent(in)  :: n    !< points per direction
  real(RDP), intent(in)  :: h    !< mesh width
  real(RDP), intent(in)  :: u(n) !< variable do convect / diffuse
  real(RDP), intent(in)  :: v(n) !< velocity
  real(RDP), intent(out) :: w(n) !< convection term

  integer :: i

  do i = 2, n-1

    w(i) =  FirstDerivative(h,u(i-1),     u(i+1)) * v(i)       &
         + SecondDerivative(h,u(i-1),u(i),u(i+1))

  end do

end subroutine Baseline

!-------------------------------------------------------------------------------
!> Convection diffusion term for diffusivity lambda = 1 -- Improved

subroutine Improved(n, h, u, v, w)
  integer,   intent(in)  :: n    !< points per direction
  real(RDP), intent(in)  :: h    !< mesh width
  real(RDP), intent(in)  :: u(n) !< variable do convect / diffuse
  real(RDP), intent(in)  :: v(n) !< velocity
  real(RDP), intent(out) :: w(n) !< convection term
  
  ! TODO

end subroutine Improved

!===============================================================================

end program Potpourri_2
