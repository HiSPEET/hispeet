!> \file      poutpouri_a2.f
!> \brief     3D interpolation operator for Spectral-Element / DG methods
!> \author    Immo Huismann
!> \date      2016/11/22
!> \copyright Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> \details
!>
!>
!===============================================================================
module Potpourri_Print_MLUPS
  use Kind_Parameters, only: RDP
  implicit none
  private

  public :: PrintMLUPS

contains

!-------------------------------------------------------------------------------
!> \brief   Prints the  Mega Lattice Updates Per Second
!> \author  Immo Huismann

subroutine PrintMLUPS(n_dof, n_repititions, start_time, stop_time)
  integer,   intent(in) :: n_dof         !< number of degrees of freedom
  integer,   intent(in) :: n_repititions !< number of repititions of test
  real(RDP), intent(in) :: start_time    !< start of time measurement
  real(RDP), intent(in) :: stop_time     !< end   of time measurement

  real(RDP) :: time_taken, time_per_application, mlups
  character(len=*), parameter :: USED_FORMAT = '(A,1X,E15.7)'
  
  time_taken           = stop_time - start_time
  time_per_application = time_taken / n_repititions
  mlups                = n_dof / time_per_application / 10**6
  
  print USED_FORMAT, 'MLUPS:', mlups

end subroutine PrintMLUPS

!===============================================================================

end module Potpourri_Print_MLUPS


module Potpourri_A2_Derivatives
  use Kind_Parameters, only: RDP
  implicit none
  private

  public ::  FirstDerivative
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

end module Potpourri_A2_Derivatives


program Poutpouri_A2
  use Kind_Parameters,          only: RDP
  use Constants,                only: ONE

  use Potpourri_Print_MLUPS,    only: PrintMLUPS
  use Potpourri_A2_Derivatives, only: FirstDerivative, SecondDerivative
  use MPI_F08
  implicit none

  integer,   parameter :: N            = 10000
  integer,   parameter :: N_DOF        = N
  integer,   parameter :: N_REPITITION = 100
  real(RDP), parameter :: H = ONE / (N - 1)

  real(RDP), allocatable :: u(:), v(:), w(:)
  real(RDP) :: start_time, stop_time
  integer   :: i_repitition

  !.............................................................................
  ! Allocation of data

  allocate(u(N), v(N), w(N))

  !.............................................................................
  ! Initialization

  u = 7
  v = 0

  call MPI_Init()

  !.............................................................................
  ! Test with baseline variant

  start_time = MPI_Wtime()

  do i_repitition = 1, N_REPITITION
    call Baseline(N,H,u,v,w)
  end do

  stop_time            = MPI_Wtime()

  call PrintMLUPS(N_DOF,N_REPITITION,start_time,stop_time)

  !.............................................................................
  ! Test with improved variant

  start_time = MPI_Wtime()

  do i_repitition = 1, N_REPITITION
    call Improved(N,H,u,v,w)
  end do

  stop_time            = MPI_Wtime()

  call PrintMLUPS(N_DOF,N_REPITITION,start_time,stop_time)

  !.............................................................................
  ! Cleanup

  call MPI_Finalize()

contains

!-------------------------------------------------------------------------------
!> \brief   Convection diffusion term for diffusivity lambda = 1 -- Baseline
!> \author  Immo Huismann

subroutine Baseline(n,h,u,v,w)
  integer,   intent(in)  :: n    !< points per direction
  real(RDP), intent(in)  :: h    !< mesh width
  real(RDP), intent(in)  :: u(n) !< variable do convect / diffuse
  real(RDP), intent(in)  :: v(n) !< velocity
  real(RDP), intent(out) :: w(n) !< convection term

  integer :: i

  do i = 2, n-1

    w(i) =  FirstDerivative(h,u(i-1),     u(i+1)) * v(i)                       &
         + SecondDerivative(h,u(i-1),u(i),u(i+1))

  end do

end subroutine Baseline

!-------------------------------------------------------------------------------
!> \brief   Convection diffusion term for diffusivity lambda = 1 -- Improved
!> \author  Immo Huismann

subroutine Improved(n,h,u,v,w)
  integer,   intent(in)  :: n    !< points per direction
  real(RDP), intent(in)  :: h    !< mesh width
  real(RDP), intent(in)  :: u(n) !< variable do convect / diffuse
  real(RDP), intent(in)  :: v(n) !< velocity
  real(RDP), intent(out) :: w(n) !< convection term

end subroutine Improved

!===============================================================================

end program Poutpouri_A2
