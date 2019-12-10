!> \file      poutpouri_a1.f
!> \brief     A small test program
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



program Poutpouri_A1
  use Kind_Parameters,       only: RDP

  use Potpourri_Print_MLUPS, only: PrintMLUPS
  use MPI_F08
  implicit none

  integer, parameter :: N            = 2000
  integer, parameter :: N_DOF        = N**2
  integer, parameter :: N_REPITITION = 100

  real(RDP), allocatable :: u(:,:), v(:,:)
  real(RDP) :: start_time, stop_time
  integer   :: i_repitition

  !.............................................................................
  ! Allocation of data

  allocate(u(N,N), v(N,N))

  !.............................................................................
  ! Initialization

  u = 7
  v = 0
  call MPI_Init()

  !.............................................................................
  ! Test with baseline variant

  start_time = MPI_Wtime()

  do i_repitition = 1, N_REPITITION
    call Baseline(N,u,v)
  end do

  stop_time            = MPI_Wtime()

  call PrintMLUPS(n_dof,n_repitition,start_time,stop_time)

  !.............................................................................
  ! Test with improved variant

  start_time = MPI_Wtime()

  do i_repitition = 1, N_REPITITION
    call Improved(N,u,v)
  end do

  stop_time            = MPI_Wtime()

  call PrintMLUPS(n_dof,n_repitition,start_time,stop_time)

  !.............................................................................
  ! Cleanup

  call MPI_Finalize()

contains

!-------------------------------------------------------------------------------
!> \brief   Copy operator -- Baseline
!> \author  Immo Huismann

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
!> \brief   Copy operator -- Improved
!> \author  Immo Huismann

subroutine Improved(n_point,u,v)
  integer,   intent(in)  :: n_point
  real(RDP), intent(in)  :: u(n_point,n_point)
  real(RDP), intent(out) :: v(n_point,n_point)

end subroutine Improved

!===============================================================================

end program Poutpouri_A1
