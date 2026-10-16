!> summary:   A small test program
!> author:    Immo Huismann
!> date:      2016/11/22
!> license:   Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Potpourri_Print_MLUPS
  use Kind_Parameters, only: RDP
  implicit none
  private

  public :: Print_MLUPS

contains

!-------------------------------------------------------------------------------
!> Prints the  Mega Lattice Updates Per Second

subroutine Print_MLUPS(n_dof, n_repititions, start_time, stop_time)
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

end subroutine Print_MLUPS

!===============================================================================

end module Potpourri_Print_MLUPS
