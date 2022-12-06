!> summary:   A small program running a stencil operation
!> author:    Immo Huismann
!> date:      2016/11/21
!> license:   Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Stencil_Baseline
  use Kind_Parameters, only: RDP
  implicit none
  private

  public :: StencilBaseline

contains

!-------------------------------------------------------------------------------
!> A baseline version of the five point stencil

subroutine StencilBaseline(n,h,a,b)
  integer,   intent(in)  :: n
  real(RDP), intent(in)  :: h
  real(RDP), intent(in)  :: a(n,n)
  real(RDP), intent(out) :: b(n,n)

  integer :: i,j

  do j = 2, n - 1
  do i = 2, n - 1
    b(i,j) = (a(i-1,j) + a(i+1,j) + a(i,j-1) + a(i,j+1)) / h**2
  end do
  end do
      

end subroutine StencilBaseline

!===============================================================================

end module Stencil_Baseline


module Stencil_Direct
  use Kind_Parameters, only: RDP
  implicit none
  private

  public :: StencilDirect

contains

!-------------------------------------------------------------------------------
!> A direct version of the five point stencil

subroutine StencilDirect(n,h,a,b)
  integer,   intent(in)  :: n
  real(RDP), intent(in)  :: h
  real(RDP), intent(in)  :: a(n,n)
  real(RDP), intent(out) :: b(n,n)

  integer   :: i,j
  real(RDP) :: s
  
  s = 1 / h**2

  do j = 2, n - 1
  do i = 2, n - 1
    b(i,j) = (a(i-1,j) + a(i+1,j) + a(i,j-1) + a(i,j+1)) * s
  end do
  end do

end subroutine StencilDirect

!===============================================================================

end module Stencil_Direct



module Stencil_Blocked
  use Kind_Parameters, only: RDP
  implicit none
  private

  public :: StencilBlocked

contains

!-------------------------------------------------------------------------------
!> A blocked version of the five point stencil

subroutine StencilBlocked(n,h,a,b)
  integer,   intent(in)  :: n
  real(RDP), intent(in)  :: h
  real(RDP), intent(in)  :: a(n,n)
  real(RDP), intent(out) :: b(n,n)

  integer   :: i,j
  real(RDP) :: s

  s = 1 / h**2

  do j = 2, n - 1
  do i = 2, n - 1
    b(i,j) = (a(i-1,j) + a(i+1,j) + a(i,j-1) + a(i,j+1)) * s
  end do
  end do

end subroutine StencilBlocked

!===============================================================================

end module Stencil_Blocked

!-------------------------------------------------------------------------------

program Stencil_Test
  use Kind_Parameters,  only: RDP
  use Constants,        only: ZERO, ONE, PI

  use Stencil_Baseline, only: StencilBaseline
  use Stencil_Direct,   only: StencilDirect
  use Stencil_Blocked,  only: StencilBlocked  

  use MPI_F08
  implicit none

  integer,   parameter   :: n = 200             ! number of intervals
  real(RDP), parameter   :: l = ONE             ! domain length
  real(RDP), parameter   :: h = l / (n - 1)     ! interval width
  integer,   parameter   :: n_repetition = 1000 ! number of repititions for test  

  real(RDP), allocatable :: u(:,:), v(:,:), x(:,:), y(:,:)
  integer :: i, j

  integer :: i_repetition

  real(RDP) :: start_time, stop_time, time_taken, time_per_application

  !.............................................................................
  ! Allocation of data
  
  allocate(u(n,n), source = ZERO)
  allocate(v(n,n), source = ZERO)
  allocate(x(n,n), source = ZERO)
  allocate(y(n,n), source = ZERO)

  !.............................................................................
  ! Initialization

  call MPI_Init()
  
  do j = 1, n
  do i = 1, n
    x(i,j) = (i - 1) * h
    y(i,j) = (j - 1) * h
  end do
  end do

  u = sin(2 * PI * x) * cos(2 * PI * x)

  print *

  !.............................................................................
  ! Baseline implementation

  start_time = MPI_Wtime()
  
  do i_repetition = 1, n_repetition
    call StencilBaseline(n, h, u, v)
  end do

  stop_time            = MPI_Wtime()

  time_taken           = stop_time - start_time
  time_per_application = time_taken / n_repetition
  
  print *
  print '(A)',          'Baseline implementation'
  print '(A,1X,E15.7)', 'Time taken:',            time_per_application
  print '(A,1X,E15.7)', 'MLUPS:     ', (n-1)**2 / time_per_application

  !.............................................................................
  ! Direct implementation

  start_time = MPI_Wtime()
  
  do i_repetition = 1, n_repetition
    call StencilDirect(n, h, u, v)
  end do

  stop_time            = MPI_Wtime()
  
  time_taken           = stop_time - start_time
  time_per_application = time_taken / n_repetition
  
  print *
  print '(A)',          'Direct implementation'
  print '(A,1X,E15.7)', 'Time taken:',            time_per_application
  print '(A,1X,E15.7)', 'MLUPS:     ', (n-1)**2 / time_per_application

  !.............................................................................
  ! Blocked variant

  start_time = MPI_Wtime()
  
  do i_repetition = 1, n_repetition
    call StencilBlocked(n, h, u, v)
  end do

  stop_time            = MPI_Wtime()

  time_taken           = stop_time - start_time
  time_per_application = time_taken / n_repetition
  
  print *
  print '(A)',          'Blocked variant'
  print '(A,1X,E15.7)', 'Time taken:',            time_per_application
  print '(A,1X,E15.7)', 'MLUPS:     ', (n-1)**2 / time_per_application

  !.............................................................................
  ! Cleanup
  
  call MPI_Finalize()

!===============================================================================

end program Stencil_Test
