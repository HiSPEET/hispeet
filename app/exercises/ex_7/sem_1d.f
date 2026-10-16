!> summary:   A 1D solver for the convection-diffusion equation
!> author:    Immo Huismann
!> date:      2016/12/07
!> license:   Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

program SEM_1D
  use Kind_Parameters, only: RDP
  use Constants,       only: ZERO, HALF, PI
  use Gauss_Jacobi
  use MPI_F08
  implicit none

  !.............................................................................
  ! field variables
  real(RDP), allocatable :: x     (:,:) !< node locations
  real(RDP), allocatable :: u_h   (:,:) !< numerical solution
  real(RDP), allocatable :: u_ex  (:,:) !< exact solution
  real(RDP), allocatable :: rhs   (:,:) !< right-hand side for time-stepping
  real(RDP), allocatable :: gm_inv(:,:) !< inverse of the global mass matrix

  !.............................................................................
  ! GLL variables
  real(RDP), allocatable :: xi(:)       !< GLL nodes in standard element
  real(RDP), allocatable :: D (:,:)     !< differentiation matrix
  real(RDP), allocatable :: L (:,:)     !< (weak) Laplace matrix
  real(RDP), allocatable :: w (:)       !< integration weights

  !.............................................................................
  ! time-stepping variables
  integer   :: i_time, n_time           !< time step, number of time steps
  real(RDP) :: time_end                 !< time at end of simulation
  real(RDP) :: sim_time                 !< time in simulation
  real(RDP) :: dt                       !< time step width

  real(RDP) :: start_time, stop_time    !< time stamps for measurements
  real(RDP), parameter :: CFL = 0.2_RDP !< Courant-Friedrich-Levy number
  real(RDP), parameter :: DCN = HALF**2 !< Diffusion Courant number

  integer   :: n_element, po            !< number of elements, polynomial degree
  real(RDP) :: lambda                   !< diffusivity
  real(RDP) :: x_min, x_max             !< min. location in mesh, max location
  real(RDP) :: h, d_min                 !< element width, min. width in mesh
  integer   :: i,j,m,e                  !< loop indices
  integer   :: uni                      !< IO unit

  character(len=*), parameter :: INPUT_FILE      = 'input_sem_1d'
  character(len=*), parameter :: OUTPUT_FILES(2) = ['output_1','output_2']

  !-----------------------------------------------------------------------------
  ! Begin of program

  !.............................................................................
  ! Read input

  open(newunit = uni, file = INPUT_FILE)
  read(uni,*) po
  read(uni,*) n_element
  read(uni,*) x_min
  read(uni,*) x_max
  read(uni,*) time_end
  read(uni,*) lambda
  close(uni)

  !.............................................................................
  ! Initialize GLL Variables

  allocate(xi(0:po),D(0:po,0:po),L(0:po,0:po),w(0:po))

  xi = LobattoPoints(po)
  w  = LobattoWeights(xi)
  D  = LobattoDiffMatrix(xi)

  do i = 0, po
  do j = 0, po

    L(i,j) = 0

    do m = 0, po
      L(i,j) = L(i,j) + D(m,i) * w(m) * D(m,j)
    end do

  end do
  end do

  !.............................................................................
  ! Initialize field variables

  allocate(x     (0:po,n_element))
  allocate(u_h   (0:po,n_element))
  allocate(u_ex  (0:po,n_element))
  allocate(rhs   (0:po,n_element))
  allocate(gm_inv(0:po,n_element))

  ! Compute grid
  h = (x_max - x_min) / n_element

  do e = 1, n_element
  do i = 0, po
    x(i,e) = (xi(i) +1 ) * h / 2 + (e - 1) * h + x_min
  end do
  end do

  ! compute the global mass matrix
  do e = 1, n_element
    gm_inv(:,e) = w * h / 2
  end do

  do e = 2, n_element
    gm_inv( 0,e  ) = gm_inv(0,e) + gm_inv(po,e-1)
    gm_inv(po,e-1) = gm_inv(0,e)
  end do

  gm_inv = 1 / gm_inv

  ! Set initial condition
  call SetSolution(x,u_ex,t = ZERO, lambda = lambda)
  u_h  = u_ex

  !.............................................................................
  ! Compute time step width
  d_min    = x(1,1) - x(0,1)

  dt       = min(d_min * CFL, d_min**2 * DCN / lambda)
  n_time   = ceiling(time_end / dt)
  dt       = time_end / n_time
  sim_time = 0

  !-----------------------------------------------------------------------------
  ! Computation

  call GenerateOutput(x,u_h,u_ex,OUTPUT_FILES(1))

  start_time = MPI_Wtime()

  time_marching_loop: do i_time = 1, n_time

    ! compute the right-hand side
    do e = 1, n_element
      rhs(:,e) = - w * matmul(D,u_h(:,e)) - lambda * matmul(L,u_h(:,e)) * (2/h)
    end do

    ! Assembly over element boundaries
    do e = 2, n_element
      rhs( 0,e  ) = (rhs(0,e) + rhs(po,e-1))
      rhs(po,e-1) = rhs(0,e)
    end do

    ! apply Dirichlet boundary conditions
    rhs( 0,        1) = 0
    rhs(po,n_element) = 0

    ! update solution variable
    u_h = u_h + dt * rhs * gm_inv

    sim_time = sim_time + dt

  end do time_marching_loop

  stop_time = MPI_Wtime()

  !.............................................................................
  ! Output & quit
  
  call SetSolution(x,u_ex,t = time_end, lambda = lambda)
  
  call GenerateOutput(x,u_h,u_ex,OUTPUT_FILES(2))

  print *, 'Time taken [s]:', stop_time - start_time
  print *, 'Max. error [-]:', maxval(abs(u_h - u_ex))
  print *, 'L2   error [-]:', sum((u_h - u_ex) * gm_inv * (u_h - u_ex))

contains

!-------------------------------------------------------------------------------
!> Generates an output

subroutine GenerateOutput(x,u_h,u_ex,output_file_name)
  real(RDP),        intent(in) :: x   (:,:)
  real(RDP),        intent(in) :: u_h (:,:)
  real(RDP),        intent(in) :: u_ex(:,:)
  character(len=*), intent(in) :: output_file_name

  integer :: i, e, uni

  open(newunit = uni, file = output_file_name)
  do e = 1, size(x,2)
  do i = 1, size(x,1)
    write(uni,'(E15.7,", ",E15.7,", ",E15.7)') x(i,e), u_h(i,e), u_ex(i,e)
  end do
  end do
  close(uni)

end subroutine GenerateOutput

!-------------------------------------------------------------------------------
!> Sets the solution for the test case

subroutine SetSolution(x,u,t,lambda)
  real(RDP), intent(in)  :: x(:,:)
  real(RDP), intent(out) :: u(:,:)
  real(RDP), intent(in)  :: t
  real(RDP), intent(in)  :: lambda

  real(RDP) :: factor, exp_factor

  factor     = 1 / sqrt(1 + 8 * lambda * t)
  exp_factor = 1 /     (1 + 8 * lambda * t)

  u = factor * exp(-2 * (x - t)**2 * exp_factor)

end subroutine SetSolution

!===============================================================================

end program SEM_1D
