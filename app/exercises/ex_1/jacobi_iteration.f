!> \file      jacobi_iteration.f
!> \brief     Program solving the Poisson equation
!> \author    Immo Huismann
!> \date      2014/11/27
!> \copyright Institute of Fluid Mechanics
!>            TU Dresden, 01062 Dresden, Germany
!>
!> \details
!>
!>
!===============================================================================

program Jacobi_Iteration
  use Kind_Parameters, only: RNP
  use Constants,       only: PI
  implicit none

  character(len=*), parameter :: OUTPUT_FILE = 'poisson_results'
  real(RNP),        parameter :: L = 2 * PI ! domain length

  real(RNP), allocatable      :: x_h(:)     ! node locations
  real(RNP), allocatable      :: u_h(:)     ! solution variable
  real(RNP), allocatable      :: u_e(:)     ! exact solution on node locations
  real(RNP), allocatable      :: f_h(:)     ! right-hand side of the equation
  real(RNP), allocatable      :: A(:,:)     ! system matrix

  integer   :: n   ! number of intervals
  integer   :: i   ! loop indices
  integer   :: uni ! index for the output unit
  real(RNP) :: h   ! mesh width

  ! set number of intervals.
  ! it would normally be read from an input file
  n = 20

  ! compute the mesh width
  h = 1

  !.............................................................................
  ! allocation of data
  allocate(x_h(0:n), u_h(0:n), f_h(0:n), u_e(0:n))
  allocate(A(n-1,n-1))

  !.............................................................................
  ! Initialization -- loop notation

  do i = 0, n
    ! compute grid
    x_h(i) = i * h
    ! exact solution
    u_e(i) = sin(x_h(i))

    ! right-hand side of the equation
    f_h(i) = - sin(x_h(i))

    ! intial value for the solution variable
    u_h(i) = 0
  end do

  !.............................................................................
  ! Initialization -- array notation

  x_h = [(i * h, i = 0, n)] ! implicit do loop
  u_e =   sin(x_h)          ! arrays in elemental functions generate arrays
  f_h = - sin(x_h)          ! arrays in elemental functions generate arrays
  u_h = 0                   ! you can set a whole array to a value

  !.............................................................................
  ! Set matrix A
  A = 0

  ! first line
  A(1,1) =  2 / h**2
  A(1,2) = -1 / h**2
  do i = 2, n - 2
    A(i,i-1) = -1 / h**2
    A(i,i  ) =  2 / h**2
    A(i,i+1) = -1 / h**2
  end do

  ! last line
  A(n-1,n-2) = -1 / h**2
  A(n-1,n-1) =  2 / h**2

  !.............................................................................
  ! Solve the system

  call JacobiIteration(A,u_h,f_h)

  !.............................................................................
  ! Output & quit

  print *, 'The maximum error is:', maxval(abs(u_h - u_e))

  ! Write the results into an output file
  open(newunit = uni, file = OUTPUT_FILE)

  do i = 0, n
    ! Fortran requires the index of the output unit (STDOUT, STDERR, some file)
    write(uni,*) x_h(i), u_h(i), u_e(i)
  end do

  close(uni)

contains

!-------------------------------------------------------------------------------
!> \brief   Jacobi iteration for homogeneous Dirichlet BC
!> \author  Immo Huismann
!>
!> \details
!> u_h is assumed to be zero at the locations 0 and n .

subroutine JacobiIteration(A,u_h,f_h)
  real(RNP), intent(in)    :: A(1:,1:) !< system matrix
  real(RNP), intent(inout) :: u_h(0:)  !< solution variable
  real(RNP), intent(in)    :: f_h(0:)  !< right-hand side

  real(RNP) :: u_new(0:ubound(u_h,1))  ! temporary storage for u

  integer, parameter :: ITER_MAX = 10000 ! maximum number of iterations

  integer   :: n, i, j
  integer   :: i_iter
  real(RNP) :: residual

  n = ubound(u_h,1)

  do i_iter = 1, ITER_MAX

    do i = 1, n - 1

      ! add the diagonal part here
      ! the loop over the whole matrix need not be disturbed
      u_new(i) = - f_h(i) + A(i,i) * u_h(i)

      do j = 1, n - 1
        u_new(i) = u_new(i) - A(i,j) * u_h(j)
      end do

      u_new(i) = u_new(i) / A(i,i)
    end do

    ! boundary conditions
    u_new(0) = 0
    u_new(n) = 0

    ! update
    u_h = u_new
  end do

  ! mind that we are solving with -f_h as RHS, so the equation system is
  ! A u = - f, and the residual vector A u + f
  residual = maxval(abs(f_h(1:n-1) + matmul(A,u_h(1:n-1))))

  print *, 'The algebraic residual is:', residual

end subroutine JacobiIteration

!===============================================================================

end program Jacobi_Iteration
