!> \file       exercise_2.f 
!> \brief      Solver for the one-dimensional Helmholtz equation
!> \author     Immo Huismann
!> \date       2013/12/02
!> \copyright  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!> \details
!> This code solves the Helmholtz equation, at least it would if it was
!> implemented correctly.
!> To get the code up and running on only one element the TODOs 1-4 have to be
!> solved.
!> To get it to run with more than one element the TODOs 2,4,5 and 6 have to be
!> implemented correctly for multiple elements.
!===============================================================================

program Exercise_2
  use Kind_Parameters,      only: RNP
  use Constants,            only: HALF, ONE, TWO, PI
  use Gauss_Jacobi,         only: GLL_Points, GLL_Weights, GLL_DiffMatrix
  use Lagrange_Interpolation, only: LagrangePolynomial
  use Execution_Control,    only: Error
  implicit none

  ! Discretization parameters
  integer   :: po     = 7   !< polynomial degree
  integer   :: n_ele  = 1   !< number of elements
  integer   :: n_plot = 200 !< number of plotting points per element
  real(RNP) :: lambda = 1   !< Helmholtz parameter

  ! Parameters for the solution process and 
  real(RNP) :: r_max    = 1.0E-10
  integer   :: iter_max = 2000
  
  ! GLL variables
  real(RNP), allocatable :: xi(:)     !< GLL standard element nodes
  real(RNP), allocatable :: wi(:)     !< GLL standard element weights
  real(RNP), allocatable :: diff(:,:) !< GLL standard element differentiation mat

  real(RNP), allocatable :: x_h (:,:) !< the node locations (elementwise)
  real(RNP), allocatable :: u_h (:,:) !< the solution at the nodes
  real(RNP), allocatable :: f_h (:,:) !< the discrete RHS at the nodes
  real(RNP), allocatable :: u_ex(:,:) !< exact solution

  ! standard element matrices, numerical system variables
  real(RNP), allocatable :: mass_mat_diag(:)  !< mass matrix, main diagonal
  real(RNP), allocatable :: stiffness_mat(:,:)!< stiffness matrix
  real(RNP) :: h ! element width

  ! plotting points and plotting variable
  real(RNP), allocatable :: x_plot  (:,:) !< plotting nodes
  real(RNP), allocatable :: u_h_plot(:,:) !< solution on plotting nodes
  real(RNP), allocatable :: u_plot  (:,:) !< exact solution on plotting nodes

  real(RNP), allocatable :: xi_plot(:)  !< standard element node locations

  ! problem dependent parameters
  integer   :: func     = 1 !< desired test function
  procedure(TestFuncOne), pointer :: test_rhs => TestFuncOne
  procedure(TestSolOne) , pointer :: test_sol => TestSolOne
  real(RNP) :: bc_n !< Neumann BC

  integer :: uni     !< io-unit
  integer :: i, j, k !< loop indices

  ! the error in different norms
  real(RNP) :: err_max !< maximum norm
  real(RNP) :: err_L2  !< L2 norm

  NAMELIST /input_e2/ po, n_ele, lambda, n_plot, r_max, func

  !.............................................................................
  ! read input, stats the test case
  
  open(newunit = uni, file = 'input.prm')
  read(uni,input_e2)
  close(uni)

  print '(A)', 'SEM+CG for: lambda * u(x) - u"(x) = f(x), u(0)=g, u''(1) = h'
  print '(2X,A,I4)', 'Using a polynomial degree of: ', po
  print '(2X,A,I4)', 'The number of elements is:    ', n_ele

  !.............................................................................
  ! select the function for the right-hand side, the solution and the Neumann BC
  
  select case(func)
  case(1)
    test_rhs => TestFuncOne
    test_sol => TestSolOne
    bc_n = 0
  case(2)
    test_rhs => TestFuncTwo
    test_sol => TestSolTwo
    bc_n = real(- PI, RNP)
  case(3)
    test_rhs => TestFuncThree
    test_sol => TestSolThree
    bc_n = 1
  case default
    call Error('Exercise 2', 'Unknown function desired')  
  end select

  ! array allocation ...........................................................

  ! the GLL variables and the system matrices
  allocate(xi(0:po), wi(0:po), diff(0:po,0:po))
  allocate(mass_mat_diag(0:po), stiffness_mat(0:po,0:po))

  ! the node locations and the coefficients for the solution and the RHS
  allocate(x_h(0:po,n_ele), u_h(0:po,n_ele), f_h(0:po,n_ele), u_ex(0:po,n_ele))

  ! the plotting points and the solution to plot
  allocate(x_plot(0:n_plot,n_ele), u_plot(0:n_plot,n_ele), &
       u_h_plot(0:n_plot,n_ele))

  !-----------------------------------------------------------------------------
  ! system setup

  ! set the element width
  h = TWO / n_ele

  ! get the standard element variables, node locations, system matrices
  call GetGLLVariables  (po,xi,wi,diff)
  call GetNodeLocations (h, xi, x_h)
  call GetSystemMatrices(wi, diff, mass_mat_diag, stiffness_mat)

  ! get the right-hand side for the given problem
  f_h = test_rhs(lambda, x_h)

  !-----------------------------------------------------------------------------
  ! solution process

  ! integrate the right-hand side for the weighted problem
  call IntegrateRHS(h,mass_mat_diag, f_h)

  ! set u_h to zero as a starting guess
  u_h = 0
  
  ! set the Dirichlet BC at x = -1
  u_h(0:0,1:1) = test_sol(x_h(0:0,1:1))
  
  ! solve the system utilizing a Conjugated Gradients method
  call ConjugatedGradients (h, lambda, bc_n, &
       mass_mat_diag, stiffness_mat, r_max, iter_max, f_h, u_h)

  print '(2X,A,2X,E15.7)', 'Achieved residual:        ', r_max
  print '(2X,A,2X,I5.5)' , 'Used number of iterations:', iter_max
  
  ! calculate the error
  u_ex = test_sol(x_h)

  ! calculate the maximum error
  err_max = maxval(abs(u_ex - u_h))

  ! calculate the L2 error
  err_L2 = 0
  do k = 1, n_ele
    err_L2 = err_L2 + h / 2 * sum(wi * (u_ex(:,k) - u_h(:,k)) ** 2)
  end do

  print '(2X,A,2X,E15.7)', 'The maximum error is: ', err_max
  print '(2X,A,2X,E15.7)', 'The L2 error is:      ', err_L2
  
  !-----------------------------------------------------------------------------
  ! output

  ! interpolation from the solution nodes onto the plot nodes
  allocate(xi_plot(0:n_plot))
  xi_plot = [(TWO / n_plot * i, i = 0, n_plot)] - 1
  
  ! distribute the plotting nodes on the domain
  do k = 1, n_ele
    x_plot(:,k) = -1 +  h * (xi_plot + 1) / 2 + (k-1) *  h
  end do
  
  ! interpolate elementwise from the solution nodes onto the plotting nodes
  do k = 1, n_ele
  do j = 0, n_plot
    u_h_plot(j,k) = sum([(u_h(i,k) * LagrangePolynomial(i,xi,xi_plot(j)), i = 0, po)])
  end do
  end do

  u_plot = test_sol(x_plot)

  open(newunit = uni, file = 'output_e2_data')
  write(uni, '(3(I0, 1X))') po, n_ele, n_plot
  do k = 1, n_ele
  do j = 0, n_plot
    write(uni, '(3(E15.7, 1X))') x_plot(j,k), u_h_plot(j,k), u_plot(j,k)
  end do
  end do
  close(uni)
  
  open(newunit = uni, file = 'output_e2_nodal')
  do k = 1, n_ele
  do i = 0, po
    write(uni,'(3(E15.7, 1X))') x_h(i,k), u_h(i,k), test_sol(x_h(i:i,k:k))
  end do
  end do
  close(uni)

  deallocate(xi, wi, diff, xi_plot)
  deallocate(x_h,u_h,f_h,x_plot,u_plot)

contains

!===============================================================================
! Setup of the numerical system

!-------------------------------------------------------------------------------
!> \brief   This routine sets the GLL Variables
!> \details This routine gets the node locations, the integration weights, and 
!>          the differentiation matrix depending on the polynomial degree p.

pure subroutine GetGLLVariables(po, xi, wi, diff)
  integer  , intent(in)  :: po          !< desired polynomial degree
  real(RNP), intent(out) :: xi  (0:)    !< node locations
  real(RNP), intent(out) :: wi  (0:)    !< weights
  real(RNP), intent(out) :: diff(0:,0:) !< differentiation matrix

  integer :: i

  ! TODO 1 - Get the GLL Variables

  ! currently equidistant nodes are set
  xi   = [(-1 + TWO / po * i, i = 0, po)]
  wi   = 0  
  diff = 0

end subroutine GetGLLVariables

!-------------------------------------------------------------------------------
!> Sets the node locations on the interval [-1, 1]

pure subroutine GetNodeLocations(h, xi, x)
  real(RNP), intent(in)  :: h         !< element width
  real(RNP), intent(in)  :: xi(0:)    !< GLL nodes on standard element ([-1, 1])
  real(RNP), intent(out) :: x (0:,1:) !< nodes locations to set

  integer :: po, n_ele
  integer :: i, k

  po    = ubound(x ,1)
  n_ele = ubound(x ,2)

  ! TODO 2 - Set the node locations for the mesh with the mesh width h.
  ! The mesh starting point is -1, the last node is located on 1.

  ! Currently an equidistant mesh is used, but GLL points inside each element
  ! are desired
  do k = 1, n_ele
    x(:,k) = h * (k - 1) + (xi + 1) / 2 * h - 1
  end do

end subroutine GetNodeLocations

!-------------------------------------------------------------------------------
!> Sets the standard element matrices (mass and stiffness)

pure subroutine GetSystemMatrices(wi, diff, mass_mat_diag, stiffness_mat)
  real(RNP), intent(in)  :: wi  (0:)             !< integration weights
  real(RNP), intent(in)  :: diff(0:,0:)          !< differentiation matrix
  real(RNP), intent(out) :: mass_mat_diag(0:)    !< mass matrix, main diagonal
  real(RNP), intent(out) :: stiffness_mat(0:,0:) !< stiffness matrix

  real(RNP) :: mass_mat(0:ubound(wi,1),0:ubound(wi,1))

  integer :: po  ! polynomial degree
  integer :: i   ! loop index

  po = ubound(wi,1)

  ! TODO 3 - Set the mass matrix and the stiffness matrix correctly
  ! Currently spline values are utilized
  mass_mat_diag     = TWO / po
  mass_mat_diag(0)  = ONE / po
  mass_mat_diag(po) = ONE / po
  mass_mat          = 0
  
  stiffness_mat = 0

end subroutine GetSystemMatrices

!-------------------------------------------------------------------------------
!> A weighted integration using a test function phi_i (application of mass mat)

pure subroutine IntegrateRHS(h, mass_mat_diag, f_h)
  real(RNP), intent(in)    :: h
  real(RNP), intent(in)    :: mass_mat_diag(0:)
  real(RNP), intent(inout) :: f_h(0:,1:)

  integer :: po, n_ele ! polynomial degree, number of elements
  integer :: k         ! loop index

  po    = ubound(f_h, 1)
  n_ele = ubound(f_h, 2)

  ! TODO 5 - Add the metric coefficient for more than one element
  do k = 1, n_ele
    f_h(:,k) = mass_mat_diag * f_h(:,k)
  end do

end subroutine IntegrateRHS

!===============================================================================
! Solver with the Helmholtz Operator and assembly process

!-------------------------------------------------------------------------------
!> A Conjugated Gradients algorithm for the helmholtz equation

pure subroutine ConjugatedGradients &
     (h,lambda,bc_n,mass_mat_diag,stiffness_mat,r_max,iter_max,f,u)
  real(RNP), intent(in)    :: h
  real(RNP), intent(in)    :: lambda
  real(RNP), intent(in)    :: bc_n    ! Neumann boundary condition at x = 1
  real(RNP), intent(in)    :: mass_mat_diag(0:)
  real(RNP), intent(in)    :: stiffness_mat(0:,0:)
  real(RNP), intent(inout) :: r_max
  integer  , intent(inout) :: iter_max
  real(RNP), intent(in)    :: f(0:,1:)
  real(RNP), intent(inout) :: u(0:,1:)

  integer   :: iter ! iteration count
  real(RNP) :: alpha, beta
  real(RNP) :: p     (0:ubound(f,1),1:ubound(f,2)) ! search direction
  real(RNP) :: q     (0:ubound(f,1),1:ubound(f,2)) ! impact on residual by p
  real(RNP) :: r     (0:ubound(f,1),1:ubound(f,2)) ! next residual
  real(RNP) :: r_0   (0:ubound(f,1),1:ubound(f,2)) ! old residual
  real(RNP) :: ml_inv(0:ubound(f,1),1:ubound(f,2)) ! mass lump inverted

  integer   :: po, n_ele

  po    = ubound(f,1)
  n_ele = ubound(f,2)

  ! The inverted mass lump eliminates double entries due to the local
  ! formulation in the scalar product
  ml_inv = 1
  call AssembleResidual(ml_inv)
  ml_inv = 1 / ml_inv

  ! Get the residual of the Helmholtz equation
  r_0 = HelmholtzResidual(h,lambda,mass_mat_diag,stiffness_mat,u,f)
  call AssembleResidual(r_0)
  
  ! A Dirichlet BC is imposed at x = -1, the residual is set to 0
  r_0(0,1) = 0
  ! A Neumann BC is imposed at x = 1, it is added to the residual
  r_0(po, n_ele) = r_0(po, n_ele) + bc_n

  ! set the first basis vector (first search direction) as the residual
  p   = r_0

  ! exit if the residual is too low to achieve a better solution
  if(sum(r_0 * r_0 * ml_inv) == 0) then
    iter_max = 0
    r_max    = 0
    return
  end if

  do iter = 1, iter_max
    ! get the impact of direction p on the residual, store in q
    q = HelmholtzOperator(h,lambda,mass_mat_diag,stiffness_mat,p)
    ! assemble q to ensure the continuity (global vs. local formulation)
    call AssembleResidual(q)
    
    ! A Dirichlet BC is imposed at x = -1, the residual is set to 0
    q(0,1) = 0

    ! Get the distance to move in direction p
    alpha = sum(r_0 * r_0 * ml_inv) / sum(q * p * ml_inv)

    ! advance distance alpha in direction p (q = H p for r)
    u     = u   + alpha * p
    r     = r_0 - alpha * q

    ! exit if the residual is too small
    if(sqrt(sum(r * r)) < r_max) exit

    ! get the residual reduction by p
    beta = sum(r * r * ml_inv) / sum(r_0 * r_0 * ml_inv)

    ! next search direction
    p = r + beta * p

    ! save the residual for the next step
    r_0 = r

  end do

  ! save the iteration count and the achieved residual
  iter_max = iter - 1
  r_max    = sqrt(sum(r*r))

end subroutine ConjugatedGradients

!-------------------------------------------------------------------------------
!> The residual of the Helmholtz equation res_e = f_e - H_e u_e

pure function HelmholtzResidual(h,lambda,mass_mat_diag,stiffness_mat,u,f) &
     result(res)
  real(RNP), intent(in)  :: h                      ! element width
  real(RNP), intent(in)  :: lambda                 ! Helmholtz parameter
  real(RNP), intent(in)  :: mass_mat_diag(0:)      ! mass matrix
  real(RNP), intent(in)  :: stiffness_mat(0:,0:)   ! stiffness matrix
  real(RNP), intent(in)  :: u(0:,1:)               ! variable to apply to
  real(RNP), intent(in)  :: f(0:,1:)               ! discrete RHS
  real(RNP)              :: res(0:ubound(u,1),1:ubound(u,2)) ! result

  ! The residual is defined as the discrete RHS minus the Helmholtz op of u
  res = f - HelmholtzOperator(h,lambda,mass_mat_diag,stiffness_mat,u)

end function HelmholtzResidual

!-------------------------------------------------------------------------------
!> The Helmholtz Operator res_e = H_e u_e

pure function HelmholtzOperator(h,lambda,mass_mat_diag,stiffness_mat,u) &
     result(res)
  real(RNP), intent(in)  :: h                      !< element width
  real(RNP), intent(in)  :: lambda                 !< Helmholtz parameter
  real(RNP), intent(in)  :: mass_mat_diag(0:)      !< mass matrix
  real(RNP), intent(in)  :: stiffness_mat(0:,0:)   !< stiffness matrix
  real(RNP), intent(in)  :: u(0:,1:)               !< variable to apply to
  real(RNP)              :: res(0:ubound(u,1),1:ubound(u,2)) ! result

  integer :: po, n_ele
  integer :: k

  po    = ubound(u, 1)
  n_ele = ubound(u, 2)

  ! TODO 4 - Implement the Helmholtz operator

  ! apply the Helmholtz operator elementwise to u
  do k = 1, n_ele
    res(:,k) = lambda * u(:,k)
  end do

end function HelmholtzOperator

!-------------------------------------------------------------------------------
!> Assembles a residual

pure subroutine AssembleResidual(res)
  real(RNP), intent(inout) :: res(0:,1:) !< residual to assemble

  integer :: po, n_ele
  integer :: k

  po    = ubound(res, 1)
  n_ele = ubound(res, 2)

  ! TODO 6 - Assemble the residual
  do k = 1, n_ele
    res(:,k) = res(:,k)
  end do

end subroutine AssembleResidual

!===============================================================================
! Test functions for the RHS

!-------------------------------------------------------------------------------
!> A simple first test function

pure function TestFuncOne(lambda, x) result(y)
  real(RNP), intent(in) :: lambda
  real(RNP), intent(in) :: x(0:,1:)
  real(RNP)             :: y(0:ubound(x,1),1:ubound(x,2))

  y = lambda * TestSolOne(x) - 1

end function TestFuncOne

!-------------------------------------------------------------------------------
!> A second simple test function

pure function TestFuncTwo(lambda, x) result(y)
  real(RNP), intent(in) :: lambda
  real(RNP), intent(in) :: x(0:,1:)
  real(RNP)             :: y(0:ubound(x,1),1:ubound(x,2))

  y = real((lambda + PI ** 2) * TestSolTwo(x), RNP)

end function TestFuncTwo

!-------------------------------------------------------------------------------
!> A test function for the exercise
!> \details
!> This test function is to be chosen in such a way that the result is |x| ^ 3.

pure function TestFuncThree(lambda, x) result(y)
  real(RNP), intent(in) :: lambda
  real(RNP), intent(in) :: x(0:,1:)
  real(RNP)             :: y(0:ubound(x,1),1:ubound(x,2))

  ! TODO 7: Correctly determine the right-hand side
  y = lambda * x

end function TestFuncThree

!===============================================================================
! corresponding solutions

!-------------------------------------------------------------------------------
!> A polynomial of degree 2, for po >= 2 the solution should match this exactly

pure function TestSolOne(x) result(y)
  real(RNP), intent(in) :: x(0:,1:)
  real(RNP)             :: y(0:ubound(x,1),1:ubound(x,2))

  y = HALF * x ** 2 - x - HALF * 3 

end function TestSolOne

!-------------------------------------------------------------------------------
!> A sinus, the solution should match this better for higher po

pure function TestSolTwo(x) result(y)
  real(RNP), intent(in) :: x(0:,1:)
  real(RNP)             :: y(0:ubound(x,1),1:ubound(x,2))

  y = real(sin(PI * x), RNP)

end function TestSolTwo

!-------------------------------------------------------------------------------
!> A third function

pure function TestSolThree(x) result(y)
  real(RNP), intent(in) :: x(0:,1:)
  real(RNP)             :: y(0:ubound(x,1),1:ubound(x,2))

  y = abs(x) ** 3

end function TestSolThree

end program Exercise_2
