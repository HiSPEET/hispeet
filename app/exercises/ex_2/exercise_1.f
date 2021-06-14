!> \file       exercise_1.f
!> \brief      Code for polynomial interpolation, integration and derivation
!> \author     Immo Huismann
!> \date       2013/12/02
!> \copyright  Institute of Fluid Mechanics &
!>             TU Dresden, 01062 Dresden, Germany
!>
!> \details
!> This source code is designed to interpolate a known function, compute its
!> numerical integral and the numerical derivative.
!>
!> But as this is an exercise there are still some holes to fill :-)
!>
!> The goal is as follows:
!> i)  The parameters in the file 'input' poly sets the used polynomial family.
!>     If poly == 1 then Lobatto-based Lagrange polynomials should be used.
!>     If poly == 2 then Gauss-based Lagrange  polynomials should be used.
!>     If poly == 3 then Lagrange polynomials based on equidistant points are to
!>     be used.
!>
!> ii) The parameter func should set the function to interpolate where one, two
!>     and three correspond to the functions FuncOne and so on.
!>
!> This code produces a printed output that states the analytical and the
!> numerical value of the integral. A graphical output can be created by
!> invoking the 'octave -q scripts/plot_exercise_1.m'.
!>
!> \todo
!> Currently the code creates an output which can be plotted using the script
!> plot_exercise_1.m in the 'scripts' folder.
!> The functionalities of the code are to be expanded as follows:
!>
!> 1)   Set the correct interpolation nodes (start with Lobatto),  and the
!>      differentiation matrix, as well as the differentiation matrix.
!>      The new nodes can be seen on the output created by the plot script, so
!>      an easy verification is present.
!>
!> 2)   Set the coefficients of the interpolant (start with only FuncOne).
!>
!> 3)   Interpolate from the interpolation nodes used for the numerics onto the
!>      plotting nodes (line 146). As the plotting programs only use a linear
!>      interpolation this step results in a smoother output. The result can be
!>      verified together with step 2 in the output.
!>
!> 4)   The numerical integration is still not working, set the analytical value
!>      first and then implement the numerical integration. The values can be
!>      compared by invoking the program (not the plot script).
!>
!> 5)   The numerical differentiation needs an implementation, the result can be
!>      seen on the output created by the plot script. Please note that the
!>      interpolation from the interpolation nodes to the plotting nodes has to
!>      work to see a change in the output.
!>
!> 6)   After the program is enabled to interpolate, differentiate and integrate
!>      the input parameters should be taken more seriously. In the details the
!>      definitions are laid out. Implement the dependence from func in 2) and
!>      5), and afterwards the dependence from poly in 1).
!>      Mind that integration weights need to be set for the equidistant
!>      Lagrange integration, these can either be derived by using those from
!>      the spline interpolation or from a coordinate transformation to the
!>      Gauss points.
!>
!> While implementing these six parts mind the comments, most problems are
!> explained in them. The formulas for interpolation, derivation and integration
!> are present in the set of slides on the CFDII homepage.
!>
!===============================================================================

program Exercise_1
  use Kind_Parameters, only: RNP
  use Constants,       only: ZERO, ONE, TWO, PI
  use Gauss_Jacobi
  use Lagrange_Interpolation
  use Execution_control
  implicit none

  integer :: po   ! Number of nodes inside the element
  integer :: n    ! Number of plotting nodes
  integer :: poly ! the kind of polynomial to choose
  integer :: func ! the function used for testing

! standard element variables
  real(RNP), allocatable :: xi  (:)   ! the node locations inside the element
  real(RNP), allocatable :: wi  (:)   ! the integration weights
  real(RNP), allocatable :: diff(:,:) ! the differentiation matrix


  real(RNP), allocatable :: u_h (:) ! coeff. of the interpolation polynomial
  real(RNP), allocatable :: du_h(:) ! coeff. of its derivative

  real(RNP), allocatable :: u_p(:), du_p(:) ! same but on the plot nodes
  real(RNP), allocatable :: x_p(:)          ! nodes to plot on
  real(RNP), allocatable :: u_e(:), du_e(:) ! exact values on plot nodes

  real(RNP) :: integ     ! integral of interpolant
  real(RNP) :: integ_ana ! analytical integral

  integer :: i, j ! loop indices
  integer :: uni  ! output unit

  NAMELIST /input_e1/ po,n,poly,func

  ! read po, n, poly and func from the input file
  open(NEWUNIT=uni,FILE='input')
  read(uni,input_e1)
  close(uni)

  ! allocate the variables for the standard element
  allocate(xi(0:po), wi(0:po), diff(0:po,0:po))

  ! allocate the coefficients of the interpolant and its derivative
  allocate(u_h(0:po), du_h(0:po))

  ! allocate the variables for plotting purposes
  allocate(x_p(0:n), u_p(0:n), du_p(0:n), u_e(0:n), du_e(0:n))

  ! compute the nodes for plotting
  do i=0, n
    x_p(i)  = -ONE + (TWO * i) / n
  end do

  ! set standard element variables .............................................
  ! The node locations, the integration weights and the differentation matrix
  ! are needed on the standard element.
  ! For the Gauss- and Lobatto-Lagrange-polynomials functions are present in
  ! gauss_jacobi.f
  ! For the general Lagrange polynomials they are in lagrange_polynomials.f

  ! TODO 1 - Set the points used for the interpolation
  do i=0, po
    xi(i) = - 1 + (TWO * i) / po
  end do
  wi   = 0
  diff = 0

  ! output of basis functions ..................................................
  open(newunit=uni,file='output/e1_nodal')
  write(uni,*) x_p
  do i=0, po
    write(uni,*) [(Lag_Polynomial(i,xi,x_p(j)), j = 0, n)]
!!!!write(uni,*) [(LagrangePolynomial(i, xi, x_p(j)), j = 0, n)]
  end do
  close(uni)

  !.............................................................................
  ! Interpolation of the funktion u(x) = func(x)

  ! Calculation of the coefficients
  ! TODO 2 - Set the coefficients of the interpolant
  u_h = ZERO

  !.............................................................................
  ! Numerical integration of func on the interval [-1,1], comparison
  ! compare it to the analytical value

  ! TODO 4 - calculate the numerical integral, set the analytical one
  integ_ana = ZERO
  integ     = ZERO

  print '(2X,A,2X,E15.7)', 'The analytical integral is:', integ_ana
  print '(2X,A,2X,E15.7)', 'The numerical integral is: ', integ
  print '(2X,A,2X,E15.7)', 'Error of integral:         ', abs(integ - integ_ana)

  ! Calculation of the derivative ..............................................


  ! TODO 5 - calculate the coefficients of the derivative
  du_h = ZERO

  ! output .....................................................................

  ! TODO 3 - interpolate onto the plotting grid
  u_p  = ZERO
  du_p = ZERO

  ! write x, du, dux
  select case(func)
  case(1)
     u_e =     FuncOne  (x_p)
    du_e = DiffFuncOne  (x_p)
  case(2)
     u_e =     FuncTwo  (x_p)
    du_e = DiffFuncTwo  (x_p)
  case(3)
     u_e =     FuncThree(x_p)
    du_e = DiffFuncThree(x_p)
  end select

  open(newunit=uni,file='output/e1_data')
  write(uni,*) po, poly, n, 0, 0
  do i = 1, n
    write(uni,*)  x_p(i), u_p(i), du_p(i), u_e(i), du_e(i)
  end do
  close(uni)

contains

!===============================================================================
! test functions

!-------------------------------------------------------------------------------
!> \brief   First test function
!> \author  Immo Huismann

elemental function FuncOne(x) result(y)
  real(RNP), intent(in) :: x
  real(RNP)             :: y

  y = real(sin(PI * x),RNP)

end function FuncOne

!-------------------------------------------------------------------------------
!> \brief   Second test function
!> \author  Immo Huismann

elemental function FuncTwo(x) result(y)
  real(RNP), intent(in) :: x
  real(RNP)             :: y

  y = 1 / (1 + 25 * x ** 2)

end function FuncTwo

!-------------------------------------------------------------------------------
!> \brief   Third test function
!> \author  Immo Huismann

elemental function FuncThree(x) result(y)
  real(RNP), intent(in) :: x
  real(RNP)             :: y

  y = abs(x)

end function FuncThree

!===============================================================================
! their integrals

!-------------------------------------------------------------------------------
!> \brief   Integral of the first test function
!> \author  Immo Huismann

elemental function IFuncOne() result(y)
  real(RNP) :: y

  y = 0

end function IFuncOne

!-------------------------------------------------------------------------------
!> \brief   Integral of the second test function
!> \author  Immo Huismann

elemental function IFuncTwo() result(y)
  real(RNP) :: y

  y = TWO / 5  * atan(5 * ONE)

end function IFuncTwo

!-------------------------------------------------------------------------------
!> \brief   Integral of the third test function
!> \author  Immo Huismann

elemental function IFuncThree() result(y)
  real(RNP) :: y

  y = 1

end function IFuncThree

!===============================================================================
! and their derivatives

!-------------------------------------------------------------------------------
!> \brief   Derivative of the first test function
!> \author  Immo Huismann

elemental function DiffFuncOne(x) result(y)
  real(RNP), intent(in) :: x
  real(RNP)             :: y

  y = real(PI * cos(PI * x),RNP)

end function DiffFuncOne

!-------------------------------------------------------------------------------
!> \brief   Derivative of the second test function
!> \author  Immo Huismann

elemental function DiffFuncTwo(x) result(y)
  real(RNP), intent(in) :: x
  real(RNP)             :: y

  y = - 1 / ((1 + 25 * x ** 2) **2) * (50 * x)

end function DiffFuncTwo

!-------------------------------------------------------------------------------
!> \brief   Derivative of the third test function
!> \author  Immo Huismann

elemental function DiffFuncThree(x) result(y)
  real(RNP), intent(in) :: x
  real(RNP)             :: y

  y = sign(ONE,x)

end function DiffFuncThree

!===============================================================================

end program Exercise_1
