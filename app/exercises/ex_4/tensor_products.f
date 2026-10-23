!> summary:  Test program for multi-dimensional interpolation
!> author:   Immo Huismann
!> date:     2013/12/02
!> license:  Institute of Fluid Mechanics TU Dresden, 01062 Dresden, Germany
!>
!> This program allows for the multi dimensional interpolation of a given
!> function u and the numerical differentiation and integration of it.
!> The solutions are interpolated onto a plotting grid and the error is
!> estimated on it.
!>
!> The program won't work at first as the TODOs 1-4 have to be solved.
!===============================================================================

program Tensor_Products
  use Kind_Parameters,        only: RNP
  use Constants,              only: ONE, TWO, PI
  use Gauss_Jacobi,           only: LobattoPoints, LobattoWeights, LobattoDiffMatrix
  use Lagrange_Interpolation, only: LagrangePolynomial
  implicit none

  ! problem parameters .........................................................
  
  integer :: po    (2) = [ 17,  13] !< polynomial degree in two dimensions
  integer :: n_plot(2) = [200, 200] !< number of plot points - 1
  integer :: func  (2) = [1, 2]     !< two functions to interpolate

  ! mind that the node locations, the approximation of u and its derivative
  ! are stored in a matrix format, due to the tensor product form
  ! x_h(i,j) stands for the coefficient of the i-th basis in the first direction 
  ! and the j-th basis in the second direction
  
  real(RNP), allocatable :: x_h(:,:,:), u_h(:,:), du_h(:,:,:)

  real(RNP) :: integ, integ_exact

  ! Lobatto variables for both dimensions
  real(RNP), allocatable :: xi_1(:), wi_1(:), diff_1(:,:)
  real(RNP), allocatable :: xi_2(:), wi_2(:), diff_2(:,:)

  ! variables for plotting
  real(RNP), allocatable ::  x_plot  (:,:,:)
  real(RNP), allocatable ::  u_plot  (:,:  )
  real(RNP), allocatable ::  du_plot (:,:,:)
  real(RNP), allocatable ::  err_plot(:,:  )

  real(RNP), allocatable :: S_1(:,:), S_2(:,:) !< interpolation matrices

  integer :: i,j  !< loop indices
  integer :: unit !< io unit
  
  interface
    
   function Funct(x) result(y)
     import
      real(RNP), intent(in) :: x(:,:)
      real(RNP)             :: y(size(x,1),size(x,2))      
    end function Funct
   
   function IntFunct() result(y)
     import
     real(RNP) :: y
   end function IntFunct    
    
  end interface

  procedure(Funct), pointer :: func_1
  procedure(Funct), pointer :: func_2
    
  procedure(IntFunct), pointer :: int_func_1
  procedure(IntFunct), pointer :: int_func_2
  
  NAMELIST /input_tensor_products/ po, n_plot, func

  ! input ......................................................................
  
  ! read the input parameters
  open(newunit = unit, file = 'tensor_products.prm')
  read(unit, input_tensor_products)
  close(unit)

  ! set the test functions accordingly
  select case(func(1))
  case(1)
    func_1     =>    FuncOne
    int_func_1 => IntFuncOne
  case(2)
    func_1     =>    FuncTwo
    int_func_1 => IntFuncTwo
  case(3)
    func_1     =>    FuncThree
    int_func_1 => IntFuncThree
  case(4)
    func_1     =>    FuncFour
    int_func_1 => IntFuncFour
  end select

  select case(func(2))
  case(1)
    func_2     =>    FuncOne
    int_func_2 => IntFuncOne
  case(2)
    func_2     =>    FuncTwo
    int_func_2 => IntFuncTwo
  case(3)
    func_2     =>    FuncThree
    int_func_2 => IntFuncThree
  case(4)
    func_2     =>    FuncFour
    int_func_2 => IntFuncFour
  end select

  ! initialization .............................................................
  
  ! allocation of field variables
  allocate( x_h(0:po(1),0:po(2),1:2))
  allocate( u_h(0:po(1),0:po(2)    ))
  allocate(du_h(0:po(1),0:po(2),1:2))

  ! allocation of Lobatto variables
  allocate(xi_1(0:po(1)), wi_1(0:po(1)), diff_1(0:po(1),0:po(1)))
  allocate(xi_2(0:po(2)), wi_2(0:po(2)), diff_2(0:po(2),0:po(2)))

  ! set the Lobatto variables for the first direction
  xi_1   = LobattoPoints    (po(1))
  wi_1   = LobattoWeights   (xi_1 )
  diff_1 = LobattoDiffMatrix(xi_1 )

  ! and for the second
  xi_2   = LobattoPoints    (po(2))
  wi_2   = LobattoWeights   (xi_2 )
  diff_2 = LobattoDiffMatrix(xi_2 )

  ! node locations .............................................................
  
  ! TODO 1 - Set the node locations

  do j = 0, po(2)
   do i = 0, po(1)
     ! set the nodes in the first  direction
     x_h(i,j,1) = i * (TWO / po(1)) - 1
     ! set the nodes in the second direction
     x_h(i,j,2) = j * (TWO / po(2)) - 1
   end do
  end do

  ! interpolate the functions ..................................................
  ! TODO 2 - Set the coefficients for the interpolation
  u_h                  = 0

  ! differentiation ............................................................
  ! TODO 3 - differentiate in the first and second direction

  ! differentiate u_h in the two directions and store it in du_h
  du_h(:,:,1) = 1
  du_h(:,:,2) = 2

  ! integration ................................................................
  integ_exact = int_func_1() * int_func_2()

  ! TODO 4 - Implement the numerical integration
  integ     = 0

  print '(A,E15.7)', 'The integral should be:  ', integ_exact
  print '(A,E15.7)', 'The calculated valus is: ', integ

  ! ouput ......................................................................

  ! allocate the plot variables and the interpolation matrices
  allocate(  x_plot(0:n_plot(1),0:n_plot(2),1:2))
  allocate(  u_plot(0:n_plot(1),0:n_plot(2)))
  allocate( du_plot(0:n_plot(1),0:n_plot(2),1:2))
  allocate(err_plot(0:n_plot(1),0:n_plot(2)))
  allocate(S_1(0:n_plot(1),0:po(1)), S_2(0:n_plot(2),0:po(2)))

  ! generate the plotting nodes
  do j = 0, n_plot(2)
   do i = 0, n_plot(1)
     ! set the nodes in the first  direction
     x_plot(i,j,1) = i * (TWO / n_plot(1)) - 1
     ! set the nodes in the second direction
    x_plot(i,j,2) = j * (TWO / n_plot(2)) - 1
   end do
  end do

  ! generate the interpolation matrices from x_h to x_plot
  do j = 0, po(1)
   do i = 0, n_plot(1)
     S_1(i,j) = LagrangePolynomial(j,x_h(:,0,1),x_plot(i,0,1))
   end do
  end do
  
  do j = 0, po(2)
   do i = 0, n_plot(2)
     S_2(i,j) = LagrangePolynomial(j,x_h(0,:,2),x_plot(0,i,2))
   end do
  end do

  ! interpolate onto the plotting nodes using a tensor product
  u_plot         = matmul(S_1,matmul( u_h       ,transpose(S_2)))
  du_plot(:,:,1) = matmul(S_1,matmul(du_h(:,:,1),transpose(S_2)))
  du_plot(:,:,2) = matmul(S_1,matmul(du_h(:,:,2),transpose(S_2)))

  err_plot = func_1(x_plot(:,:,1)) * func_2(x_plot(:,:,2)) - u_plot

  ! write the output
  open(newunit = unit, file = 'tensor_products.dat')
  write(unit, *) po(1), po(2), n_plot(1), n_plot(2), 0, 0
  do j = 0, n_plot(2)
   do i = 0, n_plot(1)
     write(unit, *) x_plot(i,j,1:2), u_plot(i,j), du_plot(i,j,1:2), err_plot(i,j)
   end do
  end do
  close(unit)

  ! cleaning ...................................................................
  
  deallocate(x_h,u_h,du_h)
  deallocate(xi_1,wi_1,diff_1)
  deallocate(xi_2,wi_2,diff_2)
  deallocate(x_plot,u_plot,du_plot)
  deallocate(S_1,S_2)

contains

!===============================================================================
! Test functions


!-------------------------------------------------------------------------------
!> Identity function

function FuncOne(x) result(y)
  real(RNP), intent(in) :: x(:,:)
  real(RNP)             :: y(size(x,1),size(x,2))

  y = x

end function FuncOne

!-------------------------------------------------------------------------------
!> Function sin(pi * x)

function FuncTwo(x) result(y)
  real(RNP), intent(in) :: x(:,:)
  real(RNP)             :: y(size(x,1),size(x,2))

  y = real(sin(PI * x), RNP)

end function FuncTwo

!-------------------------------------------------------------------------------
!> Function exp(x)

function FuncThree(x) result(y)
  real(RNP), intent(in) :: x(:,:)
  real(RNP)             :: y(size(x,1),size(x,2))

  y = exp(x)

end function FuncThree

!-------------------------------------------------------------------------------
!> Function sign(x)

function FuncFour(x) result(y)
  real(RNP), intent(in) :: x(:,:)
  real(RNP)             :: y(size(x,1),size(x,2))

  y = sign(ONE,x)

end function FuncFour

!===============================================================================
! Integrals of test functions

!-------------------------------------------------------------------------------
!> Integral of FuncOne over [-1,1]

function IntFuncOne() result(y)
  real(RNP) :: y

  y = 0

end function IntFuncOne

!-------------------------------------------------------------------------------
!> Integral of FuncTwo over [-1,1]

function IntFuncTwo() result(y)
  real(RNP) :: y

  y = 0

end function IntFuncTwo

!-------------------------------------------------------------------------------
!> Integral of FuncThree over [-1,1]

function IntFuncThree() result(y)
  real(RNP) :: y

  y = exp(ONE) - exp(-ONE)

end function IntFuncThree

!-------------------------------------------------------------------------------
!> Integral of FuncFour over [-1,1]

function IntFuncFour() result(y)
  real(RNP) :: y

  y = 0

end function IntFuncFour

end program Tensor_Products
