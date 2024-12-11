!> summary:  Test of L2 projection for a single element
!> author:   Joerg Stiller
!> date:     2024/12/10
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

program Test_L2_Projection
  use Kind_Parameters
  use Constants
  use Gauss_Jacobi
  use Standard_Element_Operators__1D
  use Embedded_Interpolation_Operator__1D
  use Projection_Operator__1D
  implicit none

  integer      :: po_b    =  4  ! basis polynomial order
  character(2) :: nodes_b = 'L' ! basis type of nodes
  integer      :: po_c    = 11  ! test function polynomial order
  character(2) :: nodes_c = 'L' ! test function type of nodes
  integer      :: n_s     = 80  ! number of sampling intervals

  namelist/input/ po_b, nodes_b, po_c, nodes_c, n_s

  type(StandardElementOperators_1D)      :: eop_b
  type(StandardElementOperators_1D)      :: eop_c
  type(ProjectionOperator_1D)            :: pop_cb
  type(EmbeddedInterpolationOperator_1D) :: iop_bs
  type(EmbeddedInterpolationOperator_1D) :: iop_cs

  real(RNP), allocatable :: x_s(:) ! sampling point coordinates
  real(RNP), allocatable :: f_s(:) ! exact     function at sampling points
  real(RNP), allocatable :: f_c(:) ! original  function at collocation points
  real(RNP), allocatable :: f_i(:) ! original  function at sampling points
  real(RNP), allocatable :: f_b(:) ! projected function at basis nodes
  real(RNP), allocatable :: f_p(:) ! projected function at sampling points
  logical :: exists
  integer :: i, io

  ! read parameters
  inquire(file='test_l2_projection.prm', exist=exists)
  if (exists) then
    open(newunit=io, file='test_l2_projection.prm')
    read(io,input)
    close(io)
  end if

  ! standard operators
  eop_b = StandardElementOperators_1d(po_b, nodes_b) ! using basis nodes
  eop_c = StandardElementOperators_1d(po_c, nodes_c) ! using collocation points

  ! exact function at sampling points
  allocate(x_s(0:n_s), source = [(real(TWO*i/n_s - ONE, RNP), i = 0, n_s)])
  allocate(f_s(0:n_s), source = real(sin(PI * x_s), RNP))

  ! projection and interpolation operators
  pop_cb = ProjectionOperator_1D(eop_b, eop_c % x, eop_c % nodes)
  iop_bs = EmbeddedInterpolationOperator_1D(eop_b, x_s)
  iop_cs = EmbeddedInterpolationOperator_1D(eop_c, x_s)

  ! original function at collocation points
  allocate(f_c(0:po_c), source = real(sin(PI * eop_c%x), RNP))

  ! projected function at basis nodes
  allocate(f_b(0:po_b), source = matmul(pop_cb%A, f_c))

  ! interpolation to sampling points
  allocate(f_i(0:n_s), source = matmul(iop_cs%A, f_c))
  allocate(f_p(0:n_s), source = matmul(iop_bs%A, f_b))

  ! control output
  write(*,'(2X,A,ES10.3)') 'interpolation error =', maxval(abs(f_i - f_s))
  write(*,'(2X,A,ES10.3)') 'projection    error =', maxval(abs(f_p - f_i))

  ! output
  open(newunit = io, file = 'test_l2_projection.dat')
  write(io,'(A)') '#    x             f             f_i           f_p'
  do i = 0, n_s
    write(io,'(4(ES12.5,2X))') x_s(i), f_s(i), f_i(i), f_p(i)
  end do
  close(io)

end program Test_L2_Projection
