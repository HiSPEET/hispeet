!> summary:  Test of L2 projection and interpolation for a single element
!> author:   Joerg Stiller
!> date:     2024/12/10
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

program Test_Projection
  use Kind_Parameters
  use Constants
  use Gauss_Jacobi
  use Erfc_Log_Filter
  use Exponential_Filter
  use Standard_Element_Operators__1D
  use Embedded_Interpolation_Operator__1D
  use Projection_Operator__1D
  implicit none

  real(RNP)    :: c_front    =   5  ! front parameter (higher --> steeper)
  integer      :: n_s        = 200  ! number of sampling intervals
  integer      :: po_c       =  16  ! test function polynomial order
  character(2) :: nodes_c    =  'L' ! test function type of nodes
  integer      :: po_p       =   8  ! polynomial order of projected function
  character(2) :: nodes_p    =  'L' ! node type of projected function
  character    :: projection =  'P' ! projection method {'I','P'}
  integer      :: filter     =   0  ! 0/1/2: none/erfc-log/exponential filter
  real(RNP)    :: po_f       =   1  ! filter order > 0

  namelist/input/ c_front, n_s, po_c, nodes_c, po_p, nodes_p, projection, &
                  filter, po_f

  type(StandardElementOperators_1D)      :: eop_c
  type(StandardElementOperators_1D)      :: eop_p
  type(ProjectionOperator_1D)            :: pop_cp
  type(EmbeddedInterpolationOperator_1D) :: iop_cp
  type(EmbeddedInterpolationOperator_1D) :: iop_cs
  type(EmbeddedInterpolationOperator_1D) :: iop_ps

  real(RNP), allocatable :: x_s (:) ! sampling point coordinates
  real(RNP), allocatable :: f_s (:) ! exact function at sampling points
  real(RNP), allocatable :: f_c (:) ! exact function at collocation points
  real(RNP), allocatable :: f_cs(:) ! interpolant of f_c at sampling points
  real(RNP), allocatable :: f_f (:) ! filtered function at collocation points
  real(RNP), allocatable :: f_fs(:) ! filtered function at sampling points
  real(RNP), allocatable :: f_p (:) ! projected function at basis nodes
  real(RNP), allocatable :: f_ps(:) ! interpolant of f_p at sampling points

  real(RNP), allocatable :: VL    (:,:) ! Legendre Vandermonde matrix
  real(RNP), allocatable :: VL_inv(:,:) ! inverse Legendre Vandermonde matrix
  real(RNP), allocatable :: A_f   (:,:) ! filtering matrix
  real(RNP) :: sigma, theta

  logical :: exists
  integer :: io
  integer :: i, j, k

  ! read parameters
  inquire(file='test_projection.prm', exist=exists)
  if (exists) then
    open(newunit=io, file='test_projection.prm')
    read(io,input)
    close(io)
  end if

  ! standard operators
  eop_c = StandardElementOperators_1d(po_c, nodes_c) ! using collocation points
  eop_p = StandardElementOperators_1d(po_p, nodes_p) ! using projection nodes

  ! exact function at sampling points
  allocate(x_s(0:n_s), source = [(real(TWO*i/n_s - ONE, RNP), i = 0, n_s)])
  allocate(f_s(0:n_s), source = real(atan(c_front * x_s), RNP))

  ! interpolation operators
  iop_cs = EmbeddedInterpolationOperator_1D(eop_c, x_s)
  iop_ps = EmbeddedInterpolationOperator_1D(eop_p, x_s)

  ! exact function at collocation points
  allocate(f_c(0:po_c), source = real(atan(c_front * eop_c%x), RNP))

  ! interpolant of f_c at sampling points
  allocate(f_cs(0:n_s), source = matmul(iop_cs%A, f_c))

  if (filter > 0 .and. po_f > 0) then

    ! filtering
    open(newunit = io, file = 'test_projection_filter.dat')
    write(io,'(A)') '#   theta       sigma'
    allocate(A_f(0:po_c,0:po_c), source = ZERO)
    allocate(VL, VL_inv,  mold = A_f)
    call eop_c % Get_Legendre_VDM(VL)
    call eop_c % Get_Inverse_Legendre_VDM(VL_inv)
    do k = 0, po_c
      theta = real(k, RNP) / po_c
      select case(filter)
      case(1)
        sigma = ErfcLogFilter(theta, po_f)
      case(2)
        sigma = ExponentialFilter(theta, po_f)
      end select
      write(io,'(2(ES12.5,2X))') theta, sigma
      do j = 0, po_c
      do i = 0, po_c
        A_f(i,j) = A_f(i,j) + VL(i,k) * sigma * VL_inv(k,j)
      end do
      end do
    end do
    close(io)

    allocate(f_f(0:po_c), source = matmul(A_f, f_c))
    allocate(f_fs(0:n_s), source = matmul(iop_cs%A, f_f))

  else

    allocate(f_f,  source = f_c )
    allocate(f_fs, source = f_cs)

  end if

  ! projected function at basis nodes
  select case(projection)
  case('I')
    ! embedded intepolation
    iop_cp = EmbeddedInterpolationOperator_1D(eop_c, eop_p % x)
    allocate(f_p(0:po_p), source = matmul(iop_cp%A, f_f))
  case default
    ! L2-projection
    pop_cp = ProjectionOperator_1D(eop_p, eop_c % x, eop_c % nodes)
    allocate(f_p(0:po_p), source = matmul(pop_cp%A, f_f))
  end select

  ! interpolant of f_p at sampling points
  allocate(f_ps(0:n_s), source = matmul(iop_ps%A, f_p))

  ! control output
  write(*,'(2X,A,ES10.3)') 'interpolation error =', maxval(abs(f_cs - f_s))
  write(*,'(2X,A,ES10.3)') 'projection    error =', maxval(abs(f_ps - f_cs))

  ! output
  open(newunit = io, file = 'test_projection.dat')
  write(io,'(A,9(8X,A,2X))') '#', '  x', '   f', 'f_cs', 'f_fs', 'f_ps', &
                             ' e_i', ' e_p', ' e_f', 'e_fp'
  do i = 0, n_s
    write(io,'(9(ES12.5,2X))') x_s(i),            &
                               f_s(i), f_cs(i),   &
                               f_fs(i), f_ps(i),  &
                               f_cs(i) - f_s(i),  &
                               f_ps(i) - f_cs(i), &
                               f_fs(i) - f_s(i),  &
                               f_fs(i) - f_cs(i)
  end do
  close(io)

  !=============================================================================

end program Test_Projection
