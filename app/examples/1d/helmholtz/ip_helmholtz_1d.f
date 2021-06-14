!> summary:  Hybridized IP/DG solver for the 1D Helmholtz equation
!> author:   Joerg Stiller
!> date:     2019/01/23
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Hybridized IP/DG solver for the 1D Helmholtz equation
!>
!> Solves the Helmholtz equation
!>
!>     lambda u - u" = f(x),    lambda = 1
!>
!> in the domain (-1,1) with u(-1) and u'(1) given. Test cases are based on the
!> exact solution
!>
!>     1)  u = sin(PI x), or
!>     2)  u = abs(x)^3
!>
!> Discretization is performed using the hybridized symmetric interior penalty
!> method with a nodal (Lobatto) basis.
!>
!===============================================================================

program IP_Helmholtz_1D
  use Kind_Parameters,  only: RNP, IXL
  use Constants,        only: ONE, TWO, ZERO
  use IP_Element_Operators__1D
  use IP_Hybridized_Solver_1D
  use Helmholtz_1D_Test_Cases

  implicit none

  ! variables ..................................................................

  ! problem parameters
  real(RNP) :: lambda  =  1   ! Helmholtz parameter
  real(RNP) :: nu      =  1   ! diffusivity
  real(RNP) :: nu_svv  = -1   ! spectral diffusivity amplitude
  integer   :: test    =  1   ! test case
  character :: bc(2)   = 'D'  ! left/right BC ('D': Dirichlet, 'N': Neumann)

  namelist /problem_parameters/ lambda, nu, nu_svv, test, bc

  ! solution parameters
  type(IP_ElementOptions_1D) :: eop_opt ! options for IP element operator
  integer                   :: ne = 10 ! number of elements

  namelist /solution_parameters/ eop_opt, ne

  ! discrete variables and operators
  type(IP_ElementOperators_1D) :: eop    ! element operators
  real(RNP), allocatable      :: x(:,:) ! mesh points
  real(RNP), allocatable      :: u(:,:) ! discrete solution
  real(RNP), allocatable      :: f(:,:) ! right hand side (RHS)
  real(RNP), allocatable      :: s(:,:) ! projected exact solution
  real(RNP), allocatable      :: e(:,:) ! error

  ! auxiliary variables
  logical      :: exists, singular
  integer      :: po, i, l, n, io
  integer(IXL) :: count0, count1, count_rate
  real(RNP)    :: dx, t_pre, t_sol

  ! initialization .............................................................

  ! greeting
  write(*,'(/,A)') 'Hybridized IP/DG-SEM solver for 1D Helmholtz equation'

  ! load parameters
  inquire(file='ip_helmholtz_1d.prm', exist=exists)
  if (exists) then
    open(newunit=io, file='ip_helmholtz_1d.prm')
    read(io, nml=problem_parameters )
    read(io, nml=solution_parameters)
    close(io)
  end if

  eop_opt%hybrid = .true.
  po = eop_opt%po
  ! set spectral diffusivity to default value if not given as parameter
  if (nu_svv < 0) nu_svv = ONE / real(po,RNP)

  call SetTestCase(test)

  ! start system clock
  call system_clock(count0, count_rate = count_rate)

  ! check if problem is singular
  singular = lambda == 0 .and. (all(bc == 'N') .or. all(bc == 'P'))

  ! workspace
  allocate(x(0:po,ne), u(0:po,ne), f(0:po,ne), s(0:po,ne), e(0:po,ne))
  u = 0
  n = size(u)

  ! standard operators
  eop = IP_ElementOperators_1D(eop_opt)

  ! mesh
  dx = TWO / ne
  do i = 1, ne
     x(:,i) = (2*i - 1 + eop%x) * dx/2 - 1
  end do

  ! right hand side
  call GetRHS(eop, dx, bc, x, f)
  if (singular) then ! project f to nullspace
    f = f - sum(f) / n
  end if

  ! read system clock
  call system_clock(count1)

  ! time consumed with preprocessing
  t_pre = (count1 - count0) / real(count_rate, RNP)

  ! solution ...................................................................

  call system_clock(count0)
  if (eop % Has_SVV()) then
    call HybridEllipticSolver(eop, dx, lambda, nu, nu_svv, bc, f, u)
  else
    call HybridEllipticSolver(eop, dx, lambda, nu        , bc, f, u)
  end if
  call system_clock(count1)

  t_sol = (count1 - count0) / real(count_rate, RNP)

  ! evaluation .................................................................

  ! error
  s = u_exact(x)
  e = u - s
  ! remove constant in singular case
  if (singular) then
    e = e - (maxval(e) + minval(e))/2
  end if
  write(*,'(/,A)') 'error'
  write(*,'(2X,A,ES12.5)') 'e_max =', maxval(abs(e))
  write(*,'(2X,A,ES12.5)') 'e_rms =', sqrt(sum(e*e)/n)

  ! timing
  write(*,'(/,A)') 'runtime'
  write(*,'(2X,A,ES12.5)') 't_pre =', t_pre
  write(*,'(2X,A,ES12.5)') 't_sol =', t_sol

  ! save results
  open(newunit=io, file='ip_helmholtz_1d.dat')
  write(io,'(10(A17,1X))') '# x', 'u', 's', 'e', 'f'
  do l = 1, ne
  do i = 0, po
    write(io,'(10(ES17.10,1X))') x(i,l), u(i,l), s(i,l), e(i,l), f(i,l)
  end do
  end do
  close(io)

contains

!-------------------------------------------------------------------------------
!> Right hand side

subroutine GetRHS(eop, dx, bc, x, f)
  class(IP_ElementOperators_1D), intent(in)  :: eop     !< IP-H element operators
  real(RNP),                    intent(in)  :: dx      !< element extension
  character,                    intent(in)  :: bc(2)   !< boundary conditions
  real(RNP),                    intent(in)  :: x(0:,:) !< mesh points
  real(RNP),                    intent(out) :: f(0:,:) !< RHS

  real(RNP), allocatable :: delta_0(:), delta_P(:)
  real(RNP) :: tau
  integer   :: l, ne

  ne = ubound(x,2)

  associate(P => eop%po, Ms => eop%w, Ds => eop%D)

    tau = 2 * eop % PenaltyFactor(dx)

    allocate(delta_0(0:P), source = ZERO)
    delta_0(0) = ONE

    allocate(delta_P(0:P), source = ZERO)
    delta_P(P) = ONE

    ! projection of the source term
    do l = 1, ne
       f(:,l) = dx/2 * Ms * (lambda*u_exact(x(:,l)) - ddu_exact(x(:,l)))
    end do

    ! left boundary
    if (bc(1) == 'D') then
      f(:,1) = f(:,1) + (2/dx * Ds(0,:) + tau * delta_0) * u_exact(x(0,1))
    else if (bc(1) == 'N') then
      f(0,1) = f(0,1) - du_exact(-ONE)
    end if

    ! right boundary
    if (bc(2) == 'D') then
      f(:,ne) = f(:,ne) + (-2/dx * Ds(P,:) + tau * delta_P) * u_exact(x(P,ne))
    else if (bc(2) == 'N') then
      f(P,ne) = f(P,ne) + du_exact(ONE)
    end if

  end associate

end subroutine GetRHS

!===============================================================================

end program IP_Helmholtz_1D
