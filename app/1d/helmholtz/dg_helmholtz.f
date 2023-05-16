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

program DG_Helmholtz
  use Kind_Parameters, only: RNP, IXL
  use Constants,       only: ZERO, ONE, TWO, HALF
  use DG__Element_Operators__1D
  use DG__Elliptic_Operator__1D
  use Helmholtz_Test_Cases

  implicit none

  ! variables ..................................................................

  ! problem parameters
  real(RNP) :: lambda  =  1   ! Helmholtz parameter
  real(RNP) :: nu_p    =  1   ! physical diffusivity
  real(RNP) :: nu_s    =  0   ! spectral diffusivity amplitude
  integer   :: test    =  1   ! test case
  character :: bc(2)   = 'D'  ! left/right BC ('D': Dirichlet, 'N': Neumann)

  namelist /problem_parameters/ lambda, nu_p, nu_s, test, bc

  ! solution parameters
  type(DG_ElementOptions_1D) :: dg_opt      ! DG element operator options
  type(DG_SchwarzOptions_1D) :: schwarz_opt ! Schwarz method options
  integer :: ne    = 10  ! number of elements

  namelist /solution_parameters/ dg_opt, schwarz_opt, ne

  ! solution method:
  !   1   hybridization + tridiagonal Gauss
  !   2   conjugate gradient method
  !   3   additive overlapping Schwarz method
  !   4   Schwarz IPCG method
  integer   :: method = 1
  integer   :: i_max  = 1000  ! max num iterations
  real(RNP) :: r_red  = 1e-8  ! targeted residual reduction
  real(RNP) :: r_max  = 1e-12 ! targeted maximum residual

  namelist /solution_parameters/ method, i_max, r_red, r_max

  ! discrete variables and operators
  type(DG_EllipticOperator_1D) :: elliptic_op
  real(RNP), allocatable       :: u(:,:) ! discrete solution
  real(RNP), allocatable       :: x(:,:) ! mesh points
  real(RNP), allocatable       :: f(:,:) ! right hand side (RHS)
  real(RNP), allocatable       :: r(:,:) ! residual
  real(RNP), allocatable       :: s(:,:) ! projected exact solution
  real(RNP), allocatable       :: e(:,:) ! error

  ! auxiliary variables
  logical      :: exists, singular
  integer      :: i, l, io, ni = 0
  integer(IXL) :: count0, count1, count_rate
  real(RNP)    :: dx, nu, r_nu_s, t_pre, t_res, t_sol
  real(RNP)    :: bv(2) = 0

  ! initialization .............................................................

  ! greeting
  write(*,'(/,A)') 'IP/DG-SEM for 1D Helmholtz equation'

  ! load parameters
  inquire(file='dg_helmholtz.prm', exist=exists)
  if (exists) then
    open(newunit=io, file='dg_helmholtz.prm')
    read(io, nml=problem_parameters )
    read(io, nml=solution_parameters)
    close(io)
  end if

  if (method == 1) then
    dg_opt % hybrid = .true. ! required with hybrid solver
  end if

  nu = nu_p + nu_s
  r_nu_s = nu_s / nu

  call SetTestCase(test)

  ! start system clock
  call system_clock(count0, count_rate = count_rate)

  ! check if problem is singular
  singular = lambda == 0 .and. (all(bc == 'N') .or. all(bc == 'P'))

  ! initialize operator
  elliptic_op = DG_EllipticOperator_1D(dg_opt, schwarz_opt, r_nu_s)

  associate(eop => elliptic_op%eop, po => elliptic_op%eop%po)

    ! workspace
    allocate(u(0:po,ne), x(0:po,ne), f(0:po,ne))
    allocate(r(0:po,ne), s(0:po,ne), e(0:po,ne))
    call random_number(u)
    u = u - HALF
!   u = 0

    ! mesh
    dx = TWO / ne
    do i = 1, ne
      x(:,i) = (2*i - 1 + eop%x) * dx/2 - 1
    end do

    ! exact solution
    s = u_exact(x)

    ! right hand side and boundary values
    call GetRHS(eop, dx, lambda, nu_p, bc, x, f, bv)

    ! read system clock
    call system_clock(count1)

    ! time consumed with preprocessing
    t_pre = (count1 - count0) / real(count_rate, RNP)

    ! consistency error ........................................................

    call system_clock(count0)
    call elliptic_op % Residual(bc, bv, dx, lambda, nu, f, s, r)
    call system_clock(count1)

    t_res = (count1 - count0) / real(count_rate, RNP)

    ! solution .................................................................

    call system_clock(count0)
    select case(method)
    case(1)
      call elliptic_op % HybridSolver(bc, bv, dx, lambda, nu, f, u)
    case(2)
      call elliptic_op % CG_Method( bc, bv, dx, lambda, nu, f, u &
                                  , i_max, r_red, r_max, ni = ni )
    case(3)
      call elliptic_op % Schwarz_Method( bc, bv, dx, lambda, nu, f, u &
                                       , i_max, r_red, r_max, ni = ni )
    case(4)
      call elliptic_op % SchwarzPCG_Method( bc, bv, dx, lambda, nu, f, u &
                                          , i_max, r_red, r_max, ni = ni )
    end select
    call system_clock(count1)

    t_sol = (count1 - count0) / real(count_rate, RNP)

    ! evaluation ...............................................................

    ! error
    e = u - s
    ! remove constant in singular case
    if (singular) then
      e = e - (maxval(e) + minval(e))/2
    end if

    write(*,*)
    write(*,'(A)', advance = 'NO') 'test case: '
    select case(test)
    case(1)
      write(*,'(A)') 'u = sin(πx)'
    case(2)
      write(*,'(A)') 'u = |x|³'
    end select

    write(*,'(A)', advance = 'NO') 'solver:    '
    select case(method)
    case(1)
      write(*,'(A)') 'hybridization + tridiagonal Gauss'
    case(2)
      write(*,'(A)') 'conjugate gradient method'
    case(3)
      write(*,'(A)') 'additive overlapping Schwarz method'
    case(4)
      write(*,'(A)') 'Schwarz IPCG method'
    end select

    write(*,'(/,A)') 'consistency'
    write(*,'(2X,A,ES12.5)') 'r_max =', maxval(abs(r))
    write(*,'(2X,A,ES12.5)') 'r_rms =', sqrt(sum(r*r)/size(r))

    if (method > 1) then
      write(*,'(/,A,I0)') 'iterations: ni = ', ni
    end if

    write(*,'(/,A)') 'error'
    write(*,'(2X,A,ES12.5)') 'e_max =', maxval(abs(e))
    write(*,'(2X,A,ES12.5)') 'e_rms =', sqrt(sum(e*e)/size(e))

    ! timing
    write(*,'(/,A)') 'runtime'
    write(*,'(2X,A,ES12.5)') 't_pre =', t_pre
    write(*,'(2X,A,ES12.5)') 't_res =', t_res
    write(*,'(2X,A,ES12.5)') 't_sol =', t_sol

    ! save results
    open(newunit=io, file='dg_helmholtz.dat')
    write(io,'(10(A17,1X))') '# x', 'u', 's', 'e', 'f'
    do l = 1, ne
    do i = 0, po
      write(io,'(10(ES17.10,1X))') x(i,l), u(i,l), s(i,l), e(i,l), f(i,l)
    end do
    end do
    close(io)

  end associate

contains

  !-----------------------------------------------------------------------------
  !> Right hand side

  subroutine GetRHS(eop, dx, lambda, nu, bc, x, f, bv)
    class(DG_ElementOperators_1D), intent(in)  :: eop !< IP-H element operators
    real(RNP), intent(in)  :: dx      !< element extension
    real(RNP), intent(in)  :: lambda  !< Helmholtz parameter
    real(RNP), intent(in)  :: nu      !< diffusivity
    character, intent(in)  :: bc(2)   !< boundary conditions
    real(RNP), intent(in)  :: x(0:,:) !< mesh points
    real(RNP), intent(out) :: f(0:,:) !< RHS
    real(RNP), intent(out) :: bv(2)   !< boundary values

    integer :: l

    associate(po => eop%po, Ms => eop%w)

      ! projection of the source term
      do l = 1, ubound(x,2)
        f(:,l) = dx/2 * Ms * (lambda * u_exact(x(:,l)) - nu*ddu_exact(x(:,l)))
      end do

      ! left boundary
      select case(bc(1))
      case('D')
        bv(1) = u_exact( x(0,1) )
      case('N')
        bv(1) = nu * du_exact( x(0,1) )
      case default
        bv(1) = ZERO
      end select

      ! right boundary
      select case(bc(2))
      case('D')
        bv(2) = u_exact( x(po,ne) )
      case('N')
        bv(2) = nu * du_exact( x(po,ne) )
      case default
        bv(2) = ZERO
      end select

    end associate

  end subroutine GetRHS

  !=============================================================================

end program DG_Helmholtz
