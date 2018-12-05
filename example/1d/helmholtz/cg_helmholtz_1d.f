!> summary:  Spectral element solver for the 1D Helmholtz equation
!> author:   Joerg Stiller
!> date:     2016/10/28
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Continuous Galerkin spectral element solver for the 1D Helmholtz equation
!>
!> Solves the Helmholtz equation
!>
!>     u - lambda u" = f(x),    lambda = 1
!>
!> in the domain (-1,1) with u(-1) and u'(1) given. Test cases are based on the
!> exact solution
!>
!>     u = sin(PI x), or
!>     u = abs(x)^3
!>
!> Discretization is performed using continuous nodal spectral elements based on
!> GLL points. Available solution methods are
!>
!>   *  Conjugate gradients (CG)
!>   *  Static condensation + tridiagonal Gauss elimination (SC+GE)
!>
!> Implementation largely follows the approach described in the course "Höhere
!> Numerische Strömungsmechanik" given at TU Dresden from 2016 on. The CG method
!> was adopted from J.R. Shewchuk, An Introduction to the Conjugate Gradient
!> Method Without the Agonizing Pain, Carnegie Mellon University, 1994.
!> For the treatment of the singular problem, i.e. Poisson with Neumann BC, see
!> E.F. Kaasschieter, Preconditioned conjugate gradients for solving singular
!> systems J. Comput. Appl. Math., 1988, 24, 265-275.
!>
!> @note
!> This program is written for readabilty, not speed! Serious optimization would
!> require to restrain from using array intrinsics in favor of simple loops or
!> tuned library routines.
!> @endnote
!>
!===============================================================================

program CG_Helmholtz_1D
  use Kind_Parameters,  only: RNP, IXL
  use Constants,        only: PI, ONE, TWO
  use CG_Element_Operators_1D
  use CG_Utilities_1D
  use CG_Condensed_Solver_1D

  implicit none

  ! variables ..................................................................

  ! problem parameters
  real(RNP) :: lambda  = 1      ! Helmholtz parameter
  integer   :: problem = 1      ! test case
  integer   :: init    = 1      ! intial conditions (0: zero, 1: random)
  character :: bc(2)   = 'D'    ! left/right BC ('D': Dirichlet, 'N': Neumann)

  namelist /problem_parameters/ lambda, problem, init, bc

  ! solution parameters
  integer   :: po     = 16      ! polynomial order
  integer   :: ne     = 10      ! number of elements
  integer   :: method = 1       ! solution method (1: CG, 2: SC+GE)
  integer   :: i_max  = huge(1) ! maximum number of CG iterations
  real(RNP) :: r_max  = 1E-12   ! maximum CG residual

  namelist /solution_parameters/ po, ne, method, r_max, i_max

  ! discrete variables and operators
  type(CG_ElementOperators1D) :: eop       ! element operators
  real(RNP), allocatable      :: x(:,:)    ! mesh points
  real(RNP), allocatable      :: u(:,:)    ! discrete solution
  real(RNP), allocatable      :: f(:,:)    ! right hand side (RHS)
  real(RNP), allocatable      :: r(:,:)    ! residual
  real(RNP), allocatable      :: s(:,:)    ! projected exact solution
  real(RNP), allocatable      :: e(:,:)    ! error
  real(RNP), allocatable      :: w(:,:)    ! node weights
  real(RNP), allocatable      :: Me(:)     ! element mass matrix
  real(RNP), allocatable      :: He(:,:)   ! element Helmholtz matrix

  ! auxiliary variables
  logical      :: exists, singular
  integer      :: i, l, n, io
  integer(IXL) :: count0, count1, count_rate
  real(RNP)    :: dx, t_pre, t_sol

  ! initialization .............................................................

  ! greeting
  write(*,'(/,A)') 'Spectral element solver for 1D Helmholtz equation'

  ! load parameters
  inquire(file='cg_helmholtz_1d.prm', exist=exists)
  if (exists) then
    open(newunit=io, file='cg_helmholtz_1d.prm')
    read(io, nml=problem_parameters)
    read(io, nml=solution_parameters)
    close(io)
  end if

  ! start system clock
  call system_clock(count0, count_rate = count_rate)

  ! workspace
  allocate( x(0:po,ne), u(0:po,ne), f(0:po,ne), r(0:po,ne), &
            s(0:po,ne), e(0:po,ne), w(0:po,ne)              )

  ! element operators
  call eop % New(po)
  if (method == 2) then
    call eop % BuildInteriorEigensystem()
  end if

  ! mesh and point weights
  call GetMeshPoints(eop, -ONE, ONE, dx, x)
  call GetPointWeights(bc, w)

  ! element operators
  allocate(Me(0:po), He(0:po,0:po))
  call GetElementOperators(eop, dx, lambda, Me, He)

  ! check if problem is singular
  singular = lambda == 0 .and. (all(bc == 'N') .or. all(bc == 'P'))

  ! right hand side
  call GetRHS(Me, x, bc, f)
  if (singular) then ! project f to nullspace
    f = f - sum(w*f) / sum(w)
  end if

  ! initial values
  if (init == 0) then ! start from zero
    u = 0
  else ! intial guess, chosen at random from [0,1]
    call random_number(u)
    call MakeContinuous(bc, u)
  end if

  ! inject Dirichlet BC
  if (bc(1) == 'D')  u( 0,  1) = u_exact(-ONE)
  if (bc(2) == 'D')  u(po, ne) = u_exact( ONE)

  ! read system clock
  call system_clock(count1)

  ! time consumed with preprocessing
  t_pre = (count1 - count0) / real(count_rate, RNP)

  ! solution ...................................................................

  call system_clock(count0)

  select case(method)
  case(1) ! CG
    call CG(He, bc, u, f, w, r_max, i_max)
  case(2) ! static condensation + Gauss elimination, nu = ONE
    call CondensedEllipticSolver(eop, dx, lambda, ONE, bc, f, u)
  end select

  call system_clock(count1)
  t_sol = (count1 - count0) / real(count_rate, RNP)

  ! evaluation .................................................................

  ! number of unknowns
  n = po * ne

  ! consistency error
  s = u_exact(x)
  call HelmholtzResidual(He, bc, s, f, r)
  write(*,'(/,A)') 'consistency error'
  write(*,'(2X,A,ES12.5)') 'c_max =', maxval(abs(r))
  write(*,'(2X,A,ES12.5)') 'c_rms =', sqrt(sum(w*r*r)/n)

  ! final residual
  call HelmholtzResidual(He, bc, u, f, r)
  write(*,'(/,A)') 'final residual'
  write(*,'(2X,A,ES12.5)') 'r_max =', maxval(abs(r))
  write(*,'(2X,A,ES12.5)') 'r_rms =', sqrt(sum(w*r*r)/n)

  ! error
  e = u - s
  ! remove constant in singular case
  if (singular) then
    e = e - (maxval(e) + minval(e))/2
  end if
  write(*,'(/,A)') 'error'
  write(*,'(2X,A,ES12.5)') 'e_max =', maxval(abs(e))
  write(*,'(2X,A,ES12.5)') 'e_rms =', sqrt(sum(w*e*e)/n)

  ! timing
  write(*,'(/,A)') 'runtime'
  write(*,'(2X,A,ES12.5)') 't_pre =', t_pre
  write(*,'(2X,A,ES12.5)') 't_sol =', t_sol

  ! save results
  open(newunit=io, file='helmholtz_sem_1d.dat')
  write(io,'(10(A17,1X))') '# x', 'u', 's', 'e', 'f', 'r'
  do l = 1, ne
  do i = 0, po
    write(io,'(10(ES17.10,1X))') x(i,l), u(i,l), s(i,l), e(i,l), f(i,l), r(i,l)
  end do
  end do
  close(io)

contains

!-------------------------------------------------------------------------------
!> Exact solution

elemental real(RNP) function u_exact(x) result(u)
  real(RNP), intent(in) :: x !< point coordinate

  select case(problem)
  case(1)
    u = sin(PI*x)
  case(2)
    u = abs(x)**3
  case default
    u = 0
  end select

end function u_exact

!-------------------------------------------------------------------------------
!> Exact first derivative

elemental real(RNP) function du_exact(x) result(du)
  real(RNP), intent(in) :: x !< point coordinate

  select case(problem)
  case(1)
    du = PI * cos(PI*x)
  case(2)
    du = 3 * abs(x) * x
  case default
    du = 0
  end select

end function du_exact

!-------------------------------------------------------------------------------
!> Exact second derivative

elemental real(RNP) function ddu_exact(x) result(ddu)
  real(RNP), intent(in) :: x !< point coordinate

  select case(problem)
  case(1)
    ddu = -PI**2 * sin(PI*x)
  case(2)
    ddu = 6 * abs(x)
  case default
    ddu = 0
  end select

end function ddu_exact

!-------------------------------------------------------------------------------
!> Element operators

subroutine GetElementOperators(eop, dx, lambda, Me, He)
  class(CG_ElementOperators1D), intent(in) :: eop !< standard operators
  real(RNP), intent(in)  :: dx        !< element length
  real(RNP), intent(in)  :: lambda    !< Helmholtz parameter
  real(RNP), intent(out) :: Me(0:)    !< element mass matrix (main diagonal)
  real(RNP), intent(out) :: He(0:,0:) !< element Helmholtz operator

  integer :: i

  associate(po => eop % po, Ms => eop % w, Ls => eop % L)

    ! element mass matrix
    Me = dx/2 * Ms

    ! element Helmhotz operator
    He = 2/dx * Ls
    if (lambda > 0) then
      do i = 0, po
        He(i,i) = He(i,i) + lambda * Me(i)
      end do
    end if

  end associate

end subroutine GetElementOperators

!-------------------------------------------------------------------------------
!> Right hand side

subroutine GetRHS(Me, x, bc, f)
  real(RNP), intent(in)  :: Me(0:)  !< element mass matrix (main diagonal)
  real(RNP), intent(in)  :: x(0:,:) !< mesh points
  character, intent(in)  :: bc(2)   !< boundary conditions
  real(RNP), intent(out) :: f(0:,:) !< RHS

  integer :: l, ne

  po = ubound(x,1)
  ne = ubound(x,2)

  ! projection of the source term
  do l = 1, ne
     f(:,l) = Me * (lambda*u_exact(x(:,l)) - ddu_exact(x(:,l)))
  end do

  ! Neumann BC
  if (bc(1) == 'N')  f( 0,  1) = f( 0,  1) - du_exact(-ONE)
  if (bc(2) == 'N')  f(po, ne) = f(po, ne) + du_exact( ONE)

  ! assemble element contributions
  call Assembly(bc, f)

end subroutine GetRHS

!------------------------------------------------------------------------------
!> Residual of a given approximate solution

subroutine HelmholtzResidual(He, bc, u, f, r)
  real(RNP),           intent(in)  :: He(0:,0:) !< element Helmholtz operator
  character,           intent(in)  :: bc(2)     !< boundary conditions
  real(RNP),           intent(in)  :: u(0:,:)   !< approximate solution
  real(RNP), optional, intent(in)  :: f(0:,:)   !< RHS, default: f = 0
  real(RNP),           intent(out) :: r(0:,:)   !< residual, r = f - Au

  ! element contributions
  r = -matmul(He, u)

  ! assembly
  call Assembly(bc, r)

  ! add RHS contribution
  if (present(f)) then
    r = r + f
  end if

  ! nullify residual at Dirichlet points
  if (bc(1) == 'D')  r( 0,  1) = 0
  if (bc(2) == 'D')  r(po, ne) = 0

end subroutine HelmholtzResidual

!------------------------------------------------------------------------------
!> Conjugate gradient method

subroutine CG(He, bc, u, f, w, r_max, i_max)
  real(RNP), intent(in)    :: He(0:,0:) !< element Helmholtz operator
  character, intent(in)    :: bc(2)     !< boundary conditions
  real(RNP), intent(inout) :: u(0:,:)   !< approximate solution
  real(RNP), intent(in)    :: f(0:,:)   !< right hand side
  real(RNP), intent(in)    :: w(0:,:)   !< node weights
  real(RNP), intent(in)    :: r_max     !< target residual
  integer,   intent(in)    :: i_max     !< max number of iterations

  real(RNP), dimension(:,:), allocatable :: r, p, q
  real(RNP) :: alpha, beta, delta, delta_old
  integer   :: i, po, ne, ni

  ! initialization .............................................................

  ! dimensions
  po = ubound(u,1)
  ne = ubound(u,2)

  ! number of iterations, with safety factor of 10
  ni = 10 * ne * po*po

  ! workspace
  allocate(r(0:po,1:ne), p(0:po,1:ne), q(0:po,1:ne))

  ! initial residual
  call HelmholtzResidual(He, bc, u, f, r)

  ! iteration ..................................................................

  delta = sum(w * r * r)
  p = r
  do i = 1, ni
    if (delta <= r_max**2 .or. i > i_max) exit
    call HelmholtzResidual(He, bc, p, r=q)
    alpha = -delta / sum(w * p * q)
    u = u + alpha * p
    if (mod(i,50) == 0) then
      call HelmholtzResidual(He, bc, u, f, r)
    else
      r = r + alpha * q
    end if
    delta_old = delta
    delta = sum(w * r * r)
    beta  = delta / delta_old
    p = r + (delta / delta_old) * p
  end do

end subroutine CG

!===============================================================================

end program CG_Helmholtz_1D
