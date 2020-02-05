!> summary:  Spectral element solver for the 1D Helmholtz equation
!> author:   Joerg Stiller
!> date:     2016/10/28
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Continuous Galerkin spectral element solver for the 1D Helmholtz equation
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
  use Constants,        only: ONE, ZERO
  use CG_Element_Operators_1D
  use CG_Utilities_1D
  use CG_Condensed_Solver_1D
  use Helmholtz_1D_Test_Cases

  implicit none

  ! variables ..................................................................

  ! problem parameters
  real(RNP) :: lambda  = 1      ! Helmholtz parameter
  real(RNP) :: nu      = 1      ! physical diffusivity
  integer   :: test    = 1      ! test case
  integer   :: init    = 1      ! intial conditions (0: zero, 1: random)
  character :: bc(2)   = 'D'    ! left/right BC ('D': Dirichlet, 'N': Neumann)

  namelist /problem_parameters/ lambda, nu, test, init, bc

  ! solution parameters
  type(CG_ElementOptions1D) :: eop_opt ! options for CG element operator
  integer   :: ne     = 10      ! number of elements
  integer   :: method = 1       ! solution method (1: CG, 2: SC+GE)
  integer   :: i_max  = huge(1) ! maximum number of CG iterations
  real(RNP) :: r_max  = 1E-12   ! maximum CG residual
  real(RNP) :: nu_svv = -1      ! spectral diffusivity amplitude

  namelist /solution_parameters/ eop_opt, ne, method, r_max, i_max, nu_svv

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
  logical      :: exists, periodic, singular
  integer      :: po, i, l, n, io
  integer(IXL) :: count0, count1, count_rate
  real(RNP)    :: dx, t_pre, t_sol

  ! initialization .............................................................

  ! greeting
  write(*,'(/,A)') 'Spectral element solver for 1D Helmholtz equation'

  ! load parameters
  inquire(file='cg_helmholtz_1d.prm', exist=exists)
  if (exists) then
    open(newunit=io, file='cg_helmholtz_1d.prm')
    read(io, nml=problem_parameters )
    read(io, nml=solution_parameters)
    close(io)
  end if

  po = eop_opt%po
  ! set spectral diffusivity to default value if not given as parameter
  if (nu_svv == -1) nu_svv = ONE / real(po,RNP)
  if (.not. svv) nu_svv = ZERO

  call SetTestCase(test)
  periodic = all(bc == 'P')

  ! start system clock
  call system_clock(count0, count_rate = count_rate)

  ! workspace
  allocate( x(0:po,ne), u(0:po,ne), f(0:po,ne), r(0:po,ne), &
            s(0:po,ne), e(0:po,ne), w(0:po,ne)              )

  ! standard operators
  eop = CG_ElementOperators1D(eop_opt)

  ! mesh and point weights
  call GetMeshPoints(eop, -ONE, ONE, dx, x)
  call GetPointWeights(w, periodic)

  ! element operators
  allocate(Me(0:po), He(0:po,0:po))
  call GetElementOperators(eop, dx, lambda, nu, nu_svv, Me, He)

  ! check if problem is singular
  singular = lambda == 0 .and. (all(bc == 'N') .or. all(bc == 'P'))

  ! right hand side
  call GetRHS(Me, bc, x, lambda, f)
  if (singular) then ! project f to nullspace
    f = f - sum(w*f) / sum(w)
  end if

  ! initial values
  if (init == 0) then ! start from zero
    u = 0
  else ! intial guess, chosen at random from [0,1]
    call random_number(u)
    call MakeContinuous(u, periodic)
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
  case(2) ! static condensation + Gauss elimination
    call CondensedEllipticSolver(eop, dx, lambda, nu, nu_svv, bc, f, u)
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
  open(newunit=io, file='cg_helmholtz_1d.dat')
  write(io,'(10(A17,1X))') '# x', 'u', 's', 'e', 'f', 'r'
  do l = 1, ne
  do i = 0, po
    write(io,'(10(ES17.10,1X))') x(i,l), u(i,l), s(i,l), e(i,l), f(i,l), r(i,l)
  end do
  end do
  close(io)

contains

!-------------------------------------------------------------------------------
!> Element operators

subroutine GetElementOperators(eop, dx, lambda, nu, nu_svv, Me, He)
  class(CG_ElementOperators1D), intent(in) :: eop !< standard operators
  real(RNP), intent(in)  :: dx        !< element length
  real(RNP), intent(in)  :: lambda    !< Helmholtz parameter
  real(RNP), intent(in)  :: nu        !< diffusivity
  real(RNP), intent(in)  :: nu_svv    !< spectral diffusivity amplitude
  real(RNP), intent(out) :: Me(0:)    !< element mass matrix (main diagonal)
  real(RNP), intent(out) :: He(0:,0:) !< element Helmholtz operator

  integer :: i, po
  real(RNP), allocatable :: Ls(:,:)

  po = eop % po

  ! standard diffusion matrix comprising regular and SVV contributions
  allocate(Ls(0:po,0:po), source = ZERO)
  if (eop % Has_SVV()) then
    call eop % Get_SVV_StandardStiffnessMatrix(Ls)
  end if
  Ls = nu * eop%L + nu_svv * Ls

  associate(Ms => eop % w)

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

subroutine GetRHS(Me, bc, x, lambda, f)
  real(RNP), intent(in)  :: Me(0:)  !< element mass matrix (main diagonal)
  character, intent(in)  :: bc(2)   !< boundary conditions
  real(RNP), intent(in)  :: x(0:,:) !< mesh points
  real(RNP), intent(in)  :: lambda  ! Helmholtz parameter
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
  call Assembly(f, periodic = bc(1)=='P')

  ! nullify RHS at Dirichlet boundaries
  if (bc(1) == 'D')  f( 0,  1) = 0
  if (bc(2) == 'D')  f(po, ne) = 0

end subroutine GetRHS

!------------------------------------------------------------------------------
!> Application of Helmholtz operator, `v = Au`, `v = 0` at Dirichlet boundaries

subroutine HelmholtzOperator(He, bc, u, v)
  real(RNP), intent(in)  :: He(0:,0:) !< element Helmholtz operator
  character, intent(in)  :: bc(2)     !< boundary conditions
  real(RNP), intent(in)  :: u(0:,:)   !< approximate solution
  real(RNP), intent(out) :: v(0:,:)   !< residual, r = f - Au

  ! element contributions
  v = matmul(He, u)

  ! assembly
  call Assembly(v, periodic = bc(1)=='P')

  ! nullify result at Dirichlet points
  if (bc(1) == 'D')  v( 0,  1) = 0
  if (bc(2) == 'D')  v(po, ne) = 0

end subroutine HelmholtzOperator

!------------------------------------------------------------------------------
!> Residual of a given approximate solution

subroutine HelmholtzResidual(He, bc, u, f, r)
  real(RNP), intent(in)  :: He(0:,0:) !< element Helmholtz operator
  character, intent(in)  :: bc(2)     !< boundary conditions
  real(RNP), intent(in)  :: u(0:,:)   !< approximate solution
  real(RNP), intent(in)  :: f(0:,:)   !< RHS
  real(RNP), intent(out) :: r(0:,:)   !< residual, r = f - Au

  call HelmholtzOperator(He, bc, u, r)
  r = f - r

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

  ! max number of iterations
  ni = min(i_max, 10 * ne * po*po)

  ! workspace
  allocate(r(0:po,1:ne), p(0:po,1:ne), q(0:po,1:ne))

  ! initial residual, providing r = 0 in Dirichlet points
  call HelmholtzResidual(He, bc, u, f, r)

  ! iteration ..................................................................

  ! Note that scalar products are weighted because r,p,q represent local arrays
  ! with interior element boundary values listed twice.

  delta = sum(w * r * r)
  p = r
  do i = 1, ni
    if (delta <= r_max**2) exit
    ! q = Ap, 0 in Dirichlet points
    call HelmholtzOperator(He, bc, p, q)
    alpha = delta / sum(w * p * q)
    delta_old = delta
    u = u + alpha * p
    if (mod(i,50) == 0) then
      ! r = f - Au, evaluated explicitly eliminate accumulated round-off errors
      call HelmholtzResidual(He, bc, u, f, r)
    else
      r = r - alpha * q
    end if
    delta = sum(w * r * r)
    beta  = delta / delta_old
    p = r + beta * p
  end do

end subroutine CG

!===============================================================================

end program CG_Helmholtz_1D
