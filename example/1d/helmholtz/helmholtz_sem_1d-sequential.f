!> summary:  Spectral element solver for the 1D Helmholtz equation (sequential)
!> author:   Joerg Stiller
!> date:     2016/10/28
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Spectral element solver for the 1D Helmholtz equation - sequential version
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

program Helmholtz_SEM_1D
  use Kind_Parameters,  only: RNP, IXL
  use Constants,        only: PI, ZERO, ONE, TWO, HALF
  use Matrix_Operators, only: Inverse
  use Linear_Equations, only: TridiagonalSolver
  use Standard_Operators_1D ! provides the SEM standard operators
  implicit none

  ! Variables ..................................................................

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
  type(StandardOperators1D) :: standard_op   ! standard element operators
  real(RNP), allocatable    :: x(:,:)        ! mesh points
  real(RNP), allocatable    :: u(:,:)        ! discrete solution
  real(RNP), allocatable    :: f(:,:)        ! right hand side (RHS)
  real(RNP), allocatable    :: r(:,:)        ! residual
  real(RNP), allocatable    :: s(:,:)        ! projected exact solution
  real(RNP), allocatable    :: e(:,:)        ! error
  real(RNP), allocatable    :: w(:,:)        ! node weights
  real(RNP), allocatable    :: Me(:)         ! element mass matrix
  real(RNP), allocatable    :: He(:,:)       ! element Helmholtz matrix

  ! variables used with the condensed solver
  real(RNP), allocatable :: HBB(:,:)         ! boundary-boundary part of He
  real(RNP), allocatable :: HBI(:,:)         ! boundary-interior part of He
  real(RNP), allocatable :: HII_inv(:,:)     ! inverse of interior-interior part
  real(RNP), allocatable :: Ac(:,:)          ! condensed system matrix
  real(RNP), allocatable :: fc(:)            ! condensed system RHS

  ! auxiliary variables
  logical      :: exists, singular
  integer      :: i, l, n, io
  integer(IXL) :: count0, count1, count_rate
  real(RNP)    :: dx, t_pre, t_sol

  ! Initialization .............................................................

  ! greeting
  write(*,'(/,A)') 'Spectral element solver for 1D Helmholtz equation'

  ! load parameters
  inquire(file='helmholtz_sem_1d.prm', exist=exists)
  if (exists) then
    open(newunit=io, file='helmholtz_sem_1d.prm')
    read(io, nml=problem_parameters)
    read(io, nml=solution_parameters)
    close(io)
  end if

  ! start system clock
  call system_clock(count0, count_rate = count_rate)

  ! standard operators
  call standard_op % New(po)

  ! element operators
  dx = TWO / ne
  allocate(Me(0:po), He(0:po,0:po))
  call GetElementOperators(standard_op, dx, lambda, Me, He)

  ! workspace
  allocate( x(0:po,ne), u(0:po,ne), f(0:po,ne), r(0:po,ne), &
            s(0:po,ne), e(0:po,ne), w(0:po,ne)              )

  ! mesh
  call GetMesh(standard_op, dx, x)

  ! initialize mesh variables

  ! node weights
  call GetNodeWeights(bc, w)

  ! check if problem is singular
  singular = lambda == 0 .and. all(bc == 'N')

  ! right hand side
  call GetRHS(Me, x, f)
  if (singular) then ! project f to nullspace
    f = f - sum(w*f) / sum(w)
  end if

  if (method == 2) then
    call BuildSubOperators(He, HBB, HBI, HII_inv)
    call BuildCondensedSystem(HBB, HBI, HII_inv, bc, f, Ac, fc)
  end if

  ! initial values
  if (init == 0) then ! start from zero
    u = 0
  else ! intial guess, chosen at random from [0,1]
    call random_number(u)
    call MakeContinuous(u)
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
    call CG(He, u, f, w, r_max, i_max)
  case(2) ! static condensation + Gauss elimination
    call SolveCondensedSystem(Ac, fc) ! yields fc = uc on output
    call SolveElementSystems(HBI, HII_inv, bc, fc, f, u)
  end select

  call system_clock(count1)
  t_sol = (count1 - count0) / real(count_rate, RNP)

  ! evaluation .................................................................

  ! number of unknowns
  n = po * ne

  ! consistency error
  s = u_exact(x)
  call HelmholtzResidual(He, s, f, r)
  write(*,'(/,A)') 'consistency error'
  write(*,'(2X,A,ES12.5)') 'c_max', maxval(abs(r))
  write(*,'(2X,A,ES12.5)') 'c_rms', sqrt(sum(w*r*r)/n)

  ! final residual
  call HelmholtzResidual(He, u, f, r)
  write(*,'(/,A)') 'final residual'
  write(*,'(2X,A,ES12.5)') 'r_max', maxval(abs(r))
  write(*,'(2X,A,ES12.5)') 'r_rms', sqrt(sum(w*r*r)/n)

  ! error
  e = u - s
  ! remove constant in singular case
  if (singular) then
    e = e - (maxval(e) + minval(e))/2
  end if
  write(*,'(/,A)') 'error'
  write(*,'(2X,A,ES12.5)') 'e_max', maxval(abs(e))
  write(*,'(2X,A,ES12.5)') 'e_rms', sqrt(sum(w*e*e)/n)

  ! timing
  write(*,'(/,A)') 'runtime'
  write(*,'(2X,A,ES12.5)') 't_pre', t_pre
  write(*,'(2X,A,ES12.5)') 't_sol', t_sol

  ! save results
  open(newunit=io, file='helmholtz_sem_1d.dat')
  write(io,'(A)') '# x, u, s, e, f, r'
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

subroutine GetElementOperators(standard_op, dx, lambda, Me, He)
  class(StandardOperators1D), intent(in) :: standard_op !< standard operators
  real(RNP), intent(in)  :: dx        !< element length
  real(RNP), intent(in)  :: lambda    !< Helmholtz parameter
  real(RNP), intent(out) :: Me(0:)    !< element mass matrix (main diagonal)
  real(RNP), intent(out) :: He(0:,0:) !< element Helmholtz operator

  integer :: i, po

  po = ubound(Me,1)

  associate( Ms => standard_op % w, &
             Ls => standard_op % L  )

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
!> Computes the mesh points

subroutine GetMesh(standard_op, dx, x)
  class(StandardOperators1D), intent(in) :: standard_op !< standard operators
  real(RNP), intent(in)  :: dx      !< element length
  real(RNP), intent(out) :: x(0:,:) !< mesh points x(0:po,1:ne)

  integer   :: l, ne
  real(RNP) :: xm

  ne = ubound(x,2)

  associate(xi => standard_op%x)
    do l = 1, ne
       xm = (l - HALF) * dx - ONE    ! element midpoint
       x(:,l) = xm + HALF * dx * xi  ! transformed GLL points
    end do
  end associate

end subroutine GetMesh

!------------------------------------------------------------------------------
!> Node weights acounting for Dirichlet BC and multiple copies

subroutine GetNodeWeights(bc, w)
  character, intent(in)  :: bc(2)   !< left/right BC
  real(RNP), intent(out) :: w(0:,:) !< node weights

  integer :: po, ne

  po = ubound(w,1)
  ne = ubound(w,2)

  w = 1

  ! element interfaces
  w(po, 1:ne-1) = HALF
  w( 0, 2:ne  ) = HALF

  ! Dirichlet BC
  if (bc(1) == 'D')  w( 0,  1) = 0
  if (bc(2) == 'D')  w(po, ne) = 0

end subroutine GetNodeWeights

!-------------------------------------------------------------------------------
!> Average variable over element boundaries

subroutine MakeContinuous(v)
  real(RNP), intent(inout) :: v(0:,:) !< element contributions

  integer   :: l, po
  real(RNP) :: va

  po = ubound(v,1)

  do l = 1, size(v,2)-1
    va = HALF * (v(po,l) + v(0,l+1))
    v(po,l  ) = va
    v(0 ,l+1) = va
  end do

end subroutine MakeContinuous

!-------------------------------------------------------------------------------
!> Assembly of element contributions

subroutine Assembly(v)
  real(RNP), intent(inout) :: v(0:,:) !< element contributions

  integer   :: l, po
  real(RNP) :: va

  po = ubound(v,1)

  do l = 1, size(v,2)-1
    va = v(po,l) + v(0,l+1)
    v(po,l  ) = va
    v(0 ,l+1) = va
  end do

end subroutine Assembly

!-------------------------------------------------------------------------------
!> Right hand side

subroutine GetRHS(Me, x, f)
  real(RNP), intent(in)  :: Me(0:)  !< element mass matrix (main diagonal)
  real(RNP), intent(in)  :: x(0:,:) !< mesh points
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
  call Assembly(f)

end subroutine GetRHS

!------------------------------------------------------------------------------
!> Residual of a given approximate solution

subroutine HelmholtzResidual(He, u, f, r)
  real(RNP),           intent(in)  :: He(0:,0:) !< element Helmholtz operator
  real(RNP),           intent(in)  :: u(0:,:)   !< approximate solution
  real(RNP), optional, intent(in)  :: f(0:,:)   !< RHS, default: f = 0
  real(RNP),           intent(out) :: r(0:,:)   !< residual, r = f - Au

  ! element contributions
  r = -matmul(He, u)

  ! assembly
  call Assembly(r)

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

subroutine CG(He, u, f, w, r_max, i_max)
  real(RNP), intent(in)    :: He(0:,0:) !< element Helmholtz operator
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
  call HelmholtzResidual(He, u, f, r)

  ! iteration ..................................................................

  delta = sum(w * r * r)
  p = r
  do i = 1, ni
    if (delta <= r_max**2 .or. i > i_max) exit
    call HelmholtzResidual(He, p, r=q)
    alpha = -delta / sum(w * p * q)
    u = u + alpha * p
    if (mod(i,50) == 0) then
      call HelmholtzResidual(He, u, f, r)
    else
      r = r + alpha * q
    end if
    delta_old = delta
    delta = sum(w * r * r)
    beta  = delta / delta_old
    p = r + (delta / delta_old) * p
  end do

end subroutine CG

!------------------------------------------------------------------------------
!> Interior-boundary decomposition of the element Helmholtz operator

subroutine BuildSubOperators(He, HBB, HBI, HII_inv)
  real(RNP), intent(in) :: He(0:,0:) !< element Helmholtz operator
  real(RNP), allocatable, intent(out) :: HBB(:,:) !< boundary-boundary part
  real(RNP), allocatable, intent(out) :: HBI(:,:) !< boundary-interior part
  real(RNP), allocatable, intent(out) :: HII_inv(:,:) !< inverse of HII

  integer :: po

  po = ubound(He,1)
  allocate(HBB(2,2), HBI(2,po-1), HII_inv(po-1,po-1))

  HBB = He([0,po],[0,po])
  HBI = He([0,po],1:po-1)

  HII_inv = Inverse(He(1:po-1,1:po-1))

end subroutine BuildSubOperators

!-------------------------------------------------------------------------------
!> Build the condensed system

subroutine BuildCondensedSystem(HBB, HBI, HII_inv, bc, f, Ac, fc)
  real(RNP), intent(in) :: HBB(:,:)              !< boundary-boundary part
  real(RNP), intent(in) :: HBI(:,:)              !< boundary-interior part
  real(RNP), intent(in) :: HII_inv(:,:)          !< inverse of HII
  character, intent(in) :: bc(2)                 !< boundary conditions
  real(RNP), intent(in) :: f(0:,:)               !< RHS of the whole system
  real(RNP), allocatable, intent(out) :: Ac(:,:) !< condensed system matrix
  real(RNP), allocatable, intent(out) :: fc(:)   !< condensed RHS

  real(RNP) :: a(2,2), b(2,ubound(He,1)-1), df(2)
  integer   :: i, i1, i2, k, po, ne

  ! intialization ..............................................................

  po = ubound(f,1)
  ne = ubound(f,2)

  ! condensed system bounds
  if (bc(1) == 'D') then
    i1 = 1
  else
    i1 = 0
  end if
  if (bc(2) == 'D') then
    i2 = ne-1
  else
    i2 = ne
  end if

  ! condensed element operators ................................................

  b = -matmul(HBI, HII_inv)
  a =  HBB + matmul(b, transpose(HBI))

  ! condensed system matrix ....................................................

  allocate(Ac(i1:i2,3), source=ZERO)

  Ac(i1+1:i2  , 1) = a(2,1)
  Ac(i1  :i2  , 2) = a(1,1) + a(2,2)
  Ac(i1  :i2-1, 3) = a(1,2)

  ! correction for Neumann BC
  if (bc(1) == 'N') Ac(i1,2) = Ac(i1,2) - a(2,2)
  if (bc(2) == 'N') Ac(i2,2) = Ac(i2,2) - a(1,1)

  ! condensed RHS ..............................................................

  allocate(fc(i1:i2))

  ! extract element boundary values
  if (bc(1) == 'N') then
    fc( 0) = f( 0, 1)
  end if
  do i = 1, ne-1
    fc(i) = f(po,i)
  end do
  if (bc(2) == 'N') then
    fc(ne) = f(po,ne)
  end if

  ! add interior contributions
  k = po - 1
  do i = 1, ne
    df = matmul(b, f(1:k,i))
    if (i >  i1) fc(i-1) = fc(i-1) + df(1)
    if (i <= i2) fc(i  ) = fc(i  ) + df(2)
  end do

  ! correction for Dirichlet BC
  if (bc(1) == 'D') fc(i1) = fc(i1) - a(2,1) * u_exact(-ONE)
  if (bc(2) == 'D') fc(i2) = fc(i2) - a(1,2) * u_exact( ONE)

end subroutine BuildCondensedSystem

!------------------------------------------------------------------------------

subroutine SolveCondensedSystem(Ac, fc)
  real(RNP), intent(in)    :: Ac(:,:) !< condensed system matrix
  real(RNP), intent(inout) :: fc(:)   !< condensed RHS (in) / solution (out)

  real(RNP), dimension(size(fc)) :: a, b, c
  logical :: regular, success
  integer :: n

  n = size(fc)

  if (n == 0) then
    return

  else if (n == 1) then
    fc(1) = fc(1) / Ac(1,2)

  else
    a = Ac(:,1)
    b = Ac(:,2)
    c = Ac(:,3)
    regular = abs(b(1)) > abs(c(1)) .or. abs(b(n)) > abs(a(n))
    call TridiagonalSolver(fc, a, b, c, regular, success)

  end if

end subroutine SolveCondensedSystem

!------------------------------------------------------------------------------

subroutine SolveElementSystems(HBI, HII_inv, bc, uc, f, u)
  real(RNP), intent(in)    :: HBI(:,:)     !< boundary-interior element op.
  real(RNP), intent(in)    :: HII_inv(:,:) !< inverse interior element op.
  character, intent(in)    :: bc(2)        !< boundary conditions
  real(RNP), intent(in)    :: uc(:)        !< condensed solution
  real(RNP), intent(in)    :: f(0:,:)      !< full RHS
  real(RNP), intent(inout) :: u(0:,:)      !< full solution

  integer :: i, k, po, ne

  ! bounds
  po = ubound(u,1)
  ne = ubound(u,2)

  ! inject condensed solution
  k = 1
  if (bc(1) == 'N') then
    u( 0, 1) = uc(k)
    k = k + 1
  end if
  do i = 1, ne-1
    u(po, i  ) = uc(k)
    u( 0, i+1) = uc(k)
    k = k + 1
  end do
  if (bc(2) == 'N') then
    u(po,ne) = uc(k)
  end if

  ! solve interior subsystems
  k = po-1
  do i = 1, ne
    u(1:k,i) = matmul(HII_inv, f(1:k,i) - (HBI(1,:)*u(0,i) + HBI(2,:)*u(po,i)))
  end do

end subroutine SolveElementSystems

!===============================================================================

end program Helmholtz_SEM_1D
