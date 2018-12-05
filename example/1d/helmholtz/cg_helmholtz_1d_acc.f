!> summary:  Spectral element solver for 1D Helmholtz equation with OpenMP/ACC
!> author:   Joerg Stiller
!> date:     2016/10/28
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Spectral element solver for 1D Helmholtz equation with OpenMP/ACC
!
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
!> Discretization is performed using continuous nodal spectral elements based
!> on GLL points. The conjugate gradient (CG) method is used as the solver.
!>
!> Implementation largely follows the approach described in the course "Höhere
!> Numerische Strömungsmechanik" given at TU Dresden in 2016/17. The CG method
!> was adopted from J.R. Shewchuk, An Introduction to the Conjugate Gradient
!> Method Without the Agonizing Pain, Carnegie Mellon University, 1994.
!> For the treatment of the singular problem, i.e. Poisson with Neumann BC, see
!> E.F. Kaasschieter, Preconditioned conjugate gradients for solving singular
!> systems J. Comput. Appl. Math., 1988, 24, 265-275.
!>
!> @note
!> With OpenACC, currently only the compute part is parallelized.
!> @endnote
!>
!===============================================================================

program CG_Helmholtz_1D_Acc
  use Kind_Parameters,  only: RNP, IXL
  use Constants,        only: PI, ZERO, ONE, TWO, HALF
  use Matrix_Operators, only: Inverse
  use Linear_Equations, only: TridiagonalSolver
  use Standard_Operators_1D ! provides the SEM standard operators
  use CG_Utilities_1D
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
  integer   :: i_max  = huge(1) ! maximum number of CG iterations
  real(RNP) :: r_max  = 1E-12   ! maximum CG residual

  namelist /solution_parameters/ po, ne, r_max, i_max

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

  ! auxiliary variables
  logical      :: exists, singular
  integer      :: j, k, n, io
  integer(IXL) :: count0, count1, count_rate
  real(RNP)    :: dx, t_pre, t_sol
  real(RNP)    :: sw, sf, df

  ! Initialization .............................................................

  ! greeting
  write(*,'(/,A)') 'Accelerated CG-SEM for 1D Helmholtz equation'

  ! load parameters
  inquire(file='cg_helmholtz_1d_acc.prm', exist=exists)
  if (exists) then
    open(newunit=io, file='cg_helmholtz_1d_acc.prm')
    read(io, nml=problem_parameters)
    read(io, nml=solution_parameters)
    close(io)
  end if

  ! start system clock
  call system_clock(count0, count_rate = count_rate)

  ! workspace
  allocate( x(0:po,ne), u(0:po,ne), f(0:po,ne), r(0:po,ne), &
            s(0:po,ne), e(0:po,ne), w(0:po,ne)              )

  ! standard operators
  call standard_op % New(po)

  ! check if problem is singular
  singular = lambda == 0 .and. all(bc == 'N')

  !$omp parallel default(shared)

  ! mesh
  call GetMeshPoints(standard_op, -ONE, ONE, dx, x)
  call GetPointWeights(bc, w)

  ! element operators
  !$omp single
  allocate(Me(0:po), He(0:po,0:po))
  call GetElementOperators(standard_op, dx, lambda, Me, He)
  !$omp end single

  ! right hand side
  call GetRHS(Me, x, f)
  if (singular) then ! project f to nullspace
    sw = 0
    sf = 0
    !$omp barrier
    !$omp do reduction(+:sw,sf)
    do k = 1, ne
    do j = 0, po
      sw = sw + w(j,k)
      sf = sf + w(j,k) * f(j,k)
    end do
    end do
    df = -sf / sw
    !$omp do
    do k = 1, ne
    do j = 0, po
      f(j,k) = f(j,k) + df
    end do
    end do
  end if

  ! initial values
  if (init == 0) then ! start from zero
    do k = 1, ne
    do j = 0, po
      u(j,k) = 0
    end do
    end do
  else ! intial guess, chosen at random from [0,1]
    call random_number(u)
    call MakeContinuous(bc, u)
  end if

  ! inject Dirichlet BC
  !$omp single
  if (bc(1) == 'D')  u( 0,  1) = u_exact(-ONE)
  if (bc(2) == 'D')  u(po, ne) = u_exact( ONE)
  !$omp end single

  !$omp end parallel

  ! read system clock
  call system_clock(count1)

  ! time consumed with preprocessing
  t_pre = (count1 - count0) / real(count_rate, RNP)

  ! solution ...................................................................

  call system_clock(count0)

  !$omp parallel default(shared)
  !$acc data copy(u) copyin(He,f,w)
  call CG(He, u, f, w, r_max, i_max)
  !$acc end data
  !$omp end parallel

  call system_clock(count1)
  t_sol = (count1 - count0) / real(count_rate, RNP)

  ! evaluation .................................................................

  !$omp parallel default(shared)

  ! number of unknowns
  n = po * ne

  ! consistency error
  !$omp workshare
  s = u_exact(x)
  !$omp end workshare
  !$acc data copy(r) copyin(He,s,f)
  call HelmholtzResidual(He, s, f, r)
  !$acc end data
  !$omp master
  write(*,'(/,A)') 'consistency error'
  write(*,'(2X,A,ES12.5)') 'c_max', maxval(abs(r))
  write(*,'(2X,A,ES12.5)') 'c_rms', sqrt(sum(w*r*r)/n)
  !$omp end master
  !$omp barrier

  ! final residual
  !$acc data copy(r) copyin(He,u,f)
  call HelmholtzResidual(He, u, f, r)
  !$acc end data
  !$omp barrier
  !$omp master
  write(*,'(/,A)') 'final residual'
  write(*,'(2X,A,ES12.5)') 'r_max', maxval(abs(r))
  write(*,'(2X,A,ES12.5)') 'r_rms', sqrt(sum(w*r*r)/n)
  !$omp end master

  ! error
  !$omp workshare
  e = u - s
  !$omp end workshare
  ! remove constant in singular case
  if (singular) then
    !$omp workshare
    e = e - (maxval(e) + minval(e))/2
    !$omp end workshare
  end if
  !$omp barrier
  !$omp master
  write(*,'(/,A)') 'error'
  write(*,'(2X,A,ES12.5)') 'e_max', maxval(abs(e))
  write(*,'(2X,A,ES12.5)') 'e_rms', sqrt(sum(w*e*e)/n)
  !$omp end master

  !$omp end parallel

  ! timing
  write(*,'(/,A)') 'runtime'
  write(*,'(2X,A,ES12.5)') 't_pre', t_pre
  write(*,'(2X,A,ES12.5)') 't_sol', t_sol

  ! save results
  open(newunit=io, file='helmholtz_sem_1d.dat')
  write(io,'(A)') '# x, u, s, e, f, r'
  do k = 1, ne
  do j = 0, po
    write(io,'(10(ES17.10,1X))') x(j,k), u(j,k), s(j,k), e(j,k), f(j,k), r(j,k), w(j,k)
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
!> Weigthed dot product

subroutine WeightedDotProduct(w, u, v, r)
  real(RNP), intent(in)  :: w(0:,:) !< node weights                     !shared!
  real(RNP), intent(in)  :: u(0:,:) !< left factor                      !shared!
  real(RNP), intent(in)  :: v(0:,:) !< right factor                     !shared!
  real(RNP), intent(out) :: r       !< weighted dot product of u and v

  integer :: j, k, po, ne
  real(RNP), save :: s

  po = ubound(w,1)
  ne = ubound(w,2)

  s = 0

  !$omp barrier
  !$omp do reduction(+:s)
  !$acc parallel loop collapse(2) reduction(+:s) present(u,v,w)
  do k = 1, ne
  do j = 0, po
    s = s + w(j,k) * u(j,k) * v(j,k)
  end do
  end do
  !$acc end parallel

  r = s

end subroutine WeightedDotProduct

!-------------------------------------------------------------------------------
!> Right hand side

subroutine GetRHS(Me, x, f)
  real(RNP), intent(in)  :: Me(0:)  !< diagonal element mass matrix     !shared!
  real(RNP), intent(in)  :: x(0:,:) !< mesh points                      !shared!
  real(RNP), intent(out) :: f(0:,:) !< RHS                              !shared!

  integer :: k, po, ne

  po = ubound(x,1)
  ne = ubound(x,2)

  ! projection of the source term
  !$omp do
  do k = 1, ne
     f(:,k) = Me * (lambda*u_exact(x(:,k)) - ddu_exact(x(:,k)))
  end do

  ! Neumann BC
  !$omp single
  if (bc(1) == 'N')  f( 0,  1) = f( 0,  1) - du_exact(-ONE)
  if (bc(2) == 'N')  f(po, ne) = f(po, ne) + du_exact( ONE)
  !$omp end single

  ! assemble element contributions
  !$acc data copy(f)
  call Assembly(bc, f)
  !$acc end data

end subroutine GetRHS

!------------------------------------------------------------------------------
!> Residual of a given approximate solution -- to be placed in parallel region

subroutine HelmholtzResidual(He, u, f, r)
  real(RNP), intent(in)  :: He(0:,0:) !< element Helmholtz operator     !shared!
  real(RNP), intent(in)  :: u(0:,:)   !< approximate solution           !shared!
  real(RNP), intent(in)  :: f(0:,:)   !< RHS, default: f = 0            !shared!
  real(RNP), intent(out) :: r(0:,:)   !< residual, r = f - Au           !shared!
  optional :: f

  integer :: i, j, k, po, ne
  real(RNP) :: s

  ! dimensions
  po = ubound(u,1)
  ne = ubound(u,2)

  ! element contributions
  !$omp do
  !$acc parallel loop collapse(2) present(He,u,r)
  do k = 1, ne
    do i = 0, po
      s = 0
      !$acc loop seq
      do j = 0, po
        s = s + He(i,j) * u(j,k)
      end do
      r(i,k) = -s
    end do
  end do
  !$acc end parallel

  ! assembly
  call Assembly(bc, r)

  ! add RHS contribution
  if (present(f)) then
    !$omp do
    !$acc parallel loop collapse(2) present(r,f)
    do k = 1, ne
    do j = 0, po
      r(j,k) = r(j,k) + f(j,k)
    end do
    end do
    !$acc end parallel
  end if

  ! nullify residual at Dirichlet points
  !$omp single
  if (bc(1) == 'D') then
    r( 0,  1) = 0
    !$acc update device(r(0:0,1:1))
  end if
  if (bc(2) == 'D') then
    r(po, ne) = 0
    !$acc update device(r(po:po,ne:ne))
  end if
  !$omp end single

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

  ! shared variables
  real(RNP), dimension(:,:), allocatable, save :: r, p, q
  logical, save :: finished

  ! private variables
  real(RNP) :: alpha, beta, delta, delta_old, pq
  integer   :: i, j, k, po, ne, ni

  ! initialization .............................................................

  ! dimensions
  po = ubound(u,1)
  ne = ubound(u,2)

  ! number of iterations, with safety factor of 10
  ni = 10 * ne * po*po

  ! workspace
  !$omp single
  allocate(r(0:po,1:ne), p(0:po,1:ne), q(0:po,1:ne))
  !$omp end single

  !$acc data create(r, p, q) present(u)

  ! initial residual
  call HelmholtzResidual(He, u, f, r)

  ! iteration ..................................................................

  call WeightedDotProduct(w, r, r, delta)

  !$omp do
  !$acc parallel loop collapse(2)
  do k = 1, ne
  do j = 0, po
    p(j,k) = r(j,k)
  end do
  end do
  !$acc end parallel

  do i = 1, ni

    ! termination check  . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    !$omp single
    finished = delta <= r_max**2 .or. i > i_max
    !$omp end single

    if (finished) exit

    ! next iteration . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    call HelmholtzResidual(He, p, r=q)

    call WeightedDotProduct(w, p, q, pq)
    alpha = -delta / pq

    !$omp do
    !$acc parallel loop collapse(2)
    do k = 1, ne
    do j = 0, po
      u(j,k) = u(j,k) + alpha * p(j,k)
    end do
    end do
    !$acc end parallel

    if (mod(i,50) == 0) then
      call HelmholtzResidual(He, u, f, r)
    else
      !$omp do
      !$acc parallel loop collapse(2)
      do k = 1, ne
      do j = 0, po
        r(j,k) = r(j,k) + alpha * q(j,k)
      end do
      end do
      !$acc end parallel
    end if

    delta_old = delta

    !$omp barrier
    call WeightedDotProduct(w, r, r, delta)

    beta = delta / delta_old

    !$omp do
    !$acc parallel loop collapse(2)
    do k = 1, ne
    do j = 0, po
      p(j,k) = r(j,k) + beta * p(j,k)
    end do
    end do
    !$acc end parallel

  end do

  !$acc end data

end subroutine CG

!===============================================================================

end program CG_Helmholtz_1D_Acc
