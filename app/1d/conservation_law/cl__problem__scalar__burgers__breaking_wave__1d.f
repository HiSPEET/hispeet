!> summary:  1D Burgers breaking wave problem
!> author:   Joerg Stiller
!> date:     2019/11/24
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module CL__Problem__Scalar__Burgers__Breaking_Wave__1D

  use Kind_Parameters,   only: RNP
  use Constants,         only: HALF, ONE, TWO, PI
  use Execution_Control
  use CL__Problem__Scalar__Burgers__1D

  implicit nONE
  private

  public :: CL_Problem_Scalar_Burgers_BreakingWave_1D

  !-----------------------------------------------------------------------------
  !> Type 1D Burgers breaking wave problem

  type, extends(CL_Problem_Scalar_Burgers_1D) :: &
    CL_Problem_Scalar_Burgers_BreakingWave_1D
  contains
    procedure :: SetProblem
    procedure :: InitialValues
    procedure :: ExactSolution
  end type CL_Problem_Scalar_Burgers_BreakingWave_1D

contains

  !-----------------------------------------------------------------------------
  !> Initialization of the Burgers problem

  subroutine SetProblem(problem, file)
    class(CL_Problem_Scalar_Burgers_BreakingWave_1D), intent(inout) :: problem
    character(len=*), optional, intent(in) :: file !< input file (*.prm)

    real(RNP) :: xb1    = -1  ! position of first boundary  (left)
    real(RNP) :: xb2    =  1  ! position of second boundary (right)
    character :: bc(2)  = 'P' ! boundary conditions
    real(RNP) :: nu_c   =  0  ! constant regular  viscosity

    namelist /burgers_breaking_wave_prm/ nu_c

    logical :: exists, opened
    integer :: prm

    ! check for input file .....................................................

    if (present(file)) then

      inquire(file=trim(file)//'.prm', exist=exists, opened=opened, number=prm)

      if (exists .and. .not. opened) then
        open(newunit=prm, file=trim(file)//'.prm', action='READ')
      else if (opened) then
        rewind(prm)
      else
        call Error('SetProblem',                                      &
                     'Input file "'//trim(file)//'.prm" not found, '//  &
                     'using defaults',                                  &
                     'CL_Problem_Scalar_Burgers_BreakingWave_1D')
      end if

    else
      exists = .false.
    end if

    ! read parameters ..........................................................

    if (exists) then
      read(prm, nml=burgers_breaking_wave_prm)
      if (.not. opened) close(prm)
    end if

    ! set parameters ...........................................................

    problem % nc   = 1
    problem % xb1  = xb1
    problem % xb2  = xb2
    problem % bc   = reshape(bc, shape = [2,1])
    problem % nu_c = nu_c

  end subroutine SetProblem

  !-----------------------------------------------------------------------------
  !> Provides the initial values u(x,0) for mesh points x

  function InitialValues(problem) result(u)
    class(CL_Problem_Scalar_Burgers_BreakingWave_1D), intent(in) :: problem
    real(RNP) :: u(0:problem%eop%po, problem%ne, problem%nc)

    u(:,:,1) = ONE + HALF * sin(PI * problem % x + TWO - PI)

  end function InitialValues

  !-----------------------------------------------------------------------------
  !> Provides the exact solution at mesh points x for Parameter T

  subroutine ExactSolution(problem, x, t, N, ou) ! t = t_end for example
    class(CL_Problem_Scalar_Burgers_BreakingWave_1D), intent(in) :: problem
    real(RNP), contiguous, intent(in) :: x(0:,:)
    real(RNP), intent(in)             :: t                   ! time
    integer,   intent(in), optional   :: N                   ! series expansion degree
    integer   :: ou              ! series expansion degree
    real(RNP) :: u(0:problem%eop%po, problem%ne, problem%nc) ! u(x,t)

  ! Note: If N are not given, the old values will be used.
  ! If not specified at all the following defaults apply:
    real(RNP) ::  nu  = ONE ! diffusivity
    real(RNP) ::  nu0 = ONE
    integer   ::  N0  = 100

  ! local data ..................................................................
    real(RNP), allocatable, save ::  a(:)
    real(RNP) ::  c, s, q, e
    integer   ::  m, i, k

  ! for use with NR Bessel function
    double precision :: rc, rm, ri, rip, rk, rkp

    nu = problem % nu_c

  ! initialize coefficients, if necessary .......................................
    if (allocated(a) .and. abs(nu-nu0) > epsilon(nu)) deallocate(a)

    if (nu0 < epsilon(nu)) then
       u = 0
       return
    end if
    if (present(N)) then
       if (allocated(a).and. N /= N0) deallocate(a)
    end if
    if (.not. allocated(a)) then
       allocate(a(0:N0))
       c = 1 / (2*pi*nu0)
       do m = 0, N0
          rc = c
          rm = m
          call Bessik(rc, rm, ri, rk, rip, rkp)
          a(m) = (-1)**m * ri
       end do
    end if

  ! compute solution ............................................................
    do i = 1, problem % ne 
    do k = 0, problem % eop % po
       c =-pi**2 * t * nu0
       s = 0
       q = 0
       do m = 1, N0
          e = exp(c * m**2)
          s = s + m * a(m) * e * sin(m*pi*x(i,k))
          q = q + a(m) * e * cos(m*pi*x(i,k))
       end do
       u(i,k,:) = 4*pi*nu0*s / (a(0) + 2*q)
    end do
    end do

  ! write solution in file ......................................................
    open(newunit=ou, file='conservation_law_reference.dat')
    write(ou,'(A)') '# x, u'
    do k = 1, problem % ne
    do i = 0, problem % eop % po
      write(ou,'(10(ES17.10,1X))') problem % x(i,k), u(i,k,:)
    end do
    end do

  end subroutine ExactSolution












! procedures adopted from Numerical Recipes ===================================

subroutine Bessik(x, xnu, ri, rk, rip, rkp)
!  use Constants, only: RNP, pi
  real(RNP), intent(in)  :: x, xnu
  real(RNP), intent(out) :: ri, rip, rk, rkp

  real(RNP), parameter :: EPS = 1e-10, FPMIN = 1e-30, XMIN = 2
  integer  , parameter :: MAXIT = 10000

  real(RNP) :: a, a1, b, c, d, del, del1, delh, dels, e, f, fact,    &
   fact2, ff, gam1, gam2, gammi, gampl, h, p, pimu, q, q1, q2, qnew, &
   ril, ril1, rimu, rip1, ripl, ritemp, rk1, rkmu, rkmup, rktemp, s, &
   sum, sum1, x2, xi, xi2, xmu, xmu2
  integer :: i, l, nl

  if (x <= 0 .or. xnu < 0) stop "bad arguments in Bessik"
  nl   = nint(xnu)
  xmu  = xnu - nl
  xmu2 = xmu * xmu
  xi   = 1 / x
  xi2  = 2 * xi
  h    = xnu * xi
  if (h < FPMIN) h = FPMIN
  b = xi2 * xnu
  d = 0
  c = h
  do i = 1, MAXIT
     b = b + xi2
     d = 1 / (b + d)
     c = b + 1 / c
     del = c * d
     h = del * h
     if (abs(del - 1) < EPS) exit
  end do
  if (i > MAXIT) stop "x too large in Bessik; try asymptotic expansion"
  ril = FPMIN
  ripl = h * ril
  ril1 = ril
  rip1 = ripl
  fact = xnu * xi
  do l = nl, 1, -1
     ritemp = fact * ril + ripl
     fact = fact - xi
     ripl = fact * ritemp + ril
     ril = ritemp
  end do
  f = ripl / ril
  if (x < XMIN) then
     x2 = x/2
     pimu = pi * xmu
     if (abs(pimu) < EPS) then
        fact = 1
     else
        fact = pimu / sin(pimu)
     end if
     d = -log(x2)
     e = xmu * d
     if (abs(e) < EPS) then
        fact2 = 1
     else
        fact2 = sinh(e) / e
     end if
     call Beschb(xmu, gam1, gam2, gampl, gammi)
     ff = fact * (gam1 * cosh(e) + gam2 * fact2 * d)
     sum = ff
     e = exp(e)
     p = e / (2 * gampl)
     q = 1 / (2 * e * gammi)
     c = 1
     d = x2 * x2
     sum1 = p
     do i = 1, MAXIT
        ff = (i * ff + p + q) / (i * i - xmu2)
        c = c * d / i
        p = p / (i - xmu)
        q = q / (i + xmu)
        del = c * ff
        sum = sum + del
        del1 = c * (p - i * ff)
        sum1 = sum1 + del1
        if (abs(del) < abs(sum) * EPS) exit
     end do
     if (i > MAXIT) stop "Bessik series failed to converge"
     rkmu = sum
     rk1 = sum1 * xi2
  else
     b = 2 * (1 + x)
     d = 1 / b
     delh = d
     h = delh
     q1 = 0
     q2 = 1
     a1 = 0.25_RNP - xmu2
     c  = a1
     q  = c
     a  = -a1
     s  = 1 + q * delh
     do i = 2, MAXIT
        a = a - 2 * (i - 1)
        c = -a * c / i
        qnew = (q1 - b * q2) / a
        q1 = q2
        q2 = qnew
        q = q + c * qnew
        b = b + 2
        d = 1 / (b + a * d)
        delh = (b * d - 1) * delh
        h = h + delh
        dels = q * delh
        s = s + dels
        if (abs(dels/s) < EPS) exit
     end do
     if (i > MAXIT) stop "Bessik series failed to converge in cf2"
     h = a1 * h
     rkmu = sqrt (pi / (2 * x) ) * exp(-x) / s
     rk1 = rkmu * (xmu + x + 0.5_RNP - h) * xi
  end if
  rkmup = xmu * xi * rkmu - rk1
  rimu = xi / (f * rkmu - rkmup)
  ri = (rimu * ril1) / ril
  rip = (rimu * rip1) / ril
  do i = 1, nl
     rktemp = (xmu + i) * xi2 * rk1 + rkmu
     rkmu = rk1
     rk1 = rktemp
  end do
  rk = rkmu
  rkp = xnu * xi * rkmu - rk1
end subroutine Bessik

!------------------------------------------------------------------------------

subroutine Beschb(x, gam1, gam2, gampl, gammi)
!  use Constants, only: RNP
  real(RNP), intent(in)  :: x
  real(RNP), intent(out) :: gam1, gam2, gampl, gammi
  integer, parameter :: NUSE1 = 5, NUSE2 = 5
  real(RNP) :: xx, ONE = 1
  real(RNP) :: c1(7) = (/ -1.142022680371168d+0, &
                           6.5165112670737d-3,   &
                           3.087090173086d-4,    &
                          -3.4706269649d-6,      &
                           6.9437664d-9,         &
                           3.67795d-11,          &
                          -1.356d-13            /)
  real(RNP) :: c2(8) = (/  1.843740587300905d0,  &
                          -7.68528408447867d-2,  &
                           1.2719271366546d-3,   &
                          -4.9717367042d-6,      &
                          -3.31261198d-8,        &
                           2.423096d-10,         &
                          -1.702d-13,            &
                          -1.49d-15             /)
  xx = 8*x*x - 1
  gam1 = Chebev ( - ONE, ONE, c1, NUSE1, xx)
  gam2 = Chebev ( - ONE, ONE, c2, NUSE2, xx)
  gampl = gam2 - x * gam1
  gammi = gam2 + x * gam1
end subroutine Beschb

!------------------------------------------------------------------------------

function Chebev(a, b, c, m, x)
!  use Constants, only: RNP
  integer,   intent(in) :: m
  real(RNP), intent(in) :: a, b, c(m), x
  real(RNP) :: Chebev

! local declarations:
  real(RNP) :: d, dd, sv, y, y2
  integer   :: j

! body:
  if ((x-a)*(x-b) > 0) stop "x not in range in Chebev"
  d  = 0
  dd = 0
  y  = (2*x - a - b) / (b - a)
  y2 = 2*y
  do j = m, 2, -1
   sv = d
   d  = y2*d - dd + c(j)
   dd = sv
  end do
  Chebev = y*d - dd + c(1)/2
end function Chebev







  !=============================================================================

end module CL__Problem__Scalar__Burgers__Breaking_Wave__1D
