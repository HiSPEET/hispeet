!> summary:  Stokes test case of Deville, Kleiser & Montigny-Rannou (1984)
!> author:   Joerg Stiller
!> date:     2018/03/08
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module ISP_Flow_Problem__Stokes_DKM

  use Kind_Parameters,   only: RNP
  use Constants,         only: ZERO, PI
  use Execution_Control
  use Array_Assignments
  use XMPI

  use ISP_Flow_Problem

  implicit none
  private

  public :: FlowProblem_Stokes_DKM

  !-----------------------------------------------------------------------------
  !> Type for defining and handling the Stokes problem

  type, extends(FlowProblem) :: FlowProblem_Stokes_DKM

    real(RNP) :: a = 2.883356_RNP !< wave length

  contains

    procedure :: SetProblem
    procedure :: GetInitialValues
    procedure :: GetBoundaryValues
    procedure :: GetBoundaryTimeDerivative
    procedure :: GetExternalSources
    procedure :: GetExactSolution
    procedure :: GetExactTimeDerivative

  end type FlowProblem_Stokes_DKM

contains

!===============================================================================
! type bound procedures

!-------------------------------------------------------------------------------
!> Initialization

subroutine SetProblem(problem, file, comm)
  class(FlowProblem_Stokes_DKM), intent(inout) :: problem
  character(len=*), optional, intent(in) :: file    !< input file
  type(MPI_Comm),   optional, intent(in) :: comm    !< MPI communicator

  ! local variables ............................................................

  logical   :: stokes = .true.
  real(RNP) :: nu     =  1  ! kinematic viscosity
  real(RNP) :: x0(3)  = -1  ! bounding box: corner nearest to -infinity
  real(RNP) :: x1(3)  =  1  ! bounding box: corner nearest to +infinity
  character, allocatable :: bc(:,:)

  namelist /parameters/ stokes, nu, x0, x1, bc

  logical :: exists
  integer :: prm, rank
  integer :: b1, b2, d

  ! preliminaries ..............................................................

  if (present(comm)) then
    call MPI_Comm_rank(comm, rank)
  else
    rank = 0
  end if

  ! check for input file
  if (rank == 0 .and. present(file)) then

    exists = len_trim(file) > 0
    if (exists) then
      inquire(file=trim(file)//'.prm', exist=exists)
    end if

    if (exists) then
      open(newunit=prm, file=trim(file)//'.prm')
    else
      call Warning('SetProblem', 'Input file "'//trim(file)//'.prm" not found')
    end if

  else
    exists = .false.
  end if

  ! parameters .................................................................

  ! default BC (pressure is ignored)
  allocate(bc(6, problem % nc))
  bc(1:4,:) = 'D'
  bc(5:6,:) = 'P'

  if (rank == 0 .and. exists) then
    read(prm, nml=parameters)
    close(prm)
  end if

  if (present(comm)) then
    call XMPI_Bcast(stokes, 0, comm)
    call XMPI_Bcast(nu    , 0, comm)
    call XMPI_Bcast(x0    , 0, comm)
    call XMPI_Bcast(x1    , 0, comm)
    call XMPI_Bcast(bc    , 0, comm)
  end if

  problem % stokes = stokes
  problem % exact_solution = .true.

  allocate(problem % nu_ref( problem%nc ), source = ZERO)
  problem % nu_ref(1:3) = nu
  problem % x0 = x0
  problem % x1 = x1

  ! pressure BC
  do d = 1, 3
    b1 = 2*d - 1
    b2 = b1  + 1
    if (all(bc(b1:b2, d) == 'P')) then
      bc(b1:b2, 4) = 'P'
    else
      bc(b1:b2, 4) = 'N'
    end if
  end do

  call move_alloc(bc, problem % bc)

end subroutine SetProblem

!-------------------------------------------------------------------------------
!> Provides the initial values u(x,0) for mesh points x

subroutine GetInitialValues(problem, x, u)
  class(FlowProblem_Stokes_DKM), intent(in) :: problem
  real(RNP), intent(in)  :: x(:,:,:,:,:) !< mesh points
  real(RNP), intent(out) :: u(:,:,:,:,:) !< flow variables

  integer :: m, n

  n = size(x(:,:,:,:,1))

  associate(nu => problem%nu_ref(1), a => problem%a)
    call GetVelocity(n, x, ZERO, nu, a, u(:,:,:,:,1:3))
    call GetPressure(n, x, ZERO, nu, a, u(:,:,:,:,4)  )
  end associate

  ! remaining variables get zero
  do m = 5, size(u,5)
    call SetArray(u(:,:,:,:,m), ZERO)
  end do

end subroutine GetInitialValues

!-------------------------------------------------------------------------------
!> Provides the values `ub = u(xb,t)` for points `xb` on boundary `b`

subroutine GetBoundaryValues(problem, b, xb, t, ub)
  class(FlowProblem_Stokes_DKM), intent(in) :: problem
  integer,   intent(in)  :: b             !< boundary ID
  real(RNP), intent(in)  :: xb(:,:,:,:)   !< mesh boundary points
  real(RNP), intent(in)  :: t             !< time
  real(RNP), intent(out) :: ub(:,:,:,:)   !< flow variables

  integer :: m, n

  n = size(xb(:,:,:,1))

  associate(nu => problem%nu_ref(1), a => problem%a)
    call GetVelocity(n, xb, t, nu, a, ub(:,:,:,1:3))
    call GetPressure(n, xb, t, nu, a, ub(:,:,:,4)  )
  end associate

  ! remaining variables get zero
  do m = 5, size(ub,4)
    call SetArray(ub(:,:,:,m), ZERO)
  end do

  ! silence the compiler ;)
  if (b < 0) return

end subroutine GetBoundaryValues

!-------------------------------------------------------------------------------
!> Provides the values `dt_ub = ∂u/∂t(xb,t)` for points `xb` on boundary `b`

subroutine GetBoundaryTimeDerivative(problem, b, xb, t, dt_ub)
  class(FlowProblem_Stokes_DKM), intent(in) :: problem
  integer,   intent(in)  :: b              !< boundary ID
  real(RNP), intent(in)  :: xb(:,:,:,:)    !< boundary points
  real(RNP), intent(in)  :: t              !< time
  real(RNP), intent(out) :: dt_ub(:,:,:,:) !< ∂u/∂t

  integer :: m, n

  n = size(xb(:,:,:,1))

  associate(nu => problem%nu_ref(1), a => problem%a)
    call GetVelocityTimeDerivative(n, xb, t, nu, a, dt_ub(:,:,:,1:3))
    call GetPressureTimeDerivative(n, xb, t, nu, a, dt_ub(:,:,:,4)  )
  end associate

  ! remaining variables get zero
  do m = 5, size(dt_ub,4)
    call SetArray(dt_ub(:,:,:,m), ZERO)
  end do

  ! silence the compiler ;)
  if (b < 0) return

end subroutine GetBoundaryTimeDerivative

!-------------------------------------------------------------------------------
!> Provides the external sources for all variables at points x and time t

subroutine GetExternalSources(problem, x, t, f)
  class(FlowProblem_Stokes_DKM), intent(in) :: problem
  real(RNP), intent(in)  :: x(:,:,:,:,:) !< mesh points
  real(RNP), intent(in)  :: t            !< time
  real(RNP), intent(out) :: f(:,:,:,:,:) !< external sources

  call SetArray(f, ZERO, multi=.true.)

end subroutine GetExternalSources

!---------------------------------------------------------------------------
!> Provides the exact solution u(x,t)

subroutine GetExactSolution(problem, x, t, u)
  class(FlowProblem_Stokes_DKM), intent(in) :: problem
  real(RNP), intent(in)  :: x(:,:,:,:,:) !< mesh points
  real(RNP), intent(in)  :: t            !< time
  real(RNP), intent(out) :: u(:,:,:,:,:) !< solution

  integer :: m, n

  n = size(x(:,:,:,:,1))

  associate(nu => problem%nu_ref(1), a => problem%a)
    call GetVelocity(n, x, t, nu, a, u(:,:,:,:,1:3))
    call GetPressure(n, x, t, nu, a, u(:,:,:,:,4)  )
  end associate

  ! remaining variables get zero
  do m = 5, size(u,5)
    call SetArray(u(:,:,:,:,m), ZERO)
  end do

end subroutine GetExactSolution

!---------------------------------------------------------------------------
!> Provides the time derivative of the exact solution, ∂u/∂t(x,t)

subroutine GetExactTimeDerivative(problem, x, t, dt_u)
  class(FlowProblem_Stokes_DKM), intent(in) :: problem
  real(RNP), intent(in)  :: x(:,:,:,:,:)    !< mesh points
  real(RNP), intent(in)  :: t               !< time
  real(RNP), intent(out) :: dt_u(:,:,:,:,:) !< time derivative

  integer :: m, n

  n = size(x(:,:,:,:,1))

  associate(nu => problem%nu_ref(1), a => problem%a)
    call GetVelocityTimeDerivative(n, x, t, nu, a, dt_u(:,:,:,:,1:3))
    call GetPressureTimeDerivative(n, x, t, nu, a, dt_u(:,:,:,:,4)  )
  end associate

  ! remaining variables get zero
  do m = 5, size(dt_u,5)
    call SetArray(dt_u(:,:,:,:,m), ZERO)
  end do

end subroutine GetExactTimeDerivative

!===============================================================================
! problem-specific procedures

!-------------------------------------------------------------------------------
!> Velocity

subroutine GetVelocity(n, x, t, nu, a, v)
  integer,   intent(in)  :: n       !< number of points
  real(RNP), intent(in)  :: x(n,3)  !< mesh points
  real(RNP), intent(in)  :: t       !< time
  real(RNP), intent(in)  :: nu      !< kinematic viscosity
  real(RNP), intent(in)  :: a       !< wave length
  real(RNP), intent(out) :: v(n,3)  !< velocity at mesh points

  real(RNP) :: lambda, b, c, x1, x2
  integer   :: i

  lambda = nu * (1 + a*a)
  b = cos(a)
  c = exp(-lambda * t)

  do i = 1, n
    x1 = x(i,1)
    x2 = x(i,2)
    v(i,1) = c * sin(x1) * (a * sin(a*x2) - b * sinh(x2))
    v(i,2) = c * cos(x1) * (    cos(a*x2) + b * cosh(x2))
    v(i,3) = ZERO
  end do

end subroutine GetVelocity

!-------------------------------------------------------------------------------
!> Velocity time derivative

subroutine GetVelocityTimeDerivative(n, x, t, nu, a, dt_v)
  integer,   intent(in)  :: n          !< number of points
  real(RNP), intent(in)  :: x(n,3)     !< mesh points
  real(RNP), intent(in)  :: t          !< time
  real(RNP), intent(in)  :: nu         !< kinematic viscosity
  real(RNP), intent(in)  :: a          !< wave length
  real(RNP), intent(out) :: dt_v(n,3)  !< velocity time derivative

  real(RNP) :: lambda, b, c, x1, x2
  integer   :: i

  lambda = nu * (1 + a*a)
  b = cos(a)
  c = -lambda * exp(-lambda * t)

  do i = 1, n
    x1 = x(i,1)
    x2 = x(i,2)
    dt_v(i,1) = c * sin(x1) * (a * sin(a*x2) - b * sinh(x2))
    dt_v(i,2) = c * cos(x1) * (    cos(a*x2) + b * cosh(x2))
    dt_v(i,3) = ZERO
  end do

end subroutine GetVelocityTimeDerivative

!-------------------------------------------------------------------------------
!> Pressure

subroutine GetPressure(n, x, t, nu, a, p)
  integer,   intent(in)  :: n       !< number of points
  real(RNP), intent(in)  :: x(n,3)  !< mesh points
  real(RNP), intent(in)  :: t       !< time
  real(RNP), intent(in)  :: nu      !< kinematic viscosity
  real(RNP), intent(in)  :: a       !< wave length
  real(RNP), intent(out) :: p(n)    !< pressure at mesh points

  real(RNP) :: lambda, c
  integer   :: i

  lambda = nu * (1 + a*a)
  c = lambda * cos(a) * exp(-lambda * t)

  do i = 1, n
    p(i) = c * cos(x(i,1)) * sinh(x(i,2))
  end do

end subroutine GetPressure

!-------------------------------------------------------------------------------
!> Pressure time derivative

subroutine GetPressureTimeDerivative(n, x, t, nu, a, dt_p)
  integer,   intent(in)  :: n       !< number of points
  real(RNP), intent(in)  :: x(n,3)  !< mesh points
  real(RNP), intent(in)  :: t       !< time
  real(RNP), intent(in)  :: nu      !< kinematic viscosity
  real(RNP), intent(in)  :: a       !< wave length
  real(RNP), intent(out) :: dt_p(n) !< pressure time derivative

  real(RNP) :: lambda, c
  integer   :: i

  lambda = nu * (1 + a*a)
  c = -lambda**2 * cos(a) * exp(-lambda * t)

  do i = 1, n
    dt_p(i) = c * cos(x(i,1)) * sinh(x(i,2))
  end do

end subroutine GetPressureTimeDerivative

!===============================================================================

end module ISP_Flow_Problem__Stokes_DKM
