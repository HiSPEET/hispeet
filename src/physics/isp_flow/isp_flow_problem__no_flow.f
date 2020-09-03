!> summary:  No-flow test case
!> author:   Joerg Stiller
!> date:     2020/08/22
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module ISP_Flow_Problem__No_Flow

  use Kind_Parameters,   only: RNP
  use Constants,         only: ZERO, ONE, HALF, PI
  use Execution_Control
  use Array_Assignments
  use XMPI

  use ISP_Flow_Problem

  implicit none
  private

  public :: FlowProblem_No_Flow

  !-----------------------------------------------------------------------------
  !> Type for defining and handling the problem

  type, extends(FlowProblem) :: FlowProblem_No_Flow

    integer   :: init     !< 0 zero, 1 harmonic, 2 random, 3 boundary
    real(RNP) :: alpha(3) !< perturbation amplitude of components 1:3
    integer   :: k(3)     !< harmonic wave number in directions 1:3
    integer   :: m        !< boundary perturbation exponent
    integer   :: cb(6)    !< boundary perturbation coefficients

  contains

    procedure :: SetProblem
    procedure :: GetInitialValues
    procedure :: GetBoundaryValues
    procedure :: GetBoundaryTimeDerivative
    procedure :: GetExternalSources
    procedure :: GetExactSolution
    procedure :: GetExactTimeDerivative

  end type FlowProblem_No_Flow

contains

!===============================================================================
! type bound procedures

!-------------------------------------------------------------------------------
!> Initialization

subroutine SetProblem(problem, file, comm)
  class(FlowProblem_No_Flow), intent(inout) :: problem
  character(len=*), optional, intent(in) :: file  !< input file
  type(MPI_Comm),   optional, intent(in) :: comm  !< MPI communicator

  ! local variables ............................................................

  logical   :: stokes   = .false.   ! switch from Navier-Stokes to Stokes
  integer   :: init     = 1         ! 0 zero, 1 harmonic, 2 random, 3 boundary
  real(RNP) :: alpha(3) = 1         ! perturbation amplitude
  integer   :: k(3)     = 1         ! harmonic wave number in directions 1:3
  integer   :: m        = 2         ! boundary perturbation exponent
  integer   :: cb(6)    = 1         ! boundary coefficients: ±1 out/in, 0 none
  real(RNP) :: nu       = 1         ! kinematic viscosity
  character, allocatable :: bc(:,:) ! boundary conditions

  namelist /parameters/ stokes, init, alpha, k, m, cb, nu, bc

  logical :: exists
  integer :: prm, rank

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

  ! default boundary conditions and jumps (pressure is ignored)
  allocate(bc(6, problem % nc), source = 'P')

  if (rank == 0 .and. exists) then
    read(prm, nml=parameters)
    close(prm)
  end if

  if (present(comm)) then
    call XMPI_Bcast(stokes, 0, comm)
    call XMPI_Bcast(init  , 0, comm)
    call XMPI_Bcast(alpha , 0, comm)
    call XMPI_Bcast(k     , 0, comm)
    call XMPI_Bcast(m     , 0, comm)
    call XMPI_Bcast(cb    , 0, comm)
    call XMPI_Bcast(nu    , 0, comm)
    call XMPI_Bcast(bc    , 0, comm)
  end if

  allocate(problem % nu_ref( problem%nc ), source = ZERO)

  problem % stokes          =  stokes
  problem % exact_solution  =  .true.

  problem % init            =  init
  problem % alpha           =  alpha
  problem % k               =  k
  problem % m               =  m
  problem % cb              =  cb

  problem % v_ref           =  ONE
  problem % nu_ref(1:3)     =  nu
  problem % x0              = -HALF
  problem % x1              =  HALF

  call move_alloc(bc, problem % bc)
  call problem % GeneratePressureBC()

end subroutine SetProblem

!-------------------------------------------------------------------------------
!> Provides the initial values u(x,0) for mesh points x

subroutine GetInitialValues(problem, x, u)
  class(FlowProblem_No_Flow), intent(in) :: problem
  real(RNP), intent(in)  :: x(:,:,:,:,:) !< mesh points
  real(RNP), intent(out) :: u(:,:,:,:,:) !< flow variables

  real(RNP) :: kappa(3)
  integer   :: m, n

  n = size(x(:,:,:,:,1))

  select case(problem % init)
  case(1) ! harmonic wave
    kappa = 2 * PI * problem % k
    call HarmonicPerturbation(n, problem % alpha, kappa, x, u(:,:,:,:,1:3))
  case(2) ! random
    call RandomPerturbation(n, problem % alpha, u(:,:,:,:,1:3))
  case(3) ! jump
    call BoundaryPerturbation(n, problem % alpha, problem % m, problem % cb, &
                              x, u(:,:,:,:,1:3))
  case default ! zero (init == 0)
    call SetArray(u(:,:,:,:,1:3), ZERO, multi = .true.)
  end select

  ! remaining variables get zero
  do m = 4, size(u,5)
    call SetArray(u(:,:,:,:,m), ZERO)
  end do

end subroutine GetInitialValues

!-------------------------------------------------------------------------------
!> Provides the values `ub = u(xb,t)` for points `xb` on boundary `b`

subroutine GetBoundaryValues(problem, b, xb, t, ub)
  class(FlowProblem_No_Flow), intent(in) :: problem
  integer,   intent(in)  :: b             !< boundary ID
  real(RNP), intent(in)  :: xb(:,:,:,:)   !< mesh boundary points
  real(RNP), intent(in)  :: t             !< time
  real(RNP), intent(out) :: ub(:,:,:,:)   !< flow variables

  call SetArray(ub(:,:,:,:), ZERO, multi = .true.)

  ! silence the compiler ;)
  if (b < 0 .or. size(xb) < 0 .or. t < 0) return

end subroutine GetBoundaryValues

!-------------------------------------------------------------------------------
!> Provides the values `dt_ub = ∂u/∂t(xb,t)` for points `xb` on boundary `b`

subroutine GetBoundaryTimeDerivative(problem, b, xb, t, dt_ub)
  class(FlowProblem_No_Flow), intent(in)  :: problem
  integer,   intent(in)  :: b              !< boundary ID
  real(RNP), intent(in)  :: xb(:,:,:,:)    !< boundary points
  real(RNP), intent(in)  :: t              !< time
  real(RNP), intent(out) :: dt_ub(:,:,:,:) !< ∂u/∂t

  call SetArray(dt_ub, ZERO, multi=.true.)

  ! silence the compiler ;)
  if (b < 0 .or. size(xb) < 0 .or. t < 0) return

end subroutine GetBoundaryTimeDerivative

!-------------------------------------------------------------------------------
!> Provides the external sources for all variables at points x and time t

subroutine GetExternalSources(problem, x, t, f)
  class(FlowProblem_No_Flow), intent(in) :: problem
  real(RNP), intent(in)  :: x(:,:,:,:,:) !< mesh points
  real(RNP), intent(in)  :: t            !< time
  real(RNP), intent(out) :: f(:,:,:,:,:) !< external sources

  call SetArray(f, ZERO, multi=.true.)

  ! silence the compiler ;)
  if (size(x) < 0 .or. t < 0) return

end subroutine GetExternalSources

!---------------------------------------------------------------------------
!> Provides the exact solution u(x,t), or zero if not available

subroutine GetExactSolution(problem, x, t, u)
  class(FlowProblem_No_Flow), intent(in) :: problem
  real(RNP), intent(in)  :: x(:,:,:,:,:) !< mesh points
  real(RNP), intent(in)  :: t            !< time
  real(RNP), intent(out) :: u(:,:,:,:,:) !< solution

  call SetArray(u, ZERO, multi = .true.)

  ! silence the compiler ;)
  if (size(x) < 0 .or. t < 0) return

end subroutine GetExactSolution

!---------------------------------------------------------------------------
!> Provides the time derivative of the exact solution, ∂u/∂t(x,t)

subroutine GetExactTimeDerivative(problem, x, t, dt_u)
  class(FlowProblem_No_Flow), intent(in) :: problem
  real(RNP), intent(in)  :: x(:,:,:,:,:)    !< mesh points
  real(RNP), intent(in)  :: t               !< time
  real(RNP), intent(out) :: dt_u(:,:,:,:,:) !< time derivative

  call SetArray(dt_u(:,:,:,:,:), ZERO, multi = .true.)

  ! silence the compiler ;)
  if (size(x) < 0 .or. t < 0) return

end subroutine GetExactTimeDerivative

!===============================================================================
! problem-specific procedures

!-------------------------------------------------------------------------------
!> Solenoidal harmonic perturbation
!>
!>  Given
!>
!>        f(x,y,z) = -1/|κ| sin(κ₁x) sin(κ₂y) sin(κ₃z)
!>
!>  where |κ| is the length of the wave vector, the perturbation is defined as
!>
!>        v₁ = α₁ (∂f/∂z - ∂f/∂y)
!>        v₂ = α₂ (∂f/∂x - ∂f/∂z)
!>        v₃ = α₃ (∂f/∂y - ∂f/∂x)
!>
!>  with
!>
!>        ∂f/∂x = κ₁/|κ| cos(κ₁x) sin(κ₂y) sin(κ₃z)
!>        ∂f/∂y = κ₂/|κ| sin(κ₁x) cos(κ₂y) sin(κ₃z)
!>        ∂f/∂z = κ₃/|κ| sin(κ₁x) sin(κ₂y) cos(κ₃z)

subroutine HarmonicPerturbation(n, alpha, kappa, x, v)
  integer,   intent(in)  :: n        !< number of points
  real(RNP), intent(in)  :: alpha(3) !< perturbation amplitude
  real(RNP), intent(in)  :: kappa(3) !< wave vector
  real(RNP), intent(in)  :: x(n,3)   !< mesh points
  real(RNP), intent(out) :: v(n,3)   !< velocity at mesh points

  real(RNP) :: c(3), k_x, k_y, k_z, f_x, f_y, f_z
  integer   :: i

  if (any(kappa /= ZERO)) then
    c = alpha * kappa / sqrt(sum(kappa**2))
  else
    c = alpha * kappa
  end if

  !$omp do
  do i = 1, n
    k_x = kappa(1) * x(i,1)
    k_y = kappa(2) * x(i,2)
    k_z = kappa(3) * x(i,3)
    f_x = c(1) * cos(k_x) * sin(k_y) * sin(k_z)
    f_y = c(2) * sin(k_x) * cos(k_y) * sin(k_z)
    f_z = c(3) * sin(k_x) * sin(k_y) * cos(k_z)
    v(i,1) = f_z - f_y
    v(i,2) = f_x - f_z
    v(i,3) = f_y - f_x
  end do

end subroutine HarmonicPerturbation

!-------------------------------------------------------------------------------
!> Discontinuous random perturbation

subroutine RandomPerturbation(n, alpha, v)
  integer,   intent(in)  :: n        !< number of points
  real(RNP), intent(in)  :: alpha(3) !< perturbation amplitude
  real(RNP), intent(out) :: v(n,3)   !< velocity at mesh points

  integer   :: i

  call SetArray(v, ZERO, multi = .true.)
  call random_number(v)

  !$omp do
  do i = 1, n
    v(i,1) = alpha(1) * (2*v(i,1) - 1)
    v(i,2) = alpha(2) * (2*v(i,2) - 1)
    v(i,3) = alpha(3) * (2*v(i,3) - 1)
  end do

end subroutine RandomPerturbation

!-------------------------------------------------------------------------------
!> Boundary perturbation

subroutine BoundaryPerturbation(n, alpha, m, cb, x, v)
  integer,   intent(in)  :: n        !< number of points
  real(RNP), intent(in)  :: alpha(3) !< perturbation amplitude
  integer,   intent(in)  :: m        !< perturbation exponent
  integer,   intent(in)  :: cb(6)    !< boundary coefficients: ±1 out/in, 0 none
  real(RNP), intent(in)  :: x(n,3)   !< mesh points
  real(RNP), intent(out) :: v(n,3)   !< velocity at mesh points

  integer :: i

  call SetArray(v, ZERO, multi = .true.)

  !$omp do
  do i = 1, n
    v(i,:) = alpha * ( - cb(1) * (x(i,1) + HALF)**m  &
                       + cb(2) * (x(i,1) - HALF)**m  &
                       - cb(3) * (x(i,2) + HALF)**m  &
                       + cb(4) * (x(i,2) - HALF)**m  &
                       - cb(5) * (x(i,3) + HALF)**m  &
                       + cb(6) * (x(i,3) - HALF)**m  )
  end do

end subroutine BoundaryPerturbation

!===============================================================================

end module ISP_Flow_Problem__No_Flow
