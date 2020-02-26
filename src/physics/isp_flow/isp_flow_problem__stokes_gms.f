!> summary:  Stokes test case of Guermond, Minev & Shen (2006)
!> author:   Joerg Stiller
!> date:     2018/03/22
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Stokes test case of Guermond, Minev & Shen (2006)
!>
!> Test case from
!>   J. Guermond, P. Minev & J. Shen,
!>   An overview of projection methods for incompressible flows,
!>   Computer Methods in Applied Mechanics and Engineering,
!>   195(44-47):6011–6045, 2006
!> on page 6018, equation (30).
!>
!===============================================================================

module ISP_Flow_Problem__Stokes_GMS

  use Kind_Parameters,   only: RNP
  use Constants,         only: ZERO, PI
  use Execution_Control
  use Array_Assignments
  use XMPI

  use ISP_Flow_Problem

  implicit none
  private

  public :: FlowProblem_Stokes_GMS

  !-----------------------------------------------------------------------------
  !> Type for defining and handling the Stokes problem

  type, extends(FlowProblem) :: FlowProblem_Stokes_GMS

  contains

    procedure :: SetProblem
    procedure :: GetInitialValues
    procedure :: GetBoundaryValues
    procedure :: GetBoundaryTimeDerivative
    procedure :: GetExternalSources
    procedure :: GetExactSolution
    procedure :: GetExactTimeDerivative

  end type FlowProblem_Stokes_GMS

contains

!===============================================================================
! type bound procedures

!-------------------------------------------------------------------------------
!> Initialization

subroutine SetProblem(problem, file, comm)
  class(FlowProblem_Stokes_GMS), intent(inout) :: problem
  character(len=*), optional, intent(in) :: file    !< input file
  type(MPI_Comm),   optional, intent(in) :: comm    !< MPI communicator

  ! local variables ............................................................

  logical   :: stokes = .true.
  real(RNP) :: nu     =  1  ! kinematic viscosity
  real(RNP) :: x0(3)  = -1  ! bounding box: corner nearest to -infinity
  real(RNP) :: x1(3)  =  1  ! bounding box: corner nearest to +infinity
  character, allocatable :: bc(:,:)

  namelist /parameters/ x0, x1, bc

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
    call XMPI_Bcast(x0, 0, comm)
    call XMPI_Bcast(x1, 0, comm)
    call XMPI_Bcast(bc, 0, comm)
  end if

  problem % stokes = stokes
  problem % exact_solution = .true.

  problem % v_ref = PI

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
  class(FlowProblem_Stokes_GMS), intent(in) :: problem
  real(RNP), intent(in)  :: x(:,:,:,:,:) !< mesh points
  real(RNP), intent(out) :: u(:,:,:,:,:) !< flow variables

  integer :: m, n

  n = size(x(:,:,:,:,1))

  call GetVelocity(n, x, ZERO, u(:,:,:,:,1:3))
  call GetPressure(n, x, ZERO, u(:,:,:,:,4)  )

  ! remaining variables get zero
  do m = 5, size(u,5)
    call SetArray(u(:,:,:,:,m), ZERO)
  end do

end subroutine GetInitialValues

!-------------------------------------------------------------------------------
!> Provides the values `ub = u(xb,t)` for points `xb` on boundary `b`

subroutine GetBoundaryValues(problem, b, xb, t, ub)
  class(FlowProblem_Stokes_GMS), intent(in) :: problem
  integer,   intent(in)  :: b             !< boundary ID
  real(RNP), intent(in)  :: xb(:,:,:,:)   !< mesh boundary points
  real(RNP), intent(in)  :: t             !< time
  real(RNP), intent(out) :: ub(:,:,:,:)   !< flow variables

  integer :: m, n

  n = size(xb(:,:,:,1))

  call GetVelocity(n, xb, t, ub(:,:,:,1:3))
  call GetPressure(n, xb, t, ub(:,:,:,4)  )

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
  class(FlowProblem_Stokes_GMS), intent(in) :: problem
  integer,   intent(in)  :: b              !< boundary ID
  real(RNP), intent(in)  :: xb(:,:,:,:)    !< boundary points
  real(RNP), intent(in)  :: t              !< time
  real(RNP), intent(out) :: dt_ub(:,:,:,:) !< ∂u/∂t

  integer :: m, n

  n = size(xb(:,:,:,1))

  call GetVelocityTimeDerivative(n, xb, t, dt_ub(:,:,:,1:3))
  call GetPressureTimeDerivative(n, xb, t, dt_ub(:,:,:,4)  )

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
  class(FlowProblem_Stokes_GMS), intent(in) :: problem
  real(RNP), intent(in)  :: x(:,:,:,:,:) !< mesh points
  real(RNP), intent(in)  :: t            !< time
  real(RNP), intent(out) :: f(:,:,:,:,:) !< external sources

  integer :: m, n

  n = size(x(:,:,:,:,1))

  if (problem % stokes) then
    call GetStokesSource(size(x(:,:,:,:,1)), x, t, f(:,:,:,:,1:3))
    do m = 4, size(f,5)
      call SetArray(f(:,:,:,:,m), ZERO)
    end do
  else
    ! should not happen !
  end if

end subroutine GetExternalSources

!---------------------------------------------------------------------------
!> Provides the exact solution u(x,t)

subroutine GetExactSolution(problem, x, t, u)
  class(FlowProblem_Stokes_GMS), intent(in) :: problem
  real(RNP), intent(in)  :: x(:,:,:,:,:) !< mesh points
  real(RNP), intent(in)  :: t            !< time
  real(RNP), intent(out) :: u(:,:,:,:,:) !< solution

  integer :: m, n

  n = size(x(:,:,:,:,1))

  call GetVelocity(n, x, t, u(:,:,:,:,1:3))
  call GetPressure(n, x, t, u(:,:,:,:,4)  )

  ! remaining variables get zero
  do m = 5, size(u,5)
    call SetArray(u(:,:,:,:,m), ZERO)
  end do

end subroutine GetExactSolution

!---------------------------------------------------------------------------
!> Provides the time derivative of the exact solution, ∂u/∂t(x,t)

subroutine GetExactTimeDerivative(problem, x, t, dt_u)
  class(FlowProblem_Stokes_GMS), intent(in) :: problem
  real(RNP), intent(in)  :: x(:,:,:,:,:)    !< mesh points
  real(RNP), intent(in)  :: t               !< time
  real(RNP), intent(out) :: dt_u(:,:,:,:,:) !< time derivative

  integer :: m, n

  n = size(x(:,:,:,:,1))

  call GetVelocityTimeDerivative(n, x, t, dt_u(:,:,:,:,1:3))
  call GetPressureTimeDerivative(n, x, t, dt_u(:,:,:,:,4)  )

  ! remaining variables get zero
  do m = 5, size(dt_u,5)
    call SetArray(dt_u(:,:,:,:,m), ZERO)
  end do

end subroutine GetExactTimeDerivative

!===============================================================================
! problem-specific procedures

!-------------------------------------------------------------------------------
!> Velocity

subroutine GetVelocity(n, x, t, v)
  integer,   intent(in)  :: n       !< number of points
  real(RNP), intent(in)  :: x(n,3)  !< mesh points
  real(RNP), intent(in)  :: t       !< time
  real(RNP), intent(out) :: v(n,3)  !< velocity at mesh points

  real(RNP) :: c, x1, x2
  integer   :: i

  c = PI * sin(t)

  do i = 1, n
    x1 = x(i,1)
    x2 = x(i,2)
    v(i,1) =  c * sin(2*PI*x2) * sin(PI*x1)**2
    v(i,2) = -c * sin(2*PI*x1) * sin(PI*x2)**2
    v(i,3) =  ZERO
  end do

end subroutine GetVelocity

!-------------------------------------------------------------------------------
!> Velocity time derivative

subroutine GetVelocityTimeDerivative(n, x, t, dt_v)
  integer,   intent(in)  :: n          !< number of points
  real(RNP), intent(in)  :: x(n,3)     !< mesh points
  real(RNP), intent(in)  :: t          !< time
  real(RNP), intent(out) :: dt_v(n,3)  !< velocity at mesh points

  real(RNP) :: c, x1, x2
  integer   :: i

  c = PI * cos(t)

  do i = 1, n
    x1 = x(i,1)
    x2 = x(i,2)
    dt_v(i,1) =  c * sin(2*PI*x2) * sin(PI*x1)**2
    dt_v(i,2) = -c * sin(2*PI*x1) * sin(PI*x2)**2
    dt_v(i,3) =  ZERO
  end do

end subroutine GetVelocityTimeDerivative

!-------------------------------------------------------------------------------
!> Pressure

subroutine GetPressure(n, x, t, p)
  integer,   intent(in)  :: n       !< number of points
  real(RNP), intent(in)  :: x(n,3)  !< mesh points
  real(RNP), intent(in)  :: t       !< time
  real(RNP), intent(out) :: p(n)    !< pressure at mesh points

  real(RNP) :: c
  integer   :: i

  c = sin(t)

  do i = 1, n
    p(i) = c * cos(PI * x(i,1)) * sin(PI * x(i,2))
  end do

end subroutine GetPressure

!-------------------------------------------------------------------------------
!> Pressure time derivative

subroutine GetPressureTimeDerivative(n, x, t, dt_p)
  integer,   intent(in)  :: n       !< number of points
  real(RNP), intent(in)  :: x(n,3)  !< mesh points
  real(RNP), intent(in)  :: t       !< time
  real(RNP), intent(out) :: dt_p(n) !< pressure time derivative

  real(RNP) :: c
  integer   :: i

  c = cos(t)

  do i = 1, n
    dt_p(i) = c * cos(PI * x(i,1)) * sin(PI * x(i,2))
  end do

end subroutine GetPressureTimeDerivative

!-------------------------------------------------------------------------------
!> Stokes source

subroutine GetStokesSource(n, x, t, f)
  integer,   intent(in)  :: n       !< number of mesh points
  real(RNP), intent(in)  :: x(n,3)  !< mesh points
  real(RNP), intent(in)  :: t       !< time
  real(RNP), intent(out) :: f(n,3)  !< source

  real(RNP) :: x1, x2
  integer   :: i

  do i = 1, n
    x1 = x(i,1)
    x2 = x(i,2)

    f(i,1) =  PI * cos(t) * ( sin(2*PI*x2) * sin(PI*x1)**2 )   &!   ∂u/∂t
           -  2 * PI**3 * sin(t)                               &! - Lapace(u)
                        * ( cos(2*PI*x1) * sin(2*PI*x2)        &!   ..
                          - 2 * sin(2*PI*x2) * sin(PI*x1)**2 ) &!   ..
           -  PI * sin(t) * sin(PI*x1) * sin(PI*x2)             ! + ∂p/∂x

    f(i,2) = -PI * cos(t) * (sin(2*PI*x1) * sin(PI*x2)**2 )    &!   ∂v/∂t
           -  2 * PI**3 * sin(t)                               &! - Lapace(v)
                        * (-sin(2*PI*x1) * cos(2*PI*x2)        &!   ..
                          + 2 * sin(2*PI*x1) * sin(PI*x2)**2 ) &!   ..
           +  PI * sin(t) * cos(PI*x1) * cos(PI*x2)             ! + ∂p/∂y

    f(i,3) = ZERO

  end do

end subroutine GetStokesSource

!===============================================================================

end module ISP_Flow_Problem__Stokes_GMS
