!> summary:  3D Transition of disturbed Taylor-Green vortex
!> author:   Joerg Stiller, Montadhar Guesmi
!> date:     2020/02/18
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module ISP_Flow_Problem__Transition_TG

  use Kind_Parameters,   only: RNP
  use Constants,         only: ZERO, ONE, HALF, PI
  use Execution_Control
  use Array_Assignments
  use XMPI

  use ISP_Flow_Problem

  implicit none
  private

  public :: FlowProblem_Transition_TG

  !-----------------------------------------------------------------------------
  !> Type for defining and handling the problem

  type, extends(FlowProblem) :: FlowProblem_Transition_TG
  contains

    procedure :: SetProblem
    procedure :: GetInitialValues
    procedure :: GetBoundaryValues
    procedure :: GetBoundaryTimeDerivative
    procedure :: GetExternalSources

  end type FlowProblem_Transition_TG

contains

!===============================================================================
! type bound procedures

!-------------------------------------------------------------------------------
!> Initialization

subroutine SetProblem(problem, file, comm)
  class(FlowProblem_Transition_TG), intent(inout) :: problem
  character(len=*), optional, intent(in) :: file  !< input file
  type(MPI_Comm),   optional, intent(in) :: comm  !< MPI communicator

  ! local variables ............................................................

  real(RNP) :: nu     = 0.01 ! kinematic viscosity, ν = 1/Re
  real(RNP) :: nu_svv = 0.0  ! spectral diffusivity amplitude
  character, allocatable :: bc(:,:)

  namelist /parameters/ nu, nu_svv

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
  allocate(bc(6, problem % nc), source = 'P')

  if (rank == 0 .and. exists) then
    read(prm, nml=parameters)
    close(prm)
  end if

  if (present(comm)) then
    call XMPI_Bcast(nu,     0, comm)
    call XMPI_Bcast(nu_svv, 0, comm)
  end if

  problem % stokes         = .false.   ! Stokes flow
  problem % exact_solution = .false.   ! exact solution is not provided

  allocate(problem % nu_ref(     problem%nc ), source = ZERO)
  allocate(problem % nu_svv_ref( problem%nc ), source = ZERO)
  problem % nu_ref(1:3)     = nu
  problem % nu_svv_ref(1:3) = nu_svv
  problem % x0 = -PI
  problem % x1 =  PI

  call move_alloc(bc, problem % bc)

end subroutine SetProblem

!-------------------------------------------------------------------------------
!> Provides the initial values u(x,0) for mesh points x

subroutine GetInitialValues(problem, x, u)
  class(FlowProblem_Transition_TG), intent(in) :: problem
  real(RNP), intent(in)  :: x(:,:,:,:,:) !< mesh points
  real(RNP), intent(out) :: u(:,:,:,:,:) !< flow variables

  integer :: m, n

  n = size(x(:,:,:,:,1))

  call GetVelocity(n, x, u(:,:,:,:,1:3))

  ! remaining variables get zero
  do m = 4, size(u,5)
    call SetArray(u(:,:,:,:,m), ZERO)
  end do

end subroutine GetInitialValues

!-------------------------------------------------------------------------------
!> Provides the values `ub = u(xb,t)` for points `xb` on boundary `b`

subroutine GetBoundaryValues(problem, b, xb, t, ub)
  class(FlowProblem_Transition_TG), intent(in) :: problem
  integer,   intent(in)  :: b             !< boundary ID
  real(RNP), intent(in)  :: xb(:,:,:,:)   !< mesh boundary points
  real(RNP), intent(in)  :: t             !< time
  real(RNP), intent(out) :: ub(:,:,:,:)   !< flow variables

  call Warning( 'GetBoundaryValues'               &
              , 'No boundary values available'    &
              , 'ISP_Flow_Problem__Transition_TG' )

  call SetArray(ub, ZERO, multi=.true.)

  ! silence the compiler ;)
  if (b < 0 .or. size(xb) < 0 .or. t < 0) return

end subroutine GetBoundaryValues

!-------------------------------------------------------------------------------
!> Provides the values `dt_ub = ∂u/∂t(xb,t)` for points `xb` on boundary `b`

subroutine GetBoundaryTimeDerivative(problem, b, xb, t, dt_ub)
  class(FlowProblem_Transition_TG), intent(in)  :: problem
  integer,   intent(in)  :: b              !< boundary ID
  real(RNP), intent(in)  :: xb(:,:,:,:)    !< boundary points
  real(RNP), intent(in)  :: t              !< time
  real(RNP), intent(out) :: dt_ub(:,:,:,:) !< ∂u/∂t

  call Warning( 'GetBoundaryTimeDerivative'       &
              , 'No boundary values available'    &
              , 'ISP_Flow_Problem__Transition_TG' )

  call SetArray(dt_ub, ZERO, multi=.true.)

  ! silence the compiler ;)
  if (b < 0 .or. size(xb) < 0 .or. t < 0) return

end subroutine GetBoundaryTimeDerivative

!-------------------------------------------------------------------------------
!> Provides the external sources for all variables at points x and time t

subroutine GetExternalSources(problem, x, t, f)
  class(FlowProblem_Transition_TG), intent(in) :: problem
  real(RNP), intent(in)  :: x(:,:,:,:,:) !< mesh points
  real(RNP), intent(in)  :: t            !< time
  real(RNP), intent(out) :: f(:,:,:,:,:) !< external sources

  call SetArray(f, ZERO, multi=.true.)

end subroutine GetExternalSources

!===============================================================================
! problem-specific procedures

!-------------------------------------------------------------------------------
!> Velocity

subroutine GetVelocity(n, x, v)
  integer,   intent(in)  :: n      !< number of points
  real(RNP), intent(in)  :: x(n,3) !< mesh points
  real(RNP), intent(out) :: v(n,3) !< velocity at mesh points

  real(RNP) :: x1, x2, x3 ! coordinates xi
  integer   :: i

  do i = 1, n
    x1 = x(i,1)
    x2 = x(i,2)
    x3 = x(i,3)
    v(i,1) =  cos(x3) * sin(x1) * cos(x2)
    v(i,2) = -cos(x3) * sin(x2) * cos(x1)
    v(i,3) =  ZERO
  end do

end subroutine GetVelocity

!===============================================================================

end module ISP_Flow_Problem__Transition_TG
