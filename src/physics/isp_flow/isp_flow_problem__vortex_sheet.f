!> summary:  Roll-up of a vortex sheet from Brown & Minion, JCP 122, 1995
!> author:   Joerg Stiller
!> date:     2018/05/20
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module ISP_Flow_Problem__Vortex_Sheet

  use Kind_Parameters,   only: RNP
  use Constants,         only: ZERO, HALF, PI
  use Execution_Control
  use Array_Assignments
  use XMPI

  use ISP_Flow_Problem

  implicit none
  private

  public :: FlowProblem_VortexSheet

  !-----------------------------------------------------------------------------
  !> Type for defining and handling the Hesthaven-Warburton vortex problem

  type, extends(FlowProblem) :: FlowProblem_VortexSheet

    real(RNP) :: kappa = 30    !< thickness parameter
    real(RNP) :: delta = 0.05  !< perturbation amplitude

  contains

    procedure :: SetProblem
    procedure :: GetInitialValues
    procedure :: GetBoundaryValues
    procedure :: GetBoundaryTimeDerivative
    procedure :: GetExternalSources

  end type FlowProblem_VortexSheet

contains

!===============================================================================
! type bound procedures

!-------------------------------------------------------------------------------
!> Initialization

subroutine SetProblem(problem, file, comm)
  class(FlowProblem_VortexSheet), intent(inout) :: problem
  character(len=*), optional, intent(in) :: file    !< input file
  type(MPI_Comm),   optional, intent(in) :: comm    !< MPI communicator

  ! local variables ............................................................

  real(RNP) :: kappa = 30    ! thickness parameter
  real(RNP) :: delta = 0.05  ! perturbation amplitude
  real(RNP) :: Re    = 1e5   ! Reynolds number

  namelist /parameters/ kappa, delta, Re

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

  if (rank == 0 .and. exists) then
    read(prm, nml=parameters)
    close(prm)
  end if

  if (present(comm)) then
    call XMPI_Bcast(kappa , 0, comm)
    call XMPI_Bcast(delta , 0, comm)
    call XMPI_Bcast(Re    , 0, comm)
  end if

  allocate(problem % nu_ref( problem%nc ), source = ZERO)

  problem % kappa = kappa
  problem % delta = delta
  problem % nu_ref(1:3) = 1 / Re
  problem % x0 = 0
  problem % x1 = 1

  allocate(problem % bc(6, problem % nc), source = 'P')

end subroutine SetProblem

!-------------------------------------------------------------------------------
!> Provides the initial values u(x,0) for mesh points x

subroutine GetInitialValues(problem, x, u)
  class(FlowProblem_VortexSheet), intent(in) :: problem
  real(RNP), intent(in)  :: x(:,:,:,:,:) !< mesh points
  real(RNP), intent(out) :: u(:,:,:,:,:) !< flow variables

  integer :: m, n

  n = size(x(:,:,:,:,1))

  call GetInitialVelocity(n, x, problem%kappa, problem%delta, u(:,:,:,:,1:3))

  ! pressure and remaining variables get zero
  do m = 4, size(u,5)
    call SetArray(u(:,:,:,:,m), ZERO)
  end do

end subroutine GetInitialValues

!-------------------------------------------------------------------------------
!> Boundary values (not used)

subroutine GetBoundaryValues(problem, b, xb, t, ub)
  class(FlowProblem_VortexSheet), intent(in) :: problem
  integer,   intent(in)  :: b             !< boundary ID
  real(RNP), intent(in)  :: xb(:,:,:,:)   !< mesh boundary points
  real(RNP), intent(in)  :: t             !< time
  real(RNP), intent(out) :: ub(:,:,:,:)   !< flow variables

  call SetArray(ub, ZERO, multi=.true.)

end subroutine GetBoundaryValues

!-------------------------------------------------------------------------------
!> Boundary time derivative (not used)

subroutine GetBoundaryTimeDerivative(problem, b, xb, t, dt_ub)
  class(FlowProblem_VortexSheet), intent(in) :: problem
  integer,   intent(in)  :: b              !< boundary ID
  real(RNP), intent(in)  :: xb(:,:,:,:)    !< boundary points
  real(RNP), intent(in)  :: t              !< time
  real(RNP), intent(out) :: dt_ub(:,:,:,:) !< ∂u/∂t

  call SetArray(dt_ub, ZERO, multi=.true.)

end subroutine GetBoundaryTimeDerivative

!-------------------------------------------------------------------------------
!> Provides the external sources for all variables at points x and time t

subroutine GetExternalSources(problem, x, t, f)
  class(FlowProblem_VortexSheet), intent(in) :: problem
  real(RNP), intent(in)  :: x(:,:,:,:,:) !< mesh points
  real(RNP), intent(in)  :: t            !< time
  real(RNP), intent(out) :: f(:,:,:,:,:) !< external sources

  call SetArray(f, ZERO, multi=.true.)

end subroutine GetExternalSources

!===============================================================================
! problem-specific procedures

!-------------------------------------------------------------------------------
!> Velocity

subroutine GetInitialVelocity(n, x, kappa, delta, v)
  integer,   intent(in)  :: n       !< number of points
  real(RNP), intent(in)  :: x(n,3)  !< mesh points
  real(RNP), intent(in)  :: kappa   !< thickness parameter
  real(RNP), intent(in)  :: delta   !< perturbation amplitude
  real(RNP), intent(out) :: v(n,3)  !< velocity at mesh points

  integer :: i

  do i = 1, n
    if (x(i,2) < HALF) then
      v(i,1) = tanh( kappa * (x(i,2) - 0.25) )
    else
      v(i,1) = tanh( kappa * (0.75 - x(i,2)) )
    end if
    v(i,2) =  delta * sin(2 * PI * x(i,1))
    v(i,3) =  ZERO
  end do

end subroutine GetInitialVelocity

!===============================================================================

end module ISP_Flow_Problem__Vortex_Sheet
