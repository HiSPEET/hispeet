!> summary:  3D Transition of disturbed Taylor-Green vortex
!> author:   
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
  !> Type for defining and handling the Minion-Saye vortex problem

  type, extends(FlowProblem) :: FlowProblem_Transition_TG

    !real(RNP) :: vt(3) = 0  !<  Transitional Velocity 
    !real(RNP) :: xt(3) = 0  !<  Initial Displacement
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

  real(RNP) :: nu     =  0.01         ! kinematic viscosity
  !real(RNP) :: vt(3)  =  0            ! velocity field 
  !real(RNP) :: xt(3)  =  0            ! initial displacement
  real(RNP) :: x0(3)  = -HALF         ! bounding box: corner nearest to -∞
  real(RNP) :: x1(3)  =  HALF         ! bounding box: corner nearest to +∞
  character, allocatable :: bc(:,:)

  namelist /parameters/  nu, x0, x1, bc ! < have to be controlled !!!

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
    call XMPI_Bcast(nu    , 0, comm)
    call XMPI_Bcast(x0    , 0, comm)   
    call XMPI_Bcast(x1    , 0, comm)   
    call XMPI_Bcast(bc    , 0, comm)
  end if

  problem % stokes         = .false.   ! < Stokes Flow ??
  problem % exact_solution = .false.   ! < Exact Solution is not provided

  allocate(problem % nu_ref( problem%nc ), source = ZERO)
  problem % nu_ref(1:3) = nu
  problem % x0  = x0
  problem % x1  = x1

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

  associate(nu => problem%nu_ref(1))
    call GetVelocity(n, x, u(:,:,:,:,1:3))
  end associate

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

  call SetArray(ub, ZERO, multi=.true.)
   
  ! silence the compiler ;)
  if (b < 0) return

end subroutine GetBoundaryValues

!-------------------------------------------------------------------------------
!> Provides the values `dt_ub = ∂u/∂t(xb,t)` for points `xb` on boundary `b`

subroutine GetBoundaryTimeDerivative(problem, b, xb, t, dt_ub)
  class(FlowProblem_Transition_TG), intent(in)  :: problem
  integer,   intent(in)  :: b              !< boundary ID
  real(RNP), intent(in)  :: xb(:,:,:,:)    !< boundary points
  real(RNP), intent(in)  :: t              !< time
  real(RNP), intent(out) :: dt_ub(:,:,:,:) !< ∂u/∂t

  call SetArray(dt_ub, ZERO, multi=.true.)
  
  ! silence the compiler ;)
  if (b < 0) return

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

!subroutine GetVelocity(n, x, v)
subroutine GetVelocity(n, x, v)
  !arguments....................................................................
  integer,   intent(in)  :: n          !< number of points
  real(RNP), intent(in)  :: x(n,3)     !< mesh points
  real(RNP), intent(out) :: v(n,3)     !< velocity at mesh points
  !local variables..............................................................
  real(RNP)              :: c, L, K    !< constans 
  real(RNP)              :: x1, x2, x3 !< Space components xi
  real(RNP)              :: phi1, phi2, phi3
  integer                :: i
  !constants...................................................................
  K =  PI/HALF !< 2π 
  c = 1.0_RNP  !< U0: Velocity scale set 1 hier  
  L = 1.0_RNP  !< L : Length scale set 1 also hier ==> in this case Re=1/nu.
  do i = 1, n
    x1   = x(i,1)
    x2   = x(i,2)
    x3   = x(i,3)
    phi1 = K * x1 / L 
    phi2 = K * x2 / L
    phi3 = K * x3 / L
    v(i,1) =  cos(phi3) * sin( phi1)      * cos(phi2)
    v(i,2) =  cos(phi3) * sin(-phi2)      * cos(phi1)
    v(i,3) =  ZERO
  end do

end subroutine GetVelocity
!in case it is needed to initialize the pressure for t=0 
! from the work of Felix S. Schranner,  Vladyslav Rozov and Nikolaus A. Adams
! p(x,y,z) = 100 + 1/16 *[(cos(2x)+cos(2y))*(2+ cos(2z))-2]
!subroutine GetPressure(n,x,v) 
subroutine GetPressure(n, x, p)
  !arguments....................................................................
  integer,   intent(in)  :: n         !< number of points 
  real(RNP), intent(in)  :: x(n,3)    !< mesh points  
  real(RNP), intent(out) :: p(n)      !< pressure at mesh points
  !local variables..............................................................
  integer    :: i
  real(RNP)  :: x1, x2, x3, K, phi1, phi2, phi3, L
  !constants....................................................................
  K = PI/HALF 
  L = 1.0_RNP !< lenght scale set to 1. 
  do i=1, n
    x1   = x(i,1)
    x2   = x(i,2)
    x3   = x(i,3)
    phi1 = K * x1 / L 
    phi2 = K * x2 / L
    phi3 = K * x3 / L
    p(i)= (1/16) * ( (cos(2*phi1) + cos(2*phi2)) * (2 + cos(2*phi3)) - 2 )
  end do 
end subroutine GetPressure 

!===============================================================================

end module ISP_Flow_Problem__Transition_TG
