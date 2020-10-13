!> summary:  Base type of one-step time integrators for incompressible flows
!> author:   Joerg Stiller
!> date:     2020/03/05
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module CART__ISP_Flow__Time_Integrator
  use, intrinsic :: ISO_Fortran_Env, only: OUTPUT_UNIT
  use Kind_Parameters, only: RNP
  use Constants,       only: HALF
  use XMPI
  use ISP_Flow_Problem
  use CART__ISP_Flow__Operators
  implicit none
  private

  public :: TimeIntegrator
  public :: TimeIntegratorOptions

  !-----------------------------------------------------------------------------
  !> Abstract type of a one-step time integrator for incompressible flow
  !>
  !> Predefined options
  !>
  !>   * `splitting` -- defines the splitting scheme used in the stages:     \n
  !>        1: standard velocity correction with χ = -1                      \n
  !>        2: rotational velocity correction with χ = -2                    \n
  !>           and F_d3/2 removed after extrapolation                        \n
  !>        3: "native" velocity correction with χ = 0                       \n
  !>        4: velocity correction with χ = -1                               \n
  !>           and F_d3 removed after extrapolation                          \n
  !>        5: velocity correction with χ = -2                               \n
  !>           and F_d3 removed after extrapolation
  !>
  !>   * `project`                                                           \n
  !>        0: no additional projection step                                 \n
  !>        1: additional projection at the end of the time step             \n
  !>        2: additional projection at the end of each stage (RK only)
  !>
  !>   * `pressure`
  !>        0: return pressure as computed                                   \n
  !>        1: recompute pressure at the end of the time step

  type, abstract :: TimeIntegrator

    class(FlowProblem),   pointer :: problem => null() !< flow problem
    class(FlowOperators), pointer :: flow_op => null() !< flow operators

    character(len=80) :: name = ''  !< time-integrator name

    ! options
    integer :: splitting  = 2 !< selector for stage splitting scheme
    integer :: project    = 0 !< switch for additional projections
    integer :: pressure   = 0 !< switch for pressure recomputation

    ! derived parameters
    real(RNP) :: chi
    real(RNP) :: cd3

  contains
    procedure, non_overridable :: Init_TimeIntegrator
    procedure, non_overridable :: Show_TimeIntegrator
    procedure :: Show => Show_TimeIntegrator
    procedure(TimeStep), deferred :: TimeStep
  end type TimeIntegrator

  abstract interface

    !---------------------------------------------------------------------------
    !> Execution of a single time step
    !>
    !> Activate the `standby` option to allow for reusing workspace and data

    subroutine TimeStep(this, t, dt, u, standby)
      import
      class(TimeIntegrator), intent(inout) :: this
      real(RNP),         intent(inout) :: t            !< time t₀ → t
      real(RNP),         intent(in)    :: dt           !< step size ∆t = t-t₀
      real(RNP),         intent(inout) :: u(:,:,:,:,:) !< u(x,t₀) → u(x,t)
      logical, optional, intent(in)    :: standby      !< reuse workspace
    end subroutine TimeStep

  end interface

  !-----------------------------------------------------------------------------
  !> Base type for providing time integrator options

  type TimeIntegratorOptions
    integer :: splitting  = 2 !< selector for stage splitting scheme
    integer :: project    = 0 !< switch for additional projections
    integer :: pressure   = 0 !< switch for pressure recomputation
  contains
    procedure :: Bcast => Bcast_TimeIntegratorOptions
  end type TimeIntegratorOptions

contains

  !=============================================================================
  ! TimeIntegrator: type-bound procedures

  !-----------------------------------------------------------------------------
  !> Initialization of TimeIntegrator object

  subroutine Init_TimeIntegrator(this, problem, flow_op, opt)
    class(TimeIntegrator),        intent(inout) :: this
    class(FlowProblem),   target, intent(in)    :: problem    !< flow problem
    class(FlowOperators), target, intent(in)    :: flow_op    !< flow operators
    class(TimeIntegratorOptions), optional, intent(in) :: opt !< options

    this % problem => problem
    this % flow_op => flow_op

    ! options
    if (present(opt)) then
      this % splitting  = opt % splitting
      this % project    = opt % project
      this % pressure   = opt % pressure
    end if

    ! derived parameters
    select case(this % splitting)
    case(1)
      this % chi = -1
      this % cd3 =  0
    case(2)
      this % chi = -2
      this % cd3 =  HALF
    case(3)
      this % chi =  0
      this % cd3 =  0
    case(4)
      this % chi = -1
      this % cd3 =  1
    case(5)
      this % chi = -2
      this % cd3 =  1
    end select

  end subroutine Init_TimeIntegrator

  !-----------------------------------------------------------------------------
  !> Output of TimeIntegrator settings

  subroutine Show_TimeIntegrator(this, unit)
    class(TimeIntegrator), intent(in) :: this
    integer,     optional, intent(in) :: unit  !< output unit

    integer :: io

    if (present(unit)) then
      io = unit
    else
      io = OUTPUT_UNIT
    end if

    write(io,'(/,A)')       'TimeIntegrator settings'
    write(io,'(A,/)')       repeat('=',80)
    write(io,'(2X,A,T15,A)')  'name:'     , trim(this % name)
    write(io,'(2X,A,T15,I0)') 'splitting:', this % splitting
    write(io,'(2X,A,T15,I0)') 'project:'  , this % project
    write(io,'(2X,A,T15,I0)') 'pressure:' , this % pressure
    write(io,*)

    ! append further settings in corresponding routines of derived types

  end subroutine Show_TimeIntegrator

  !=============================================================================
  ! TimeIntegratorOptions: type-bound procedures

  !-----------------------------------------------------------------------------
  !> MPI broadcasting of time-integrator options

  subroutine Bcast_TimeIntegratorOptions(this, root, comm)
    class(TimeIntegratorOptions), intent(inout) :: this
    integer,        intent(in) :: root !< rank of broadcast root
    type(MPI_Comm), intent(in) :: comm !< MPI communicator

    type(MPI_Request) :: request(3)
    type(MPI_Status)  :: stat(size(request))
    integer :: n

    n = 1
    call XMPI_Ibcast( this % splitting, root, comm, request(n) );  n = n + 1
    call XMPI_Ibcast( this % project  , root, comm, request(n) );  n = n + 1
    call XMPI_Ibcast( this % pressure , root, comm, request(n) )

    call MPI_Waitall(n, request, stat)

  end subroutine Bcast_TimeIntegratorOptions

  !=============================================================================

end module CART__ISP_Flow__Time_Integrator
