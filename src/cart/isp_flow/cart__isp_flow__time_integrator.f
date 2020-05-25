!> summary:  Base type of one-step time integrators for incompressible flows
!> author:   Joerg Stiller
!> date:     2020/03/05
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module CART__ISP_Flow__Time_Integrator
  use Kind_Parameters, only: RNP
  use XMPI
  use ISP_Flow_Problem
  use CART__ISP_Flow__Operators
  implicit none
  private

  public :: TimeIntegrator
  public :: TimeIntegratorOptions

  !-----------------------------------------------------------------------------
  !> Abstract type of a one-step time integrator for incompressible flow

  type, abstract :: TimeIntegrator
    class(FlowProblem),   pointer :: problem => null() !< flow problem
    class(FlowOperators), pointer :: flow_op => null() !< flow operators
  contains
    procedure :: Init_TimeIntegrator
    procedure(TimeStep), deferred :: TimeStep
  end type TimeIntegrator

  abstract interface

    !---------------------------------------------------------------------------
    !> Execution of a single time step

    subroutine TimeStep(this, t, dt, u)
      import
      class(TimeIntegrator), intent(inout) :: this
      real(RNP), intent(inout) :: t             !< time t₀ → t
      real(RNP), intent(in)    :: dt            !< step size ∆t = t-t₀
      real(RNP), intent(inout) :: u (:,:,:,:,:) !< u(x,t₀) → u(x,t)
    end subroutine TimeStep

  end interface

  !-----------------------------------------------------------------------------
  !> Base type for providing time integrator options

  type TimeIntegratorOptions
  contains
    procedure :: Bcast => TimeIntegratorOptions_Bcast
  end type TimeIntegratorOptions

contains

  !=============================================================================
  ! TimeIntegrator: type-bound procedures

  !-----------------------------------------------------------------------------
  !> Initialization of TimeIntegrator object

  subroutine Init_TimeIntegrator(this, problem, flow_op)
    class(TimeIntegrator),        intent(inout) :: this
    class(FlowProblem),   target, intent(in)    :: problem    !< flow problem
    class(FlowOperators), target, intent(in)    :: flow_op    !< flow operators

    this % problem => problem
    this % flow_op => flow_op

  end subroutine Init_TimeIntegrator

  !=============================================================================
  ! TimeIntegratorOptions: type-bound procedures

  !-----------------------------------------------------------------------------
  !> MPI broadcasting of time-integrator options

  subroutine TimeIntegratorOptions_Bcast(this, root, comm)
    class(TimeIntegratorOptions), intent(inout) :: this
    integer,        intent(in) :: root !< rank of broadcast root
    type(MPI_Comm), intent(in) :: comm !< MPI communicator

    ! nothing to broadcast, so far
    if (root == 0 .or. comm % mpi_val == 0) return

  end subroutine TimeIntegratorOptions_Bcast

  !=============================================================================

end module CART__ISP_Flow__Time_Integrator
