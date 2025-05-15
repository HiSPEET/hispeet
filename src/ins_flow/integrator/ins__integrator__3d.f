!> summary:  Base type of one-step time integrators for incompressible flows
!> author:   Joerg Stiller
!> date:     2022/09/09
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module INS__Integrator__3D
  use Kind_Parameters
  use Constants
  use XMPI
  use Boundary_Variable__3D
  use INS__Problem__3D
  use INS__Operator__3D
  implicit none
  private

  public :: INS_Integrator_3D
  public :: INS_IntegratorOptions_3D

  !-----------------------------------------------------------------------------
  !> Abstract type of a one-step time integrator for incompressible flow

  type, abstract :: INS_Integrator_3D

    class(INS_Problem_3D),  pointer :: problem => null() !< flow problem
    class(INS_Operator_3D), pointer :: ins_op  => null() !< INS operators

    character(len=80) :: name = ''  !< time-integrator name

  contains
    procedure, non_overridable    :: Init_INS_Integrator_3D
    procedure(TimeStep), deferred :: TimeStep
  end type INS_Integrator_3D

  !=============================================================================
  ! deferred module procedures

  abstract interface

    !---------------------------------------------------------------------------
    !> Execution of a single time step
    !>
    !> Activate the `standby` option to allow for reusing workspace and data.

    subroutine TimeStep(this, t, dt, u, standby)
      import
      class(INS_Integrator_3D), intent(inout) :: this
      real(RNP), intent(inout) :: t  !< time t₀ → t
      real(RNP), intent(in)    :: dt !< step size ∆t = t-t₀
      real(RNP), contiguous, intent(inout) :: u(:,:,:,:,:) !< u(x,t₀) → u(x,t)
      logical, optional, intent(in) :: standby !< reuse workspace T/F [F]
    end subroutine TimeStep

  end interface

  !-----------------------------------------------------------------------------
  !> Base type for providing time integrator options

  type INS_IntegratorOptions_3D
  contains
    procedure :: Bcast => Bcast_INS_IntegratorOptions_3D
  end type INS_IntegratorOptions_3D

contains

  !=============================================================================
  ! INS_Integrator_3D: type-bound procedures

  !-----------------------------------------------------------------------------
  !> Initialization of Integrator object

  subroutine Init_INS_Integrator_3D(this, problem, ins_op, opt)
    class(INS_Integrator_3D),        intent(inout) :: this
    class(INS_Problem_3D),   target, intent(in)    :: problem
    class(INS_Operator_3D),  target, intent(in)    :: ins_op
    class(INS_IntegratorOptions_3D), intent(in)    :: opt

    this % problem => problem
    this % ins_op  => ins_op

  end subroutine Init_INS_Integrator_3D

  !=============================================================================
  ! IntegratorOptions: type-bound procedures

  !-----------------------------------------------------------------------------
  !> MPI broadcasting of time-integrator options

  subroutine Bcast_INS_IntegratorOptions_3D(this, root, comm)
    class(INS_IntegratorOptions_3D), intent(inout) :: this
    integer,        intent(in) :: root !< rank of broadcast root
    type(MPI_Comm), intent(in) :: comm !< MPI communicator
  end subroutine Bcast_INS_IntegratorOptions_3D

  !=============================================================================

end module INS__Integrator__3D
