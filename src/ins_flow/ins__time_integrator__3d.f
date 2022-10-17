!> summary:  Base type of one-step time integrators for incompressible flows
!> author:   Joerg Stiller
!> date:     2022/09/09
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module INS__Time_Integrator__3D
  use Kind_Parameters
  use Constants
  use XMPI
  use SEM__Boundary_Variable__3D
  use INS__Problem__3D
  use INS__Operator__3D
  implicit none
  private

  public :: INS_TimeIntegrator_3D
  public :: INS_TimeIntegratorOptions_3D

  !-----------------------------------------------------------------------------
  !> Abstract type of a one-step time integrator for incompressible flow

  type, abstract :: INS_TimeIntegrator_3D

    class(INS_Problem_3D),  pointer :: problem => null() !< flow problem
    class(INS_Operator_3D), pointer :: ins_op  => null() !< INS operators

    character(len=80) :: name = ''  !< time-integrator name

  contains
    procedure, non_overridable    :: Init_INS_TimeIntegrator_3D
    procedure                     :: ProjectionStep
    procedure(TimeStep), deferred :: TimeStep
  end type INS_TimeIntegrator_3D

  !=============================================================================
  ! external module procedures

  interface

    !---------------------------------------------------------------------------
    !> extrapolation-projection-diffusion step for incompressible flow

    module subroutine ProjectionStep( this, tau, v_0, F_c, F_d, Q, bv_u, u &
                                    , i_max_p, i_max_v, r_red, r_max       )

      class(INS_TimeIntegrator_3D),   intent(in)    :: this
      real(RNP),                      intent(in)    :: tau
      real(RNP), contiguous,          intent(in)    :: v_0(:,:,:,:,:)
      real(RNP), contiguous,          intent(in)    :: F_c(:,:,:,:,:)
      real(RNP), contiguous,          intent(in)    :: F_d(:,:,:,:,:)
      real(RNP), contiguous,          intent(in)    :: Q(:,:,:,:,:)
      class(SEM_BoundaryVariable_3D), intent(inout) :: bv_u(:)
      real(RNP), contiguous,          intent(out)   :: u(:,:,:,:,:)
      integer,                        intent(in)    :: i_max_p
      integer,                        intent(in)    :: i_max_v
      real(RNP),            optional, intent(in)    :: r_red
      real(RNP),            optional, intent(in)    :: r_max

    end subroutine ProjectionStep

  end interface

  !=============================================================================
  ! deferred module procedures

  abstract interface

    !---------------------------------------------------------------------------
    !> Execution of a single time step
    !>
    !> Activate the `standby` option to allow for reusing workspace and data.

    subroutine TimeStep(this, t, dt, u, standby)
      import
      class(INS_TimeIntegrator_3D), intent(inout) :: this
      real(RNP), intent(inout) :: t  !< time t₀ → t
      real(RNP), intent(in)    :: dt !< step size ∆t = t-t₀
      real(RNP), contiguous, intent(inout) :: u(:,:,:,:,:) !< u(x,t₀) → u(x,t)
      logical, optional, intent(in) :: standby !< reuse workspace T/F [F]
    end subroutine TimeStep

  end interface

  !-----------------------------------------------------------------------------
  !> Base type for providing time integrator options

  type INS_TimeIntegratorOptions_3D
  contains
    procedure :: Bcast => Bcast_INS_TimeIntegratorOptions_3D
  end type INS_TimeIntegratorOptions_3D

contains

  !=============================================================================
  ! INS_TimeIntegrator_3D: type-bound procedures

  !-----------------------------------------------------------------------------
  !> Initialization of TimeIntegrator object

  subroutine Init_INS_TimeIntegrator_3D(this, problem, ins_op, opt)
    class(INS_TimeIntegrator_3D),                  intent(inout) :: this
    class(INS_Problem_3D),                 target, intent(in)    :: problem
    class(INS_Operator_3D),                target, intent(in)    :: ins_op
    class(INS_TimeIntegratorOptions_3D), optional, intent(in)    :: opt

    this % problem => problem
    this % ins_op  => ins_op

    ! options
    if (present(opt)) then
      return ! nothing, so far
    end if

  end subroutine Init_INS_TimeIntegrator_3D

  !=============================================================================
  ! TimeIntegratorOptions: type-bound procedures

  !-----------------------------------------------------------------------------
  !> MPI broadcasting of time-integrator options

  subroutine Bcast_INS_TimeIntegratorOptions_3D(this, root, comm)
    class(INS_TimeIntegratorOptions_3D), intent(inout) :: this
    integer,        intent(in) :: root !< rank of broadcast root
    type(MPI_Comm), intent(in) :: comm !< MPI communicator

    type(MPI_Request) :: request(0)
    type(MPI_Status)  :: stat(size(request))
    integer :: n

    n = 1
!   call XMPI_Ibcast( this % <opt1>, root, comm, request(n) );  n = n + 1
!   call XMPI_Ibcast( this % <opt2>, root, comm, request(n) )
    if (size(request) == 0) return

    call MPI_Waitall(n, request, stat)

  end subroutine Bcast_INS_TimeIntegratorOptions_3D

  !=============================================================================

end module INS__Time_Integrator__3D
