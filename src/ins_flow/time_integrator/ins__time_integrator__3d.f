!> summary:  Base type of one-step time integrators for incompressible flows
!> author:   Joerg Stiller
!> date:     2022/09/09
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module INS__Time_Integrator__3D
  use Kind_Parameters
  use Constants
  use XMPI
  use Boundary_Variable__3D
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

    integer   :: i_krylov !< max num Krylov iterations
    integer   :: i_max_p  !< max num p-iterations in projection solver
    integer   :: i_max_v  !< max num v-iterations in projection solver
    integer   :: i_pre_p  !< max num p-iterations in projection preconditioner
    integer   :: i_pre_v  !< max num v-iterations in projection preconditioner
    real(RNP) :: r_red    !< min residual reduction, if > 0
    real(RNP) :: r_max    !< max residual to reach,  if > 0

    character(len=80) :: name = ''  !< time-integrator name

  contains
    procedure, non_overridable    :: Init_INS_TimeIntegrator_3D
    procedure                     :: FGMRES_Step
    procedure                     :: ProjectionStep
    procedure(TimeStep), deferred :: TimeStep
  end type INS_TimeIntegrator_3D

  !=============================================================================
  ! external module procedures

  interface

    !---------------------------------------------------------------------------
    !> extrapolation-projection-diffusion step for incompressible flow

    module subroutine ProjectionStep( this, tau, t, v_0, F_c, F_d, Q &
                                    , bv_u, mu, nu, u                &
                                    , i_max_p, i_max_v, r_red, r_max &
                                    , freeze                         )

      class(INS_TimeIntegrator_3D),    intent(in)    :: this
      real(RNP),                       intent(in)    :: tau
      real(RNP),                       intent(in)    :: t
      real(RNP), contiguous,           intent(in)    :: v_0(:,:,:,:,:)
      real(RNP), contiguous,           intent(in)    :: F_c(:,:,:,:,:)
      real(RNP), contiguous,           intent(in)    :: F_d(:,:,:,:,:)
      real(RNP), contiguous,           intent(in)    :: Q(:,:,:,:,:)
      class(BoundaryVariable_3D),      intent(in)    :: bv_u(:)
      real(RNP), contiguous, optional, intent(inout) :: mu(:,:,:,:)
      real(RNP), contiguous, optional, intent(inout) :: nu(:,:,:,:)
      real(RNP), contiguous,           intent(inout) :: u(:,:,:,:,:)
      integer,               optional, intent(in)    :: i_max_p
      integer,               optional, intent(in)    :: i_max_v
      real(RNP),             optional, intent(in)    :: r_red
      real(RNP),             optional, intent(in)    :: r_max
      logical,               optional, intent(in)    :: freeze

    end subroutine ProjectionStep

    !---------------------------------------------------------------------------
    !> FGMRES for Stokes part using projection step as a preconditioner

    module subroutine FGMRES_Step( this                                      &
                                 , tau, t, v_0, F_c, F_d, Q, bv_u, mu, nu, u )

      class(INS_TimeIntegrator_3D),    intent(in)    :: this
      real(RNP),                       intent(in)    :: tau
      real(RNP),                       intent(in)    :: t
      real(RNP), contiguous,           intent(in)    :: v_0(:,:,:,:,:)
      real(RNP), contiguous,           intent(in)    :: F_c(:,:,:,:,:)
      real(RNP), contiguous,           intent(in)    :: F_d(:,:,:,:,:)
      real(RNP), contiguous,           intent(in)    :: Q(:,:,:,:,:)
      class(BoundaryVariable_3D),      intent(in)    :: bv_u(:)
      real(RNP), contiguous, optional, intent(inout) :: mu(:,:,:,:)
      real(RNP), contiguous, optional, intent(inout) :: nu(:,:,:,:)
      real(RNP), contiguous,           intent(inout) :: u(:,:,:,:,:)

    end subroutine FGMRES_Step

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
    integer   :: i_krylov =     0 !< max num Krylov iterations
    integer   :: i_max_p  =  1000 !< max num p-iterations in projection solver
    integer   :: i_max_v  =   200 !< max num v-iterations in projection solver
    integer   :: i_pre_p  =     2 !< max num p-iterations in projection precon
    integer   :: i_pre_v  =     1 !< max num v-iterations in projection precon
    real(RNP) :: r_red    = 1e-08 !< min residual reduction, if > 0
    real(RNP) :: r_max    = 1e-12 !< max residual to reach,  if > 0
  contains
    procedure :: Bcast => Bcast_INS_TimeIntegratorOptions_3D
  end type INS_TimeIntegratorOptions_3D

contains

  !=============================================================================
  ! INS_TimeIntegrator_3D: type-bound procedures

  !-----------------------------------------------------------------------------
  !> Initialization of TimeIntegrator object

  subroutine Init_INS_TimeIntegrator_3D(this, problem, ins_op, opt)
    class(INS_TimeIntegrator_3D),        intent(inout) :: this
    class(INS_Problem_3D),       target, intent(in)    :: problem
    class(INS_Operator_3D),      target, intent(in)    :: ins_op
    class(INS_TimeIntegratorOptions_3D), intent(in)    :: opt

    this % problem => problem
    this % ins_op  => ins_op

    this % i_krylov = opt % i_krylov
    this % i_max_p  = opt % i_max_p
    this % i_max_v  = opt % i_max_v
    this % i_pre_p  = opt % i_pre_p
    this % i_pre_v  = opt % i_pre_v
    this % r_red    = opt % r_red
    this % r_max    = opt % r_max

  end subroutine Init_INS_TimeIntegrator_3D

  !=============================================================================
  ! TimeIntegratorOptions: type-bound procedures

  !-----------------------------------------------------------------------------
  !> MPI broadcasting of time-integrator options

  subroutine Bcast_INS_TimeIntegratorOptions_3D(this, root, comm)
    class(INS_TimeIntegratorOptions_3D), intent(inout) :: this
    integer,        intent(in) :: root !< rank of broadcast root
    type(MPI_Comm), intent(in) :: comm !< MPI communicator

    type(MPI_Request) :: request(7)
    integer :: n

    n = 1
    call XMPI_Ibcast( this % i_krylov, root, comm, request(n) );  n = n + 1
    call XMPI_Ibcast( this % i_max_p , root, comm, request(n) );  n = n + 1
    call XMPI_Ibcast( this % i_max_v , root, comm, request(n) );  n = n + 1
    call XMPI_Ibcast( this % i_pre_p , root, comm, request(n) );  n = n + 1
    call XMPI_Ibcast( this % i_pre_v , root, comm, request(n) );  n = n + 1
    call XMPI_Ibcast( this % r_red   , root, comm, request(n) );  n = n + 1
    call XMPI_Ibcast( this % r_max   , root, comm, request(n) )

    call MPI_Waitall( n, request, MPI_STATUSES_IGNORE )

  end subroutine Bcast_INS_TimeIntegratorOptions_3D

  !=============================================================================

end module INS__Time_Integrator__3D
