!> summary:  Base type of SDC correctors for incompressible flows
!> author:   Joerg Stiller
!> date:     2020/05/10
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module CART__ISP_Flow__SDC_Corrector
  use, intrinsic :: ISO_Fortran_Env, only: OUTPUT_UNIT
  use Kind_Parameters, only: RNP
  use Constants,       only: HALF
  use XMPI
  use ISP_Flow_Problem
  use CART__ISP_Flow__Operators
  implicit none
  private

  public :: SDC_Corrector
  public :: SDC_CorrectorOptions

  !-----------------------------------------------------------------------------
  !> Abstract type of a SDC corrector for incompressible flow

  type, abstract :: SDC_Corrector

    class(FlowProblem),   pointer :: problem => null() !< flow problem
    class(FlowOperators), pointer :: flow_op => null() !< flow operators

    character(len=80) :: name = ''  !< time-integrator name

    ! options
    integer :: splitting  = 0 !< selector for splitting/div-stabilization
    integer :: project    = 0 !< switch for additional projections

    ! derived parameters
    real(RNP) :: chi
    real(RNP) :: cd3

  contains
    procedure, non_overridable :: Init_SDC_Corrector
    procedure, non_overridable :: Show_SDC_Corrector
    procedure :: Show => Show_SDC_Corrector
    procedure(GetCorrectorRHS), deferred :: GetCorrectorRHS
    procedure(CorrectionStep),  deferred :: CorrectionStep
  end type SDC_Corrector

  abstract interface

    !---------------------------------------------------------------------------
    !> Computes F_ex and F_im as defined in the corrector

    subroutine GetCorrectorRHS(this, t, nu_c, nu_s, u, F_ex, F_im, F)
      import

      class(SDC_Corrector), intent(in) :: this
      real(RNP), intent(in)  :: t               !< time
      real(RNP), intent(in)  :: nu_c(:,:,:,:,:) !< variable ν for corrector RHS
      real(RNP), intent(in)  :: nu_s(:,:,:,:,:) !< variable ν for subintegral RHS
      real(RNP), intent(in)  :: u   (:,:,:,:,:) !< u
      real(RNP), intent(out) :: F_ex(:,:,:,:,:) !< explicit RHS for corrector
      real(RNP), intent(out) :: F_im(:,:,:,:,:) !< implicit RHS for corrector
      real(RNP), intent(out) :: F   (:,:,:,:,:) !< RHS for subintegrals

      optional :: nu_c, nu_s

    end subroutine GetCorrectorRHS

    !---------------------------------------------------------------------------
    !> Execution of a single correction step

    subroutine CorrectionStep( this, t, dt                   &
                             , F_ex_0_old, F_ex_0, F_ex_old  &
                             , F_im_0_old, F_im_0, F_im_old  &
                             , S, u_0, u, nu                 )
      import

      class(SDC_Corrector), intent(inout) :: this
      real(RNP), intent(inout) :: t                      !< time t₀ → t
      real(RNP), intent(in)    :: dt                     !< step size ∆t = t-t₀
      real(RNP), intent(in)    :: F_ex_0_old (:,:,:,:,:) !< F^ex (t₀)ᵏ
      real(RNP), intent(in)    :: F_ex_0     (:,:,:,:,:) !< F^ex (t₀)ᵏ⁺¹
      real(RNP), intent(in)    :: F_ex_old   (:,:,:,:,:) !< F^ex (t )ᵏ
      real(RNP), intent(in)    :: F_im_0_old (:,:,:,:,:) !< F^im (t₀)ᵏ
      real(RNP), intent(in)    :: F_im_0     (:,:,:,:,:) !< F^im (t₀)ᵏ⁺¹
      real(RNP), intent(in)    :: F_im_old   (:,:,:,:,:) !< F^im (t )ᵏ
      real(RNP), intent(in)    :: S          (:,:,:,:,:) !< S    (t₀)ᵏ
      real(RNP), intent(in)    :: u_0        (:,:,:,:,:) !< u    (t₀)ᵏ⁺¹
      real(RNP), intent(inout) :: u          (:,:,:,:,:) !< u    (t )ᵏ   → (t)ᵏ⁺¹
      real(RNP), intent(inout) :: nu         (:,:,:,:,:) !< v    (t₀)ᵏ⁺¹ → ?

      optional :: nu

    end subroutine CorrectionStep

  end interface

  !-----------------------------------------------------------------------------
  !> Base type for providing time SDC-corrector options

  type SDC_CorrectorOptions
    integer :: splitting  = 2 !< selector for stage splitting scheme
    integer :: project    = 0 !< switch for additional projections
  contains
    procedure :: Bcast => Bcast_SDC_Corrector_Options
  end type SDC_CorrectorOptions

contains

  !=============================================================================
  ! SDC_Corrector: type-bound procedures

  !-----------------------------------------------------------------------------
  !> Initialization of SDC_Corrector object

  subroutine Init_SDC_Corrector(this, problem, flow_op, opt)
    class(SDC_Corrector),         intent(inout) :: this
    class(FlowProblem),   target, intent(in)    :: problem    !< flow problem
    class(FlowOperators), target, intent(in)    :: flow_op    !< flow operators
    class(SDC_CorrectorOptions), optional, intent(in) :: opt !< options

    this % problem => problem
    this % flow_op => flow_op

    ! options
    if (present(opt)) then
      this % splitting  = opt % splitting
      this % project    = opt % project
    end if

    ! derived parameters
    select case(this % splitting)
    case(1)
      this % chi = -1
      this % cd3 =  1
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

  end subroutine Init_SDC_Corrector

  !-----------------------------------------------------------------------------
  !> Output of SDC_Corrector settings

  subroutine Show_SDC_Corrector(this, unit)
    class(SDC_Corrector), intent(in) :: this
    integer,    optional, intent(in) :: unit  !< output unit

    integer :: io

    if (present(unit)) then
      io = unit
    else
      io = OUTPUT_UNIT
    end if

    write(io,'(/,A)')       'SDC_Corrector settings'
    write(io,'(A,/)')       repeat('=',80)
    write(io,'(2X,A,T15,A)')  'name:'     , trim(this % name)
    write(io,'(2X,A,T15,I0)') 'splitting:', this % splitting
    write(io,'(2X,A,T15,I0)') 'project:'  , this % project

    ! append further settings in corresponding routines of derived types

  end subroutine Show_SDC_Corrector

  !=============================================================================
  ! SDC_CorrectorOptions: type-bound procedures

  !-----------------------------------------------------------------------------
  !> MPI broadcasting of SDC-corrector options

  subroutine Bcast_SDC_Corrector_Options(this, root, comm)
    class(SDC_CorrectorOptions), intent(inout) :: this
    integer,        intent(in) :: root !< rank of broadcast root
    type(MPI_Comm), intent(in) :: comm !< MPI communicator

    type(MPI_Request) :: request(2)
    type(MPI_Status)  :: stat(size(request))
    integer :: n

    n = 1
    call XMPI_Ibcast( this % splitting, root, comm, request(n) );  n = n + 1
    call XMPI_Ibcast( this % project  , root, comm, request(n) )

    call MPI_Waitall(n, request, stat)

  end subroutine Bcast_SDC_Corrector_Options

  !=============================================================================

end module CART__ISP_Flow__SDC_Corrector
