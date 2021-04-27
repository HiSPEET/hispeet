!> summary:  SDC method for incompressible flows
!> author:   Joerg Stiller
!> date:     2020/05/14
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module CART__ISP_Flow__SDC_Method

  use, intrinsic :: ISO_Fortran_Env, only: OUTPUT_UNIT

  use Kind_Parameters, only: RNP
  use Constants,       only: ZERO, ONE, HALF
  use Spectral_Deferred_Correction
  use Gauss_Jacobi
  use Array_Assignments
  use XMPI

  use ISP_Flow_Problem

  use CART__ISP_Flow__Operators
  use CART__ISP_Flow__Pressure
  use CART__ISP_Flow__Time_Derivative

  use CART__ISP_Flow__Time_Integrator
  use CART__ISP_Flow__Time_Integrator__Euler
  use CART__ISP_Flow__Time_Integrator__Runge_Kutta

  implicit none
  private

  public :: Flow_SDC_Options
  public :: Flow_SDC_Method

  !-----------------------------------------------------------------------------
  !> Type for providing SDC options
  !>
  !> Note
  !>   - `project` and `pressure` do not affect corresponding predictor settings

  type, extends(SDC_Options) :: Flow_SDC_Options
    integer :: project  = 0 !< additional projections: 0 none, 1 final, 2 substep
    integer :: pressure = 0 !< pressure recomputation: 0 none, 1 final, 2 substep
  contains
    procedure :: Bcast => Bcast_Flow_SDC_Options
  end type Flow_SDC_Options

  !-----------------------------------------------------------------------------
  !> Modular SDC method with flexible choice of predictor and corrector

  type, abstract, extends(SDC_Method) :: Flow_SDC_Method

    class(FlowProblem),    pointer     :: problem => null() !< flow problem
    class(FlowOperators),  pointer     :: flow_op => null() !< flow operators
    class(TimeIntegrator), allocatable :: predictor         !< predictor method

    ! corrector options
    character(len=80) :: corrector_name = ''
    integer :: project  !< additional projections: 0 none, 1 final, 2 substep
    integer :: pressure !< pressure recomputation: 0 none, 1 final, 2 substep

  contains

    procedure :: TimeStep

    procedure :: Init_Flow_SDC_Method !< common initialization
    procedure :: Show_Flow_SDC_Method !< common settings output

    procedure(CorrectorRHS ), deferred, private :: CorrectorRHS
    procedure(CorrectorStep), deferred, private :: CorrectorStep

  end type Flow_SDC_Method

  abstract interface

    !---------------------------------------------------------------------------
    !> Computes RHS at the SDC nodes
    !>
    !> Input
    !>
    !>   - `t_sdc   (  0:n_sub)`    SDC time nodes
    !>   - `nu_sdc  (…,0:n_sub)`    variable diffusivity
    !>   - `u_sdc   (…,0:n_sub)`    solution variables
    !>   - `m0                 `    time node to start with [0]
    !>
    !> Output
    !>
    !>   - `F_sdc    (…,m0:n_sub)`  RHS evaluated at SDC time nodes
    !>   - `F_sdc_ex (…,m0:n_sub)`  explicit RHS of low-order scheme
    !>   - `F_sdc_im (…,m0:n_sub)`  semi-/implicit RHS of low-order scheme
    !>
    !> Notes
    !>
    !>   - in case of constant diffusivity omit `nu_sdc`
    !>   - the low-order RHS is computed only if `F_lo_ex` and `F_lo_im`
    !>     are both present
    !>   - the optional argument `m0` can be used to avoid the unnecessary
    !>     recomputation of the RHS at `t_sdc(0)`

    subroutine CorrectorRHS( this, t_sdc, nu_sdc, u_sdc &
                           , F_sdc, F_sdc_ex, F_sdc_im  &
                           , m0                         )
      import

      class(Flow_SDC_Method),      intent(in)    :: this
      real(RNP),                   intent(in)    :: t_sdc    (0:)
      real(RNP), target, optional, intent(in)    :: nu_sdc   (:,:,:,:,:,0:)
      real(RNP),                   intent(in)    :: u_sdc    (:,:,:,:,:,0:)
      real(RNP),                   intent(inout) :: F_sdc    (:,:,:,:,:,0:)
      real(RNP),         optional, intent(inout) :: F_sdc_ex (:,:,:,:,:,0:)
      real(RNP),         optional, intent(inout) :: F_sdc_im (:,:,:,:,:,0:)
      integer  ,         optional, intent(in)    :: m0

    end subroutine CorrectorRHS

    !---------------------------------------------------------------------------
    !> Execution of a single correction step
    !>
    !> Assuming that the correction k+1 is performed on subinterval `m`, the
    !> meaning of the arguments is as follows:
    !>
    !> Input
    !>
    !>   - `m`                         SDC interval to process
    !>   - `t_sdc        (  0:n_sub)`  SDC time nodes
    !>   - `nu_sdc       (…,0:m-1)  `  diffusivity   after  k+1  corrections
    !>   - `nu_sdc       (…,m:n_sub)`  diffusivity   after  k    corrections
    !>   - `u_sdc        (…,0:m-1)  `  solution      after  k+1  corrections
    !>   - `u_sdc        (…,m:n_sub)`  solution      after  k    corrections
    !>   - `F_sdc        (…,0:n_sub)`  RHS           after  k    corrections
    !>   - `F_sdc_ex     (…,0:n_sub)`  explicit RHS  after  k    corrections
    !>   - `F_sdc_im     (…,0:n_sub)`  implicit RHS  after  k    corrections
    !>   - `F_sdc_ex_new (…,0:m-1)  `  explicit RHS  after  k+1  corrections
    !>   - `F_sdc_im_new (…,0:m-1)  `  implicit RHS  after  k+1  corrections
    !>
    !> Output
    !>
    !>   - `nu_sdc       (…,m)`        diffusivity   after  k+1  corrections
    !>   - `u_sdc        (…,m)`        solution      after  k+1  corrections
    !>   - `F_sdc_ex_new (…,m)`        explicit RHS  after  k+1  corrections
    !>   - `F_sdc_im_new (…,m)`        implicit RHS  after  k+1  corrections

    subroutine CorrectorStep( this, m, t_sdc, nu_sdc, u_sdc  &
                            , F_sdc, F_sdc_ex, F_sdc_im      &
                            , F_sdc_ex_new, F_sdc_im_new     )
      import

      class(Flow_SDC_Method),      intent(inout) :: this
      integer  ,                   intent(in)    :: m
      real(RNP),                   intent(in)    :: t_sdc        (0:)
      real(RNP), target, optional, intent(inout) :: nu_sdc       (:,:,:,:,:,0:)
      real(RNP),                   intent(inout) :: u_sdc        (:,:,:,:,:,0:)
      real(RNP),                   intent(in)    :: F_sdc        (:,:,:,:,:,0:)
      real(RNP),                   intent(in)    :: F_sdc_ex     (:,:,:,:,:,0:)
      real(RNP),                   intent(in)    :: F_sdc_im     (:,:,:,:,:,0:)
      real(RNP),                   intent(inout) :: F_sdc_ex_new (:,:,:,:,:,0:)
      real(RNP),                   intent(inout) :: F_sdc_im_new (:,:,:,:,:,0:)

    end subroutine CorrectorStep

  end interface

  !-----------------------------------------------------------------------------
  ! Constants

  real(RNP), parameter :: chi = -ONE  !< scaling parameter of F_d3 = χ∇(ν∇⋅v)

contains

  !=============================================================================
  ! SDC_Option: type-bound procedures

  !-----------------------------------------------------------------------------
  !> Broadcasting the IDP Flow SDC options

  subroutine Bcast_Flow_SDC_Options(this, root, comm)
    class(Flow_SDC_Options), intent(inout) :: this
    integer,        intent(in) :: root !< rank of broadcast root
    type(MPI_Comm), intent(in) :: comm !< MPI communicator

    call this % SDC_Options % Bcast(root, comm)

    call XMPI_Bcast(this % project , root, comm)
    call XMPI_Bcast(this % pressure, root, comm)

  end subroutine Bcast_Flow_SDC_Options

  !=============================================================================
  ! Flow_SDC_Method: type-bound procedures

  !-----------------------------------------------------------------------------
  !> Initialization of Flow_SDC_Method object

  subroutine Init_Flow_SDC_Method(this, problem, flow_op, pre_opt, sdc_opt)

    ! arguments ................................................................

    class(Flow_SDC_Method), intent(inout) :: this

    class(FlowProblem),   target, intent(in) :: problem !< flow problem
    class(FlowOperators), target, intent(in) :: flow_op !< flow operators
    class(TimeIntegratorOptions), intent(in) :: pre_opt !< predictor options
    class(Flow_SDC_Options),      intent(in) :: sdc_opt !< SDC options

    ! options ..................................................................

    ! initialize base type
    call this % Init_SDC_Method(sdc_opt)

    this % project  = sdc_opt % project
    this % pressure = sdc_opt % pressure

    ! problem and flow operators ...............................................

    this % problem => problem
    this % flow_op => flow_op

    ! predictor ................................................................

    select type(pre_opt)
    class is (TimeIntegrator_Euler_Options)
      this % predictor = TimeIntegrator_Euler(problem, flow_op, pre_opt)
    class is (TimeIntegrator_RungeKutta_Options)
      this % predictor = TimeIntegrator_RungeKutta(problem, flow_op, pre_opt)
    end select

  end subroutine Init_Flow_SDC_Method

  !-----------------------------------------------------------------------------
  !> Output of Flow_SDC_Method settings

  subroutine Show_Flow_SDC_Method(this, unit)
    class(Flow_SDC_Method), intent(in) :: this
    integer, optional, intent(in) :: unit  !< output unit

    integer :: io

    if (present(unit)) then
      io = unit
    else
      io = OUTPUT_UNIT
    end if

    call this % SDC_Method % Show(unit)

    call this % predictor % Show(unit)

    write(io,'(/,A)') 'Corrector settings'
    write(io,'(A,/)') repeat('=',80)
    write(io,'(2X,A,T15,A)')  'name:'     , trim(this % corrector_name)
    write(io,'(2X,A,T15,I0)') 'project:'  , this % project
    write(io,'(2X,A,T15,I0)') 'pressure:' , this % pressure

  end subroutine Show_Flow_SDC_Method

  !-----------------------------------------------------------------------------
  !> SDC time step

  subroutine TimeStep(this, t, dt, u, standby)

    ! arguments ................................................................

    class(Flow_SDC_Method), intent(inout) :: this         !< (flow_op may change)
    real(RNP),             intent(inout) :: t            !< time
    real(RNP),             intent(in)    :: dt           !< time step size
    real(RNP),             intent(inout) :: u(:,:,:,:,:) !< solution
    logical,     optional, intent(in)    :: standby      !< reuse workspace [F]

    ! local variables  .........................................................

    real(RNP), dimension(:),           allocatable, save :: t_sdc, dt_sdc
    real(RNP), dimension(:,:,:,:,:,:), allocatable, save :: u_sdc, nu_sdc
    real(RNP), dimension(:,:,:,:,:,:), allocatable, save :: F_sdc
    real(RNP), dimension(:,:,:,:,:,:), allocatable, save :: F_sdc_ex
    real(RNP), dimension(:,:,:,:,:,:), allocatable, save :: F_sdc_im
    real(RNP), dimension(:,:,:,:,:,:), allocatable, save :: F_sdc_ex_new
    real(RNP), dimension(:,:,:,:,:,:), allocatable, save :: F_sdc_im_new

    real(RNP) :: t_0
    integer   :: np, ne, nc, n_sub, n_sweep
    integer   :: m, n

    !---------------------------------------------------------------------------
    ! initialization

    np = size(u,1)
    ne = size(u,4)
    nc = size(u,5)

    n_sub   = this % n_sub
    n_sweep = this % n_sweep

    !$omp single

    if (allocated(u_sdc)) then
      if (any(shape(u_sdc) /= [np,np,np,ne,nc,n_sub+1])) then
        deallocate(t_sdc, dt_sdc)
        deallocate(u_sdc, F_sdc, F_sdc_ex, F_sdc_im, F_sdc_ex_new, F_sdc_im_new)
        if (allocated(nu_sdc)) deallocate(nu_sdc)
      end if
    end if

    if (.not. allocated(u_sdc)) then

      allocate(t_sdc  (0:n_sub))
      allocate(dt_sdc (1:n_sub))

      allocate(u_sdc        (np,np,np,ne,nc,0:n_sub))
      allocate(F_sdc        (np,np,np,ne,nc,0:n_sub))
      allocate(F_sdc_ex     (np,np,np,ne,nc,0:n_sub))
      allocate(F_sdc_im     (np,np,np,ne,nc,0:n_sub))
      allocate(F_sdc_ex_new (np,np,np,ne,nc,0:n_sub))
      allocate(F_sdc_im_new (np,np,np,ne,nc,0:n_sub))

      if (this % problem % HasVariableProperties()) then ! or SGS model ?
        allocate(nu_sdc(np,np,np,ne,nc,0:n_sub))
      end if

    end if

    t_sdc (0:n_sub) = this % IntermediateTimes(t, dt)
    dt_sdc(1:n_sub) = t_sdc(1:n_sub) - t_sdc(0:n_sub-1)

    !$omp end single

    !---------------------------------------------------------------------------
    ! predictor

    ! u⁰(t_0) = u(t_0)
    call SetArray(u_sdc(:,:,:,:,:,0), u(:,:,:,:,:), multi=.true.)

    do m = 1, n_sub
      t_0 = t_sdc(m-1)
      call SetArray(u_sdc(:,:,:,:,:,m), u_sdc(:,:,:,:,:,m-1), multi=.true.)
      call this % predictor % TimeStep(t_0, dt_sdc(m), u_sdc(:,:,:,:,:,m))
    end do

    !---------------------------------------------------------------------------
    ! corrector

    Corrector: if (n_sweep > 0) then

      ! prerequisites ..........................................................

      do m = 0, n_sub

        ! initialize variable diffusivity
        if (allocated(nu_sdc)) then
          if (this % problem % HasVariableProperties()) then

            call this % problem %                      &
                   GetDiffusivity( this%flow_op%x      &
                                 , t_sdc (m)           &
                                 , u_sdc (:,:,:,:,:,m) &
                                 , nu_sdc(:,:,:,:,:,m) )
          end if

        end if

        call this % CorrectorRHS(t_sdc, nu_sdc, u_sdc, F_sdc, F_sdc_ex, F_sdc_im)

      end do

      ! these will never change ;)
      call SetArray(F_sdc_ex_new(:,:,:,:,:,0), F_sdc_ex(:,:,:,:,:,0), multi=.true.)
      call SetArray(F_sdc_im_new(:,:,:,:,:,0), F_sdc_im(:,:,:,:,:,0), multi=.true.)

      ! correction sweeps ......................................................

      Sweeps: do n = 1, n_sweep

        do m = 1, n_sub
          call this % CorrectorStep( m, t_sdc, nu_sdc, u_sdc    &
                                   , F_sdc, F_sdc_ex, F_sdc_im  &
                                   , F_sdc_ex_new, F_sdc_im_new )
        end do

        if (n == n_sweep) exit

        do m = 1, n_sub
          call SetArray(F_sdc_ex(:,:,:,:,:,m), F_sdc_ex_new(:,:,:,:,:,m), multi=.true.)
          call SetArray(F_sdc_im(:,:,:,:,:,m), F_sdc_im_new(:,:,:,:,:,m), multi=.true.)
        end do

        call this % CorrectorRHS(t_sdc, nu_sdc, u_sdc, F_sdc, m0=1)

      end do Sweeps

    end if Corrector

    ! result ...................................................................

    !$omp single
    t = t + dt
    !$omp end single

    call SetArray(u, u_sdc(:,:,:,:,:,n_sub), multi=.true.)

    ! optional pressure recomputation
    if (this % pressure == 1) then
      associate(p => u(:,:,:,:,4), F_v => F_sdc(:,:,:,:,:,0))
        if (allocated(nu_sdc)) then
          call TimeDerivative( this % problem, this % flow_op, t, u    &
                             , nu = nu_sdc(:,:,:,:,:,n_sub), chi = chi &
                             , F = F_v                                 )
        else
          call TimeDerivative( this % problem, this % flow_op, t, u    &
                             , chi = chi, F = F_v                      )
        end if
        call ComputePressure( this % problem, this % flow_op, t, F_v, p &
                            , w = F_sdc(:,:,:,:,:,1))
      end associate
    end if

    !---------------------------------------------------------------------------
    ! clean-up

    ! keep workspace in case of standby
    if (present(standby)) then
      if (standby) return
    end if

    !$omp barrier
    !$omp single
    if (allocated(u_sdc)) then
      deallocate(t_sdc, dt_sdc, u_sdc)
      deallocate(F_sdc, F_sdc_ex, F_sdc_im, F_sdc_ex_new, F_sdc_im_new)
      if (allocated(nu_sdc)) then
        deallocate(nu_sdc)
      end if
    end if
    !$omp end single

  end subroutine TimeStep

  !=============================================================================

end module CART__ISP_Flow__SDC_Method
