!> summary:  Spectral deferred correction method for incompressible flow
!> author:   Joerg Stiller
!> date:     2017/08/22
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module CART__ISP_Flow__SDC

  use Kind_Parameters, only: RNP
  use Constants,       only: ZERO, ONE, HALF
  use Gauss_Jacobi
  use Array_Assignments
  use XMPI
  use ISP_Flow_Problem
  use CART__DG_Weak_Gradient
  use CART__ISP_Flow__Operators
  use CART__ISP_Flow__Pressure
  use CART__ISP_Flow__Time_Derivative

  implicit none
  private

  public :: SDC_Method3D
  public :: SDC_Options3D

  !-----------------------------------------------------------------------------
  !> Spectral deferred correction parameters and procedures

  type SDC_Method3D

    procedure(Propagator), pointer, nopass :: Predictor => null()
    procedure(Propagator), pointer, nopass :: Corrector => null()

    integer   :: n_sub   = -1         !< number of subintervals
    integer   :: n_sweep = -1         !< max num correction sweeps

    real(RNP), allocatable :: xi(:)   !< GLL points
    real(RNP), allocatable :: d0(:)   !< GLL diff operator for starting time
    real(RNP), allocatable :: xs(:,:) !< subinterval GLL points
    real(RNP), allocatable :: ws(:,:) !< subinterval GLL weights

  contains

    procedure :: Init_SpectralDeferredCorrection  =>  Init_SDC
    procedure :: NumberOfSubintervals
    procedure :: IntermediateTimes
    procedure :: TimeStep

  end type SDC_Method3D

  ! constructor
  interface SDC_Method3D
    module procedure New_SDC
  end interface

  !-----------------------------------------------------------------------------
  !> Type bundling spectral deferred correction options

  type SDC_Options3D

    integer :: n_sub   = -1  !< number of subintervals
    integer :: n_sweep = -1  !< max number of correction sweeps

  contains

    procedure :: Bcast => Bcast_SDC_Options

  end type SDC_Options3D

  abstract interface

    !---------------------------------------------------------------------------
    !> Single-step time integration, optionally returning the time derivative

    subroutine Propagator(problem, flow_op, t, dt, u_0, u, F, H, S)
      import

      class(FlowProblem),   intent(in)    :: problem !< flow problem
      class(FlowOperators), intent(inout) :: flow_op !< flow operators
      real(RNP),            intent(inout) :: t       !< time
      real(RNP),            intent(in)    :: dt      !< time step size
      real(RNP),            intent(in)    :: u_0     !< solution u(t)
      real(RNP),            intent(inout) :: u       !< solution u(t+dt)
      real(RNP),  optional, intent(out)   :: F       !< ∂u/∂t(t+dt)
      real(RNP),  optional, intent(inout) :: H       !< δu/δt(t) → δu/δt(t+dt)
      real(RNP),  optional, intent(in)    :: S       !< ∫∂u/dt over (t,t+dt)

      dimension :: u_0 (:,:,:,:,:)
      dimension :: u   (:,:,:,:,:)
      dimension :: F   (:,:,:,:,:)
      dimension :: H   (:,:,:,:,:)
      dimension :: S   (:,:,:,:,:)

    end subroutine Propagator

  end interface

contains

!===============================================================================
! SDC_Options3D: type-bound procedures

subroutine Bcast_SDC_Options(this, root, comm)
  class(SDC_Options3D), intent(inout) :: this
  integer,        intent(in) :: root !< rank of broadcast root
  type(MPI_Comm), intent(in) :: comm !< MPI communicator

  type(MPI_Request)  :: request(2)
  type(MPI_Status)   :: stat(size(request))
  integer :: n

  n = 1
  call XMPI_Ibcast( this % n_sub  , root, comm, request(n) );  n = n + 1
  call XMPI_Ibcast( this % n_sweep, root, comm, request(n) )

  call MPI_Waitall(n, request, stat)

end subroutine Bcast_SDC_Options

!===============================================================================
! SDC_Method3D: type-bound procedures

!-------------------------------------------------------------------------------
!> Constructor

function New_SDC(Predictor, Corrector, opt) result(this)
  procedure(Propagator) :: Predictor  !< predictor method
  procedure(Propagator) :: Corrector  !< corrector method
  type(SDC_Options3D), intent(in) :: opt  !< SDC options

  type(SDC_Method3D) :: this

  call Init_SDC(this, Predictor, Corrector, opt)

end function New_SDC

!-------------------------------------------------------------------------------
!> Initialization of spectral deferred correction

subroutine Init_SDC(sdc, Predictor, Corrector, opt)

  ! arguments ..................................................................

  class(SDC_Method3D), intent(inout) :: sdc

  procedure(Propagator) :: Predictor  !< predictor method
  procedure(Propagator) :: Corrector  !< corrector method
  type(SDC_Options3D), intent(in) :: opt  !< SDC options

  ! local variables ............................................................

  real(RNP), allocatable :: x(:), w(:), D(:,:)
  real(RNP) :: tmp
  integer :: i, j, k, n_sub

  ! initialization .............................................................

  n_sub = opt % n_sub

  ! safeguard
  if (allocated(sdc % xs)) deallocate(sdc % xs)
  if (allocated(sdc % ws)) deallocate(sdc % ws)

  ! GLL points, weights and diff matrix
  allocate(x(0:n_sub), source = GLL_Points(n_sub))
  allocate(w(0:n_sub), source = GLL_Weights(x))
  allocate(D(0:n_sub, 0:n_sub), source = GLL_DiffMatrix(x))

  ! storage
  allocate(sdc % xs(0:n_sub, n_sub))
  allocate(sdc % ws(0:n_sub, n_sub))

  ! SDC components .............................................................

  sdc % Predictor => Predictor
  sdc % Corrector => Corrector

  sdc % n_sub = opt % n_sub

  if (opt % n_sweep >= 0) then
    sdc % n_sweep = opt % n_sweep
  else
    sdc % n_sweep = sdc % n_sub
  end if

  ! subinterval GLL points
  do j = 1, n_sub
  do i = 0, n_sub
    sdc % xs(i,j) = x(j-1) + HALF * (x(j) - x(j-1)) * (x(i) + 1)
  end do
  end do

  ! subinterval weights
  do j = 1, n_sub
  do i = 0, n_sub
    tmp = 0
    do k = 0, n_sub
      tmp = tmp + w(k) * GLL_Polynomial(i, x, sdc%xs(k,j))
    end do
    sdc % ws(i,j) = tmp
  end do
  end do

  ! GLL points
  call move_alloc(x, sdc % xi)

  ! GLL diff operator for starting time
  allocate(sdc % d0(0:n_sub), source = D(0,:))

end subroutine Init_SDC

!-------------------------------------------------------------------------------
!> Returns the number of subintervals

pure integer function NumberOfSubintervals(sdc) result(n_sub)
  class(SDC_Method3D), intent(in) :: sdc

  n_sub = sdc % n_sub

end function NumberOfSubintervals

!-------------------------------------------------------------------------------
!> Returns the intermediate times within a given time interval

pure function IntermediateTimes(sdc, t, dt) result(ti)
  class(SDC_Method3D), intent(in) :: sdc
  real(RNP), intent(in)  :: t               !< start of the time interval
  real(RNP), intent(in)  :: dt              !< length of the time interval
  real(RNP)              :: ti(0:sdc%n_sub) !< intermediate times

  real(RNP) :: c
  integer   :: m, n

  c = dt / 2
  n = sdc % n_sub

  ti(0) = t
  do m = 1, n - 1
    ti(m) = t + c * (sdc % xi(m) + 1)
  end do
  ti(n) = t + dt

end function IntermediateTimes

!-------------------------------------------------------------------------------
!> SDC time step

subroutine TimeStep( sdc, problem, flow_op, t, dt, u, F, first, last)

  ! arguments ..................................................................

  class(SDC_Method3D), intent(in) :: sdc

  class(FlowProblem),   intent(in)    :: problem !< flow problem
  class(FlowOperators), intent(inout) :: flow_op !< flow operators

  real(RNP), intent(inout) :: t                  !< time
  real(RNP), intent(in)    :: dt                 !< time step size
  real(RNP), intent(inout) :: u(:,:,:,:,:)       !< solution
  real(RNP), intent(inout) :: F(:,:,:,:,:)       !< solution time derivative
  logical,   intent(in)    :: first              !< indicates first time step
  logical,   intent(in)    :: last               !< indicates last time step

  ! local variables  ...........................................................

  real(RNP), dimension(:,:,:,:,:,:), allocatable, save :: ui, Fi, Hi, Si
  real(RNP), dimension(:),           allocatable, save :: ti, dti

  real(RNP) :: tau
  integer   :: ni
  integer   :: m, n

  ! initialization .............................................................

  ni = sdc % NumberOfSubintervals()

  call InitializeWorkspace()

  ! predictor ..................................................................

  tau = ti(0)

  do m = 1, ni

    call SetArray(ui(:,:,:,:,:,m), ui(:,:,:,:,:,m-1), multi=.true.)

    call sdc % Predictor( problem, flow_op, tau, dti(m) &
                        , u_0   = ui(:,:,:,:,:,m-1)     &
                        , u     = ui(:,:,:,:,:,m  )     &
                        , F     = Fi(:,:,:,:,:,m  )     &
                        , H     = Hi(:,:,:,:,:,m  )     )
  end do

  ! SDC iterations .............................................................

  do n = 1, sdc % n_sweep

    tau = ti(0)

    do m = 1, ni
      call SubIntegral(sdc, m, dti(m), Fi, Si(:,:,:,:,:,m))
    end do

    do m = 1, ni

      call sdc % Corrector( problem, flow_op, tau, dti(m) &
                          , u_0   = ui(:,:,:,:,:,m-1)     &
                          , u     = ui(:,:,:,:,:,m  )     &
                          , F     = Fi(:,:,:,:,:,m  )     &
                          , H     = Hi(:,:,:,:,:,m  )     &
                          , S     = Si(:,:,:,:,:,m  )     )
    end do

  end do

  ! finalization ...............................................................

  call SetArray(u, ui(:,:,:,:,:,ni), multi=.true.)
  call SetArray(F, Fi(:,:,:,:,:,ni), multi=.true.)

  !$omp single
  t = t + dt
  !$omp end single

  if (last) then
    call FreeWorkspace()
  end if

contains

  !-----------------------------------------------------------------------------
  !> Allocation and initialization of workspace

  subroutine InitializeWorkspace()

    integer :: np, ne, nc
    integer :: i

    np = size(u,1)
    ne = size(u,4)
    nc = size(u,5)

    !$omp single
    if (first) then
      allocate(  ti(0:ni                  ) )
      allocate( dti(  ni                  ) )
      allocate(  ui(  np,np,np,ne,nc,0:ni ) )
      allocate(  Fi(  np,np,np,ne,nc,0:ni ) )
      allocate(  Hi(  np,np,np,ne,nc,1:ni ) )
      allocate(  Si(  np,np,np,ne,nc,1:ni ) )
    end if
    ti  = sdc % IntermediateTimes(t, dt)
    dti = ti(1:ni) - ti(0:ni-1)
    !$omp end single

    if (first) then
      ! compute RHS without pressure contribution
      call TimeDerivative(problem, flow_op, t, u=u, F=F)
    end if

    call SetArray(ui(:,:,:,:,:,0), u, multi=.true.)
    call SetArray(Fi(:,:,:,:,:,0), F, multi=.true.)
    do i = 1, ni
      call SetArray(ui(:,:,:,:,:,i), ZERO, multi=.true.)
      call SetArray(Fi(:,:,:,:,:,i), ZERO, multi=.true.)
      call SetArray(Hi(:,:,:,:,:,i), ZERO, multi=.true.)
      call SetArray(Si(:,:,:,:,:,i), ZERO, multi=.true.)
    end do

  end subroutine InitializeWorkspace

  !-----------------------------------------------------------------------------
  !> Release of workspace

  subroutine FreeWorkspace()

    !$omp barrier
    !$omp single
    deallocate(  ti )
    deallocate( dti )
    deallocate(  ui )
    deallocate(  Fi )
    deallocate(  Hi )
    deallocate(  Si )
    !$omp end single

  end subroutine FreeWorkspace

end subroutine TimeStep

!-------------------------------------------------------------------------------
!> Evaluation of subintervals for arrays of 3D mesh variables

subroutine SubIntegral(sdc, m, dt, f, r)
  class(SDC_Method3D), intent(in) :: sdc
  integer,   intent(in)  :: m               !< interval ID, 0 < m <= sdc%n_sub
  real(RNP), intent(in)  :: dt              !< length of the time interval
  real(RNP), intent(in)  :: f(:,:,:,:,:,0:) !< integrand at intermediate times
  real(RNP), intent(out) :: r(:,:,:,:,:)    !< result

  integer :: i

  ! r = dt/2 * ws(0,m) * f(*,0)
  call MergeArrays(ZERO, r, dt/2 * sdc%ws(0,m), f(:,:,:,:,:,0), multi=.true.)

  ! r = r + dt/2 * sum(ws(1:,m) * f(*,1:))
  do i = 1, sdc%n_sub
    call MergeArrays(ONE, r, dt/2 * sdc%ws(i,m), f(:,:,:,:,:,i), multi=.true.)
  end do

end subroutine SubIntegral

!===============================================================================

end module CART__ISP_Flow__SDC
