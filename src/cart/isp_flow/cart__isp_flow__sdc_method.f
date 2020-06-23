!> summary:  SDC method for incompressible flows
!> author:   Joerg Stiller
!> date:     2020/05/14
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module CART__ISP_Flow__SDC_Method

  use, intrinsic :: ISO_Fortran_Env, only: OUTPUT_UNIT

  use Kind_Parameters, only: RNP
  use Constants,       only: ZERO, ONE, HALF
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

  use CART__ISP_Flow__SDC_Corrector
  use CART__ISP_Flow__SDC_Corrector__Euler

  implicit none
  private

  public :: SDC_Options
  public :: SDC_Method

  !-----------------------------------------------------------------------------
  !> Type for providing SDC options

  type SDC_Options
    integer :: n_sub    =  1      !< number of subintervals
    integer :: n_sweep  = -1      !< max number of correction sweeps
    integer :: pressure =  0      !< switch for pressure recomputation
  contains
    procedure :: Bcast => Bcast_SDC_Options
  end type SDC_Options

  !-----------------------------------------------------------------------------
  !> Modular SDC method with flexible choice of predictor and corrector

  type SDC_Method

    integer :: n_sub    = -1          !< number of subintervals
    integer :: n_sweep  = -1          !< max num correction sweeps
    integer :: pressure =  0          !< switch for pressure recomputation

    real(RNP), allocatable :: xi(:)   !< SDC points in [-1,1]
    real(RNP), allocatable :: ws(:,:) !< subinterval quadrature weights

    class(FlowProblem),    pointer     :: problem => null() !< flow problem
    class(FlowOperators),  pointer     :: flow_op => null() !< flow operators
    class(TimeIntegrator), allocatable :: predictor         !< predictor method
    class(SDC_Corrector),  allocatable :: corrector         !< corrector method

  contains

    procedure :: Init_SDC_Method
    procedure :: Show => Show_SDC_Method
    procedure :: NumberOfSubintervals
    procedure :: IntermediateTimes
    procedure :: TimeStep

  end type SDC_Method

  ! overloading the constructor
  interface SDC_Method
    module procedure New_SDC_Method
  end interface

contains

  !=============================================================================
  ! SDC_Option: type-bound procedures

  subroutine Bcast_SDC_Options(this, root, comm)
    class(SDC_Options), intent(inout) :: this
    integer,        intent(in) :: root !< rank of broadcast root
    type(MPI_Comm), intent(in) :: comm !< MPI communicator

    type(MPI_Request)  :: request(3)
    type(MPI_Status)   :: stat(size(request))
    integer :: n

    n = 1
    call XMPI_Ibcast( this % n_sub   , root, comm, request(n) );  n = n + 1
    call XMPI_Ibcast( this % n_sweep , root, comm, request(n) );  n = n + 1
    call XMPI_Ibcast( this % pressure, root, comm, request(n) )

    call MPI_Waitall(n, request, stat)

  end subroutine Bcast_SDC_Options

  !=============================================================================
  ! SDC_Method: type-bound procedures

  !-----------------------------------------------------------------------------
  !> Constructor for objects of type SDC_Method

  function New_SDC_Method(problem, flow_op, pre_opt, cor_opt, sdc_opt) &
      result(this)

    class(FlowProblem),    target, intent(in) :: problem !< flow problem
    class(FlowOperators),  target, intent(in) :: flow_op !< flow operators
    class(TimeIntegratorOptions),  intent(in) :: pre_opt !< predictor options
    class(SDC_Corrector_Options),  intent(in) :: cor_opt !< corrector method
    class(SDC_Options),            intent(in) :: sdc_opt !< SDC options
    type(SDC_Method) :: this

    call Init_SDC_Method(this, problem, flow_op, pre_opt, cor_opt, sdc_opt)

  end function New_SDC_Method

  !-----------------------------------------------------------------------------
  !> Initialization of SDC_Method object

  subroutine Init_SDC_Method(this, problem, flow_op, pre_opt, cor_opt, sdc_opt)

    ! arguments ................................................................

    class(SDC_Method), intent(inout) :: this

    class(FlowProblem),    target, intent(in) :: problem !< flow problem
    class(FlowOperators),  target, intent(in) :: flow_op !< flow operators
    class(TimeIntegratorOptions),  intent(in) :: pre_opt !< predictor options
    class(SDC_Corrector_Options),  intent(in) :: cor_opt !< corrector method
    class(SDC_Options),            intent(in) :: sdc_opt !< SDC options

    ! local variables ..........................................................

    real(RNP), allocatable :: x (:), w (:)
    real(RNP), allocatable :: xs(:), ws(:,:)
    integer :: i, j, k, n_sub

    ! prerequisites ............................................................

    n_sub = max(1, sdc_opt % n_sub)

    ! GLL points and weights
    allocate(x(0:n_sub), source = GLL_Points(n_sub))
    allocate(w(0:n_sub), source = GLL_Weights(x))

    ! workspace
    allocate(xs(0:n_sub), ws(0:n_sub, n_sub))

    ! subintervals and sweeps ..................................................

    ! number of subintervals
    this % n_sub = sdc_opt % n_sub

    ! max number of correction sweeps
    if (sdc_opt % n_sweep >= 0) then
      this % n_sweep = sdc_opt % n_sweep
    else
      this % n_sweep = 2 * n_sub - 1
    end if

    ! pressure recomputation
    this % pressure = sdc_opt % pressure

    ! points and weights .......................................................

    ! subinterval weights
    do j = 1, n_sub

      ! scaled GLL points
      do i = 0, n_sub
        xs(i) = x(j-1) + HALF * (x(j) - x(j-1)) * (x(i) + 1)
      end do

      ! scaled GLL weights
      do i = 0, n_sub
        ws(i,j) = 0
        do k = 0, n_sub
          ws(i,j) = ws(i,j) + w(k) * GLL_Polynomial(i, x, xs(k))
        end do
      end do

    end do

    call move_alloc(x , this % xi)
    call move_alloc(ws, this % ws)

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

    ! corrector ................................................................

    select type(cor_opt)
    class is (SDC_Corrector_Euler_Options)
      this % corrector = SDC_Corrector_Euler(problem, flow_op, cor_opt)
    end select

  end subroutine Init_SDC_Method

  !-----------------------------------------------------------------------------
  !> Output of TimeIntegrator_Euler settings

  subroutine Show_SDC_Method(this, unit)
    class(SDC_Method), intent(in) :: this
    integer, optional, intent(in) :: unit  !< output unit

    integer :: io

    if (present(unit)) then
      io = unit
    else
      io = OUTPUT_UNIT
    end if

    write(io,'(/,A)')       'SDC_Method settings'
    write(io,'(A,/)')       repeat('≡',80)
    write(io,'(2X,A,T15,I0)') 'n_sub'     , this % n_sub
    write(io,'(2X,A,T15,I0)') 'n_sweeps:' , this % n_sweep

    call this % predictor % Show(unit)
    call this % corrector % Show(unit)

  end subroutine Show_SDC_Method

  !-----------------------------------------------------------------------------
  !> Returns the number of subintervals

  pure integer function NumberOfSubintervals(this) result(n_sub)
    class(SDC_Method), intent(in) :: this

    n_sub = this % n_sub

  end function NumberOfSubintervals

  !-----------------------------------------------------------------------------
  !> Returns the intermediate times within a given time interval

  pure function IntermediateTimes(this, t, dt) result(ti)
    class(SDC_Method), intent(in) :: this
    real(RNP), intent(in)  :: t                !< start of the time interval
    real(RNP), intent(in)  :: dt               !< length of the time interval
    real(RNP)              :: ti(0:this%n_sub) !< intermediate times

    real(RNP) :: c
    integer   :: m, n

    c = dt / 2
    n = this % n_sub

    ti(0) = t
    do m = 1, n - 1
      ti(m) = t + c * (this % xi(m) + 1)
    end do
    ti(n) = t + dt

  end function IntermediateTimes

  !-----------------------------------------------------------------------------
  !> SDC time step

  subroutine TimeStep(this, t, dt, u, standby)

    ! arguments ................................................................

    class(SDC_Method), intent(inout) :: this         !< (flow_op may change)
    real(RNP),         intent(inout) :: t            !< time
    real(RNP),         intent(in)    :: dt           !< time step size
    real(RNP),         intent(inout) :: u(:,:,:,:,:) !< solution
    logical, optional, intent(in)    :: standby      !< reuse workspace [F]

    ! local variables  .........................................................

    real(RNP), dimension(:),           allocatable, save :: t_, dt_
    real(RNP), dimension(:,:,:,:,:,:), allocatable, save :: u_
    real(RNP), dimension(:,:,:,:,:,:), allocatable, save :: F_ex_, F_im_, F_d3_
    real(RNP), dimension(:,:,:,:,:,:), allocatable, save :: F_, S_
    real(RNP), dimension(:,:,:,:,:),   allocatable, save :: F_ex_0_old
    real(RNP), dimension(:,:,:,:,:),   allocatable, save :: F_im_0_old
    real(RNP), dimension(:,:,:,:,:),   allocatable, save :: F_d3_0_old

    real(RNP), allocatable, target, save :: nu_(:,:,:,:,:,:)
    real(RNP), contiguous, pointer, save :: nu_i(:,:,:,:,:) => null()

    real(RNP) :: t_0
    integer   :: np, ne, nc, n_sub, n_sweep
    integer   :: i, n

    !---------------------------------------------------------------------------
    ! initialization

    np = size(u,1)
    ne = size(u,4)
    nc = size(u,5)

    n_sub   = this % n_sub
    n_sweep = this % n_sweep

    !$omp barrier
    !$omp single

    if (allocated(u_)) then
      if (any(shape(u_) /= [np,np,np,ne,nc,n_sub])) then
        deallocate(t_, dt_)
        deallocate(u_, F_ex_, F_im_, F_d3_, F_, S_)
        deallocate(F_ex_0_old, F_im_0_old, F_d3_0_old)
        if (allocated(nu_)) deallocate(nu_)
      end if
    end if

    if (.not. allocated(u_)) then

      allocate(t_  (0:n_sub))
      allocate(dt_ (1:n_sub))

      allocate(u_         (np,np,np,ne,nc,0:n_sub))
      allocate(F_ex_      (np,np,np,ne,nc,0:n_sub))
      allocate(F_im_      (np,np,np,ne,nc,0:n_sub))
      allocate(F_d3_      (np,np,np,ne,nc,0:n_sub))
      allocate(F_         (np,np,np,ne,nc,0:n_sub))
      allocate(S_         (np,np,np,ne,nc,1:n_sub))
      allocate(F_ex_0_old (np,np,np,ne,nc))
      allocate(F_im_0_old (np,np,np,ne,nc))
      allocate(F_d3_0_old (np,np,np,ne,nc))
      !$acc enter data create(t_,dt_,u_,...)

      if (this % problem % HasVariableProperties()) then
        allocate(nu_(np,np,np,ne,nc,0:n_sub))
        !$acc enter data create(nu_)
      end if

    end if

    t_  = this % IntermediateTimes(t, dt)
    dt_ = t_(1:n_sub) - t_(0:n_sub-1)

    !$omp end single

    !---------------------------------------------------------------------------
    ! predictor

    ! u⁰(t_0) = u(t_0)
    call SetArray(u_(:,:,:,:,:,0), u(:,:,:,:,:), multi=.true.)

    do i = 1, n_sub
      t_0 = t_(i-1)
      call SetArray(u_(:,:,:,:,:,i), u_(:,:,:,:,:,i-1), multi=.true.)
      call this % predictor % TimeStep(t_0, dt_(i), u_(:,:,:,:,:,i))
    end do

    !---------------------------------------------------------------------------
    ! corrector

    Corrector: if (n_sweep > 0) then

      ! prerequisites ..........................................................

      do i = 0, n_sub

        ! initialize variable diffusivity
        if (allocated(nu_)) then
          nu_i => nu_(:,:,:,:,:,i)
          call this % problem %                  &
                 GetDiffusivity( this%flow_op%x  &
                               , t_(i)           &
                               , u_(:,:,:,:,:,i) &
                               , nu_i            )
        end if

        ! RHS for corrector and subintegrals
        call this % corrector %                           &
               GetCorrectorRHS( t_(i)                     &
                              , nu_i                      &
                              , u    = u_   (:,:,:,:,:,i) &
                              , F_ex = F_ex_(:,:,:,:,:,i) &
                              , F_im = F_im_(:,:,:,:,:,i) &
                              , F_d3 = F_d3_(:,:,:,:,:,i) &
                              , F    = F_   (:,:,:,:,:,i) )

      end do

      ! correction sweeps ......................................................

      Sweeps: do n = 1, n_sweep

        ! subintegrals
        do i = 1, n_sub
          call SubIntegral(this, i, dt_(i), F_, S_(:,:,:,:,:,i))
        end do

        call SetArray(F_ex_0_old, F_ex_(:,:,:,:,:,0), multi=.true.)
        call SetArray(F_im_0_old, F_im_(:,:,:,:,:,0), multi=.true.)
        call SetArray(F_d3_0_old, F_d3_(:,:,:,:,:,0), multi=.true.)

        do i = 1, n_sub

          t_0 = t_(i-1)

          ! correction
          call this % corrector %                                  &
                 CorrectionStep( t          = t_0                  &
                               , dt         = dt_(i)               &
                               , F_ex_0_old = F_ex_0_old           &
                               , F_ex_0     = F_ex_(:,:,:,:,:,i-1) &
                               , F_ex       = F_ex_(:,:,:,:,:,i)   &
                               , F_im_0_old = F_im_0_old           &
                               , F_im_0     = F_im_(:,:,:,:,:,i-1) &
                               , F_im       = F_im_(:,:,:,:,:,i)   &
                               , F_d3_0_old = F_d3_0_old           &
                               , F_d3_0     = F_d3_(:,:,:,:,:,i-1) &
                               , S          = S_   (:,:,:,:,:,i)   &
                               , u_0        = u_   (:,:,:,:,:,i-1) &
                               , u          = u_   (:,:,:,:,:,i)   &
                               , nu         = nu_i                 )

          ! save old RHS
          if (i < n_sub) then
            call SetArray(F_ex_0_old, F_ex_(:,:,:,:,:,i), multi=.true.)
            call SetArray(F_im_0_old, F_im_(:,:,:,:,:,i), multi=.true.)
            call SetArray(F_d3_0_old, F_d3_(:,:,:,:,:,i), multi=.true.)
          end if

          ! update RHS
          call this % corrector %                           &
                 GetCorrectorRHS( t_(i)                     &
                                , nu_i                      &
                                , u    = u_   (:,:,:,:,:,i) &
                                , F_ex = F_ex_(:,:,:,:,:,i) &
                                , F_im = F_im_(:,:,:,:,:,i) &
                                , F_d3 = F_d3_(:,:,:,:,:,i) &
                                , F    = F_   (:,:,:,:,:,i) )

          ! update variable diffusivity
          if (associated(nu_i)) then
            call this % problem %                  &
                   GetDiffusivity( this%flow_op%x  &
                                 , t_(i)           &
                                 , u_(:,:,:,:,:,i) &
                                 , nu_i            )
          end if

        end do

      end do Sweeps

    end if Corrector

    ! result ...................................................................

    t = t + dt
    call SetArray(u, u_(:,:,:,:,:,n_sub), multi=.true.)

    if (this % pressure > 0) then
      associate( F_v => u_(:,:,:,:,:,0) &
               , w   => u_(:,:,:,:,:,1) &
               , p   => u (:,:,:,:,4)   )
        call TimeDerivative( this % problem, this % flow_op, t, u, u      &
                           , nu = nu_i, chi = this%corrector%chi, F = F_v )
        call PressureSolver( this % problem, this % flow_op, F_v, t, p, w )
      end associate
    end if

    ! clean-up .................................................................

    ! keep workspace in case of standby
    if (present(standby)) then
      if (standby) return
    end if

    !$omp barrier
    !$omp single
    if (allocated(u_)) then
      !$acc exit data delete(...)
      deallocate(t_, dt_)
      deallocate(u_, F_ex_, F_im_, F_d3_, F_, S_)
      deallocate(F_ex_0_old, F_im_0_old, F_d3_0_old)
      if (allocated(nu_)) then
        deallocate(nu_)
      end if
      nu_i => null()
    end if
    !$omp end single

  end subroutine TimeStep

  !=============================================================================
  ! SDC_Method: auxiliary procedures

  !-----------------------------------------------------------------------------
  !> Evaluation of subinterval integrals for arrays of 3D mesh variables

  subroutine SubIntegral(sdc, m, dt, f, r)
    class(SDC_Method), intent(in) :: sdc
    integer,   intent(in)  :: m               !< interval ID, 0 < m <= sdc%n_sub
    real(RNP), intent(in)  :: dt              !< length of the time interval
    real(RNP), intent(in)  :: f(:,:,:,:,:,0:) !< integrand at intermediate times
    real(RNP), intent(out) :: r(:,:,:,:,:)    !< result

    integer :: i

    ! r = 0
    call SetArray(r, ZERO, multi=.true.)

    ! r = dt/2 * sum(ws(:,m) * f(*,:))
    do i = 0, sdc%n_sub
      call MergeArrays(ONE, r, dt/2 * sdc%ws(i,m), f(:,:,:,:,:,i), multi=.true.)
    end do

  end subroutine SubIntegral

  !=============================================================================

end module CART__ISP_Flow__SDC_Method
