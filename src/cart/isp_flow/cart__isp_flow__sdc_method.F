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

  use CART__ISP_Flow__SDC_Corrector
  use CART__ISP_Flow__SDC_Corrector__Euler

  implicit none
  private

  public :: SDC_FlowOptions
  public :: SDC_FlowMethod

  !-----------------------------------------------------------------------------
  !> Type for providing SDC options

  type, extends(SDC_Options) :: SDC_FlowOptions
    integer :: pressure  =  0  !< switch for pressure recomputation
  contains
    procedure :: Bcast => Bcast_SDC_FlowOptions
  end type SDC_FlowOptions

  !-----------------------------------------------------------------------------
  !> Modular SDC method with flexible choice of predictor and corrector

  type, extends(SDC_Method) :: SDC_FlowMethod

    integer :: pressure  =  0  !< switch for pressure recomputation

    class(FlowProblem),    pointer     :: problem => null() !< flow problem
    class(FlowOperators),  pointer     :: flow_op => null() !< flow operators
    class(TimeIntegrator), allocatable :: predictor         !< predictor method
    class(SDC_Corrector),  allocatable :: corrector         !< corrector method

  contains

    procedure :: Init_SDC_FlowMethod
    procedure :: Show => Show_SDC_FlowMethod
    procedure :: TimeStep

  end type SDC_FlowMethod

  ! overloading the constructor
  interface SDC_FlowMethod
    module procedure New_SDC_FlowMethod
  end interface

  !-----------------------------------------------------------------------------
  ! Constants

  real(RNP), parameter :: chi = -ONE  !< scaling parameter of F_d3 = χ∇(ν∇⋅v)

contains

  !=============================================================================
  ! SDC_Option: type-bound procedures

  !-----------------------------------------------------------------------------
  !> Broadcasting the IDP Flow SDC options

  subroutine Bcast_SDC_FlowOptions(this, root, comm)
    class(SDC_FlowOptions), intent(inout) :: this
    integer,        intent(in) :: root !< rank of broadcast root
    type(MPI_Comm), intent(in) :: comm !< MPI communicator

    call this % SDC_Options % Bcast(root, comm)

    call XMPI_Bcast(this % pressure, root, comm)

  end subroutine Bcast_SDC_FlowOptions

  !=============================================================================
  ! SDC_FlowMethod: type-bound procedures

  !-----------------------------------------------------------------------------
  !> Constructor for objects of type SDC_FlowMethod

  function New_SDC_FlowMethod(problem, flow_op, pre_opt, cor_opt, sdc_opt) &
      result(this)

    class(FlowProblem),   target, intent(in) :: problem !< flow problem
    class(FlowOperators), target, intent(in) :: flow_op !< flow operators
    class(TimeIntegratorOptions), intent(in) :: pre_opt !< predictor options
    class(SDC_CorrectorOptions),  intent(in) :: cor_opt !< corrector method
    class(SDC_FlowOptions),       intent(in) :: sdc_opt !< SDC options
    type(SDC_FlowMethod) :: this

    call Init_SDC_FlowMethod(this, problem, flow_op, pre_opt, cor_opt, sdc_opt)

  end function New_SDC_FlowMethod

  !-----------------------------------------------------------------------------
  !> Initialization of SDC_FlowMethod object

  subroutine Init_SDC_FlowMethod( this, problem, flow_op, pre_opt, cor_opt, &
                                  sdc_opt                                   )

    ! arguments ................................................................

    class(SDC_FlowMethod), intent(inout) :: this

    class(FlowProblem),   target, intent(in) :: problem !< flow problem
    class(FlowOperators), target, intent(in) :: flow_op !< flow operators
    class(TimeIntegratorOptions), intent(in) :: pre_opt !< predictor options
    class(SDC_CorrectorOptions),  intent(in) :: cor_opt !< corrector method
    class(SDC_FlowOptions),       intent(in) :: sdc_opt !< SDC options

    ! options ..................................................................

    ! initialize base type
    call this % Init_SDC_Method(sdc_opt)

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

    ! corrector ................................................................

    select type(cor_opt)
    class is (SDC_CorrectorOptions_Euler)
      this % corrector = SDC_Corrector_Euler(problem, flow_op, cor_opt)
    end select

  end subroutine Init_SDC_FlowMethod

  !-----------------------------------------------------------------------------
  !> Output of SDC_FlowMethod settings

  subroutine Show_SDC_FlowMethod(this, unit)
    class(SDC_FlowMethod), intent(in) :: this
    integer, optional, intent(in) :: unit  !< output unit

    integer :: io

    if (present(unit)) then
      io = unit
    else
      io = OUTPUT_UNIT
    end if

    call this % SDC_Method % Show(unit)

    write(io,'(2X,A,T15,I0)') 'pressure:' , this % pressure

    call this % predictor % Show(unit)
    call this % corrector % Show(unit)

  end subroutine Show_SDC_FlowMethod

  !-----------------------------------------------------------------------------
  !> SDC time step

  subroutine TimeStep(this, t, dt, u, standby)

    ! arguments ................................................................

    class(SDC_FlowMethod), intent(inout) :: this         !< (flow_op may change)
    real(RNP),             intent(inout) :: t            !< time
    real(RNP),             intent(in)    :: dt           !< time step size
    real(RNP),             intent(inout) :: u(:,:,:,:,:) !< solution
    logical,     optional, intent(in)    :: standby      !< reuse workspace [F]

    ! local variables  .........................................................

    real(RNP), dimension(:),           allocatable, save :: t_, dt_
    real(RNP), dimension(:,:,:,:,:,:), allocatable, save :: u_
    real(RNP), dimension(:,:,:,:,:,:), allocatable, save :: F_ex_, F_im_
    real(RNP), dimension(:,:,:,:,:,:), allocatable, save :: F_, S_
    real(RNP), dimension(:,:,:,:,:),   allocatable, save :: F_ex_0_old
    real(RNP), dimension(:,:,:,:,:),   allocatable, save :: F_im_0_old

    real(RNP), allocatable, target, save :: nu_(:,:,:,:,:,:)
    real(RNP), contiguous, pointer, save :: nu_0(:,:,:,:,:) => null()
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

    !$omp single

    if (allocated(u_)) then
      if (any(shape(u_) /= [np,np,np,ne,nc,n_sub])) then
        deallocate(t_, dt_)
        deallocate(u_, F_ex_, F_im_, F_, S_)
        deallocate(F_ex_0_old, F_im_0_old)
        if (allocated(nu_)) deallocate(nu_)
        nu_0 => null()
        nu_i => null()
      end if
    end if

    if (.not. allocated(u_)) then

      allocate(t_  (0:n_sub))
      allocate(dt_ (1:n_sub))

      allocate(u_         (np,np,np,ne,nc,0:n_sub))
      allocate(F_ex_      (np,np,np,ne,nc,0:n_sub))
      allocate(F_im_      (np,np,np,ne,nc,0:n_sub))
      allocate(F_         (np,np,np,ne,nc,0:n_sub))
      allocate(S_         (np,np,np,ne,nc,1:n_sub))
      allocate(F_ex_0_old (np,np,np,ne,nc))
      allocate(F_im_0_old (np,np,np,ne,nc))

      if (this % problem % HasVariableProperties()) then
        allocate(nu_(np,np,np,ne,nc,0:n_sub))
      end if

    end if

    t_ (0:n_sub) = this % IntermediateTimes(t, dt)
    dt_(1:n_sub) = t_(1:n_sub) - t_(0:n_sub-1)

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
          !$omp single
          nu_0 => nu_(:,:,:,:,:,max(i-1,0))
          nu_i => nu_(:,:,:,:,:,i)
          !$omp end single
          call this % problem %                  &
                 GetDiffusivity( this%flow_op%x  &
                               , t_(i)           &
                               , u_(:,:,:,:,:,i) &
                               , nu_i            )
        end if

        ! RHS for corrector and subintegrals
        call this % corrector %                           &
               GetCorrectorRHS( t_(i)                     &
                              , nu_c = nu_0               &
                              , nu_s = nu_i               &
                              , u    = u_   (:,:,:,:,:,i) &
                              , F_ex = F_ex_(:,:,:,:,:,i) &
                              , F_im = F_im_(:,:,:,:,:,i) &
                              , F    = F_   (:,:,:,:,:,i) )
      end do

      ! correction sweeps ......................................................

      Sweeps: do n = 1, n_sweep

        ! subintegrals
        do i = 1, n_sub
          call SubIntegral(this, i, dt, F_, S_(:,:,:,:,:,i))
        end do

        call SetArray(F_ex_0_old, F_ex_(:,:,:,:,:,0), multi=.true.)
        call SetArray(F_im_0_old, F_im_(:,:,:,:,:,0), multi=.true.)

        do i = 1, n_sub

          t_0 = t_(i-1)

          ! initialize diffusivity with v(tᵢ₋₁)ᵏ⁺¹
          if (allocated(nu_)) then
            !$omp single
            nu_0 => nu_(:,:,:,:,:,i-1)
            nu_i => nu_(:,:,:,:,:,i)
            !$omp end single
            call SetArray(nu_i, nu_0)
          end if

          ! correction
          call this % corrector %                                  &
                 CorrectionStep( t          = t_0                  &
                               , dt         = dt_(i)               &
                               , F_ex_0_old = F_ex_0_old           &
                               , F_ex_0     = F_ex_(:,:,:,:,:,i-1) &
                               , F_ex_old   = F_ex_(:,:,:,:,:,i)   &
                               , F_im_0_old = F_im_0_old           &
                               , F_im_0     = F_im_(:,:,:,:,:,i-1) &
                               , F_im_old   = F_im_(:,:,:,:,:,i)   &
                               , S          = S_   (:,:,:,:,:,i)   &
                               , u_0        = u_   (:,:,:,:,:,i-1) &
                               , u          = u_   (:,:,:,:,:,i)   &
                               , nu         = nu_i                 )

          ! save old RHS
          if (i < n_sub) then
            call SetArray(F_ex_0_old, F_ex_(:,:,:,:,:,i), multi=.true.)
            call SetArray(F_im_0_old, F_im_(:,:,:,:,:,i), multi=.true.)
          end if

          ! update variable diffusivity
          if (associated(nu_i)) then
            call this % problem %                  &
                   GetDiffusivity( this%flow_op%x  &
                                 , t_(i)           &
                                 , u_(:,:,:,:,:,i) &
                                 , nu_i            )
          end if

          ! update RHS
          call this % corrector %                           &
                 GetCorrectorRHS( t_(i)                     &
                                , nu_c = nu_0               &
                                , nu_s = nu_i               &
                                , u    = u_   (:,:,:,:,:,i) &
                                , F_ex = F_ex_(:,:,:,:,:,i) &
                                , F_im = F_im_(:,:,:,:,:,i) &
                                , F    = F_   (:,:,:,:,:,i) )

        end do

      end do Sweeps

    end if Corrector

    ! result ...................................................................

    !$omp single
    t = t + dt
    !$omp end single

    call SetArray(u, u_(:,:,:,:,:,n_sub), multi=.true.)

    if (this % pressure > 0) then
      associate( F_v => u_(:,:,:,:,:,0) &
               , w   => u_(:,:,:,:,:,1) &
               , p   => u (:,:,:,:,4)   )
        call TimeDerivative( this % problem, this % flow_op, t, u         &
                           , nu = nu_i, chi = chi, F = F_v )
        call ComputePressure( this % problem, this % flow_op, t, F_v, p, w )
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
    if (allocated(u_)) then
      deallocate(t_, dt_)
      deallocate(u_, F_ex_, F_im_, F_, S_)
      deallocate(F_ex_0_old, F_im_0_old)
      if (allocated(nu_)) then
        deallocate(nu_)
      end if
      nu_0 => null()
      nu_i => null()
    end if
    !$omp end single

  end subroutine TimeStep

  !=============================================================================
  ! SDC_FlowMethod: auxiliary procedures

  !-----------------------------------------------------------------------------
  !> Computes the integrals of the time derivative F over subintervals [tᵢ₋₁,tᵢ]

  subroutine SubIntegral(sdc, m, dt, F, S)
    class(SDC_FlowMethod), intent(in) :: sdc
    integer,   intent(in)  :: m               !< interval ID, 0 < m <= sdc%n_sub
    real(RNP), intent(in)  :: dt              !< length of the full time step
    real(RNP), intent(in)  :: F(:,:,:,:,:,0:) !< integrand at intermediate times
    real(RNP), intent(out) :: S(:,:,:,:,:)    !< result

    integer :: i

    ! S = 0
    call SetArray(S, ZERO, multi=.true.)

    ! S = dt * sum(ws(:,m) * F(*,:))
    do i = 0, sdc%n_sub
      call MergeArrays(ONE, S, dt * sdc%ws(i,m), F(:,:,:,:,:,i), multi=.true.)
    end do

  end subroutine SubIntegral

  !=============================================================================

end module CART__ISP_Flow__SDC_Method
