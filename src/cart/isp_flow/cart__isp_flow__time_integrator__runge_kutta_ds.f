!> summary:  Runge-Kutta method for incompressible flows with dual splitting
!> author:   Joerg Stiller, ...
!> date:     2020/03/05
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> @todo
!>
!>   * write down algorithm for
!>       - IMEX Runge-Kutta step
!>       - generic Runge-Kutta stage, derived from dual-split IMEX Euler
!>
!>   * initialization of `class(TimeIntegrator_RungeKuttaDS)` objects
!>       - `problem` and `flow_op` as with EulerDS
!>       - `imex_rk` using `IMEX_RK_Method()` constructor  
!>
!>   * constructor for `type(TimeIntegrator_RungeKuttaDS)` objects
!>     based on initialization routine
!>
!>   * TimeStep
!>       - develop skeleton using ideas from
!>         `app/examples/1d/conv-diff/cg_convdiff_1d__imex_rk.f` and
!>         `examples/helmholtz_3d/flow_runge_kutta_method.f` @ HiBASE/first-flow
!>       - workspace
!>           - which auxiliary variables are needed
!>           - provide these in `TimeStep`
!>           - `save` attribute necessary to allow for OpenMP parallelization
!>       - note that with the considered RK schemes c(1) = 0, i.e. u1 = u
!>       - stages > 1 are similar and executed with generic procedure
!>       - think about
!>           - how F_ex and F_im are defined and
!>           - how they can be computed
!>           - this should be equivalent to the SDC procedure
!>
!>   * RungeKuttaStage
!>       - generic procedure for stages 2 - n_stage
!>       - should be derived from and hence similar to EulerDS time step
!>
!===============================================================================

module CART__ISP_Flow__Time_Integrator__Runge_Kutta_DS
  use Kind_Parameters,   only: RNP
  use Constants,         only: ONE, ZERO
  use Execution_Control, only: Error

  use Array_Assignments
  use IMEX_Runge_Kutta_Method

  use ISP_Flow_Problem
 
  use CART__ISP_Flow__Boundary_Values
  use CART__ISP_Flow__Diffusion
  use CART__ISP_Flow__Operators
  use CART__ISP_Flow__Pressure
  use CART__ISP_Flow__Projection
  use CART__ISP_Flow__Time_Derivative
  use CART__ISP_Flow__Time_Integrator

  implicit none             
  private

  public :: TimeIntegrator_RungeKuttaDS

  !-----------------------------------------------------------------------------
  !> IMEX Runge-Kutta method for incompressible flows with dual splitting

  type, extends(TimeIntegrator) :: TimeIntegrator_RungeKuttaDS
    type(IMEX_RK_Method) :: imex_rk !< IMEX Runge-Kutta method
  contains
    procedure :: Init_TimeIntegrator_RungeKuttaDS => Init_from_File
    procedure :: TimeStep
  end type TimeIntegrator_RungeKuttaDS

  ! overloading the constructor
  interface TimeIntegrator_RungeKuttaDS
    module procedure New_from_File
  end interface

contains

  !-----------------------------------------------------------------------------
  !> Constructor for objects of type TimeIntegrator_RungeKuttaDS

  function New_from_File(problem, flow_op, file) result(this)
    class(FlowProblem),   target, intent(in) :: problem !< flow problem
    class(FlowOperators), target, intent(in) :: flow_op !< flow operators
    character(len=*),             intent(in) :: file    !< input file
    type(TimeIntegrator_RungeKuttaDS) :: this
    
    call Init_from_File(this, problem, flow_op, file)

  end function New_from_File

  !-----------------------------------------------------------------------------
  !> Initialization of Init_TimeIntegrator_RungeKuttaDS object from file
  !>
  !> If the file is open, it will be read from the current position.
  !> Otherwise, it will be opened for read and closed afterwards.

  subroutine Init_from_File(this, problem, flow_op, file)
    class(TimeIntegrator_RungeKuttaDS), intent(inout) :: this
    class(FlowProblem),   target, intent(in) :: problem !< flow problem
    class(FlowOperators), target, intent(in) :: flow_op !< flow operators
    character(len=*),             intent(in) :: file    !< input file

    integer  :: ns      = 0 !< number of stages
    integer  :: method  = 1 !< RK scheme
    logical   :: opened, exists
    integer   :: io
    namelist /runge_kutta_ds_parameters/ ns, method
    
    inquire(file=file, opened=opened, exist=exists, number=io)

    if (.not. opened) then
      if (exists) then
        open(newunit=io, file=file, action='READ')
      else
        call Error('New_from_File',                                    &
                   'input file "' // trim(file) // '" does not exist', &
                   'CART__ISP_Flow__Time_Integrator'                   )
      end if
    end if

    read(io, nml=runge_kutta_ds_parameters)

    ! intialize parent type
    call this % Init_TimeIntegrator(problem, flow_op)

    ! initialize RK method
    call this % imex_rk % Init_IMEX_RK_Method(ns, method)
    ! close IO unit if file was closed on entry
    if (.not. opened) then
    close(io)
    end if

  end subroutine Init_from_File
  
  !-----------------------------------------------------------------------------
  !> Performs a single IMEX Runge-Kutta step


   ! TODO : IMEX Runge-Kutta step

  subroutine TimeStep(this, t, dt, u, nu)
    class(TimeIntegrator_RungeKuttaDS), intent(inout) :: this

    real(RNP),           intent(inout) :: t             !< time t₀ → t
    real(RNP),           intent(in)    :: dt            !< step size ∆t = t-t₀
    real(RNP), optional, intent(in)    :: nu(:,:,:,:,:) !< variable ν(x,t₀)
    real(RNP),           intent(inout) :: u (:,:,:,:,:) !< u(x,t₀) → u(x,t)

    ! local variables  .........................................................

    !type(IMEX_RK_Method)                                 :: imex
    real(RNP), dimension(:,:,:,:,:,:), allocatable, save :: K_ex, K_im !impl. & expl. K
    real(RNP), dimension(:,:,:,:,:,:), allocatable, save :: F_d, F_d1, F_d3
    real(RNP), dimension(:,:,:,:,:,:), allocatable, save :: F_s
    real(RNP), dimension(:)          , allocatable, save :: ts         ! node times
    real(RNP), dimension(:,:,:,:,:),   allocatable, save :: u_i, u_0

    ! auxiliary
    integer   :: i, k  ! i = 1,..,s
    real(RNP) :: tau    

    associate(problem => this % problem            &
             , flow_op => this % flow_op           &
             , mesh    => this % flow_op % mesh    &
             , x       => this % flow_op % x       &
             , p       => u  (:,:,:,:,4)           &
             , a_im    => this % imex_rk % a_im    &
             , a_ex    => this % imex_rk % a_ex    &
             , c       => this % imex_rk % c       &
             , b       => this % imex_rk % b       &
             , ns      => this % imex_rk % n_stage &
             , eop     => this % flow_op % eop_u  )

      
      ! workspace ................................................................
      !$omp single

      allocate(K_im ( size(u,1), size(u,2), size(u,3), &
                    & size(u,4), size(u,5), ns ))

      allocate(K_ex ( size(u,1), size(u,2), size(u,3), &
                    & size(u,4), size(u,5), ns ))

      allocate(F_d  ( size(u,1), size(u,2), size(u,3), &
                    & size(u,4), size(u,5), ns ))

      allocate(F_d1 ( size(u,1), size(u,2), size(u,3), &
                    & size(u,4), size(u,5), ns ))

      allocate(F_d3 ( size(u,1), size(u,2), size(u,3), &
                    & size(u,4), size(u,5), ns ))
      allocate(F_s  ( size(u,1), size(u,2), size(u,3), &
                    & size(u,4), size(u,5), ns ))

      allocate(u_i , mold = u)
      allocate(u_0 , mold = u)
      allocate (ts(ns))

      !$omp end single

      ! node times ...............................................................
      do k = 1, ns
        ts(k) = t + c(k) * dt 
      end do
      ! stages 1 to ns ...........................................................
      ! u₀ = u(x,t₀)
      call SetArray(u_0, u, multi=.true.)
      ! stage 1
      call TimeDerivative( problem, flow_op, ts(1)     &
                          , u_c  = u_0                 &
                          , u_d  = u_0                 &
                          , nu   = nu                  &
                          , F    = K_im (:,:,:,:,:,1)  &
                          , F_c  = K_ex (:,:,:,:,:,1)  &
                          , F_d  = F_d  (:,:,:,:,:,1)  &
                          , F_d1 = F_d1 (:,:,:,:,:,1)  &
                          , F_d3 = F_d3 (:,:,:,:,:,1)  &
                          , F_s  = F_s  (:,:,:,:,:,1)  &
                         )    
      ! K_ex(:,:,:,:,:,i) = F_c   
      ! K_im(:,:,:,:,:,i) = F - F_c 

      call MergeArrays(ONE , K_im(:,:,:,:,:,1),-ONE &
                     , K_ex(:,:,:,:,:,1), multi = .true.)

      ! stage 2,..,ns  
      do i=2, ns
        
        ! u(i) = u(t) + dt*Σ a_ex(i,:i-1) * K_ex(*,: i-1)&
        !              + dt*Σ a_im(i,:i) * K_im(*,:i)
        ! p(i) 
        call RungeKuttaStage(this , i, ts(i), dt, u_0 , &
                            nu, K_im, K_ex, F_s, F_d, F_d1, F_d3, u_i)

      end do
      
      ! compute result............................................................
      ! u(t+dt) = u(t) + dt * Σ b(i) * ( K_im(:,:ns) + K_ex(:,:ns) )

      call SetArray(u, u_i, multi = .true.)

      !do i=1, ns
        !if (a_ex(ns,i) /= ZERO) then 
        !  tau = - dt * a_ex(ns,i) 
         ! call MergeArrays(ONE, u, tau, K_ex(:,:,:,:,:,i), multi = .true. )
        !end if

        !if ( b(i) /= ZERO) then 
        !  tau = dt * b(i)
        !  call MergeArrays(ONE, u, tau, K_ex(:,:,:,:,:,i), multi = .true. )
        !end if
     ! end do
      ! prepare next timestep .....................................................
      t = t + dt

      ! clean up .................................................................

      !$omp barrier
      !$omp master
      if(allocated(K_im ))   deallocate(K_im)
      if(allocated(K_ex ))   deallocate(K_ex)
      if(allocated(F_d  ))   deallocate(F_d )
      if(allocated(F_d1 ))   deallocate(F_d1)
      if(allocated(F_d3 ))   deallocate(F_d3)
      if(allocated(F_s  ))   deallocate(F_s )
      if(allocated(u_i  ))   deallocate(u_i )
      if(allocated(u_0  ))   deallocate(u_0 )
      if(allocated(ts   ))   deallocate(ts  )
      !$omp end master
    end associate

    end subroutine TimeStep

  !-----------------------------------------------------------------------------
  !> Executes one IMEX Runge-Kutta stage

  ! TODO : generic Runge-Kutta stage, derived from dual-split IMEX Euler

  subroutine RungeKuttaStage(this, i,t, dt,u_0, nu, &           
                             K_im, K_ex, F_s, F_d, F_d1, F_d3, u_i)

    class(TimeIntegrator_RungeKuttaDS), intent(inout) :: this
    integer,              intent(in) :: i               !< stage i=2,..,s
    real(RNP),            intent(in) :: t               !< stage time tᵢ = t₀ + cᵢ∆t
    real(RNP),            intent(in) :: dt              !< time step width
    real(RNP), optional,  intent(in) :: nu  (:,:,:,:,:) !< variable ν(x,t₀)
    real(RNP),         intent(in   ) :: u_0 (:,:,:,:,:) !< u(x,t₀)
    real(RNP),         intent(inout) :: u_i (:,:,:,:,:) !< u'|u"|u(x,tᵢ)
    
    real(RNP), dimension(:,:,:,:,:,:), intent(inout) :: K_im, K_ex
    real(RNP), dimension(:,:,:,:,:,:), intent(inout) :: F_d, F_d1, F_d3, F_s

    ! local variables.........................................................
    real(RNP)              :: tau
    integer                :: j              ! j=1,..,i
    real(RNP), allocatable :: f(:,:,:,:,:)   ! rhs for 3rd step
    real(RNP), allocatable :: w(:,:,:,:,:)   ! workspace for 2nd and 3rd step
    real(RNP), allocatable :: dp(:,:,:,:)    ! workspace for final step


    associate( problem => this % problem          &
             , flow_op => this % flow_op          &
             , mesh    => this % flow_op % mesh   &
             , x       => this % flow_op % x      &
             , p       => u_i(:,:,:,:,4)          &
             , a_ex    => this % imex_rk % a_ex   &
             , a_im    => this % imex_rk % a_im   )

    ! allocate workspaces .....................................................
    !$omp single
    allocate(f   , mold = u_i)
    allocate(w   , mold = u_i)
    allocate(dp  , mold = p  )
    !$omp end single
    
    !step 1, exptrapolation: u_i→u'(x,,tᵢ) ....................................
    ! u'(i) = u(t) + dt*Σ a_ex(i,:i-1) * (K_ex(*,: i-1) + F_d(*,:i-1)) &
    !              + dt*Σ a_im(i,:i) * F_s(*,:i)

    call SetArray(u_i, u_0, multi=.true.)
    do j=1, i-1
      
      ! apply ERK scheme to F_ex and F_im
      tau = dt * a_ex(i,j)
      call MergeArrays(ONE, u_i, tau, K_ex(:,:,:,:,:,j), multi =.true.)
      call MergeArrays(ONE, u_i, tau, F_d (:,:,:,:,:,j), multi =.true.)
      
     ! upgrade F_s to DIRK scheme
      tau = dt * a_im(i,j)
      call MergeArrays(ONE, u_i, tau, F_s(:,:,:,:,:,j), multi = .true.)

    end do
    ! + for j = i
    call this % problem % GetExternalSources( x, t, F_s(:,:,:,:,:,i))
    tau = dt * a_im(i,i)
    call MergeArrays(ONE, u_i, tau, F_s(:,:,:,:,:,i), multi =.true. )

    ! boundary conditions.......................................................
    call GetBoundaryValues(problem, mesh, flow_op % bv_x, t, flow_op % bv_u)

    !step 2, projection: u_i → u"(x,tᵢ) p → p''(i)..............................
    ! solve for p''(i) Poisson eq. ∇²p''(i) = ∇ u'(i)/dt
    ! then,  u''(i) = u'(i) - dt * ∇ p''(i)
    call PressureSolver(problem, flow_op, dt, u_i, p, w)
    call ProjectionStep(problem, flow_op, dt, p, u_i, w)

    ! step 3, diffusion:  u_i → u(x,tᵢ) ........................................
    !u'''(i) = u''(i) + dt * a_im(i,i)*∇.ν(i)∇u'''(i) &
    !                 + dt * Σ(j=1,..,i-1) a_im(i,j)*F_d1(j)
    !                 - dt * Σ(j=1,..,i-1) a_ex(i,j)*(Fd_1(j)+Fd_3(j))
    
    ! f = u"(x,tᵢ) - ∆t ∑ a_ex(i,:i-1) F_13(*,:i-1)
    call SetArray(f, u_i, multi = .true.)
    do j=1, i-1
      tau = dt * a_im(i,j)
      call MergeArrays(ONE, f, tau, F_d1(:,:,:,:,:,j), multi = .true.)
      !call MergeArrays(ONE, f, tau, F_d3(:,:,:,:,:,j), multi = .true.) ! Variante 02.

      tau = -dt * a_ex(i,j) 
      call MergeArrays(ONE, f, tau, F_d1(:,:,:,:,:,j), multi = .true.)
      call MergeArrays(ONE, f, tau, F_d3(:,:,:,:,:,j), multi = .true.)
    end do

    ! solve diffusion equation
    tau = dt * a_im(i,i)
    call DiffusionStep(problem, flow_op, tau, f=f, u=u_i, w=w, nu=nu)

    ! step4, final projection in case the viscosity is variable.................
    
    if (flow_op % control % div_final) then

      ! solve for p = p" + dp
      call SetArray(dp, ZERO)
      call PressureSolver(problem, flow_op, dt, u_i, dp, w) ! u=u(x,t₀+∆t)
      call MergeArrays(ONE, p, ONE, dp)

      ! v = v''' - 1/∆t ∇p - J(v)
      call ProjectionStep(problem, flow_op, dt, dp, u_i, w)

    end if

    ! now determine the values of K_im, K_ex, F_d, F_d1, F_d3 and F_s......

    call TimeDerivative( problem, flow_op, t         &
                        , u_c  = u_i                 &
                        , u_d  = u_i                 &
                        , nu   = nu                  &
                        , F    = K_im (:,:,:,:,:,i)  &
                        , F_c  = K_ex (:,:,:,:,:,i)  &
                        , F_d  = F_d  (:,:,:,:,:,i)  &
                        , F_d1 = F_d1 (:,:,:,:,:,i)  &
                        , F_d3 = F_d3 (:,:,:,:,:,i)  &
                        , F_s  = F_s  (:,:,:,:,:,i)  &
                        )
    ! K_im
    call MergeArrays(ONE, K_im(:,:,:,:,:,i), -ONE, K_ex(:,:,:,:,:,i), multi=.true.)

    !clean up .............................................................
    !!$omp barrier
    !!$omp master
    !$omp single
    if (allocated(dp  )) deallocate(dp )
    if (allocated(w   )) deallocate(w  )
    if (allocated(f   )) deallocate(f  )
    !$omp end single
    !!$omp end master
    end associate
    
 end subroutine RungeKuttaStage

  !=============================================================================

end module CART__ISP_Flow__Time_Integrator__Runge_Kutta_DS

