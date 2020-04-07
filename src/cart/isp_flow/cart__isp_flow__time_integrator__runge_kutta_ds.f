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
  use Kind_Parameters, only: RNP
  use Constants,       only: ONE, ZERO
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
    type(IMEX_RK_Method) :: imex_rk !< imex_rk using IMEX_RK_Method() constructor
  contains
    procedure :: Init_TimeIntegrator_RungeKuttaDS
    procedure :: TimeStep
  end type TimeIntegrator_RungeKuttaDS

  ! overloading the constructor
  interface TimeIntegrator_RungeKuttaDS
    module procedure New_TimeIntegrator_RungeKuttaDS
  end interface
 
contains

  !-----------------------------------------------------------------------------
  !> Constructor for objects of type TimeIntegrator_RungeKuttaDS

  function New_TimeIntegrator_RungeKuttaDS(problem, flow_op) result(this)
    class(FlowProblem),   target, intent(in) :: problem !< flow problem
    class(FlowOperators), target, intent(in) :: flow_op !< flow operators
    type(TimeIntegrator_RungeKuttaDS) :: this

    call Init_TimeIntegrator_RungeKuttaDS(this, problem, flow_op)

  end function New_TimeIntegrator_RungeKuttaDS

  !-----------------------------------------------------------------------------
  !> Initialization of Init_TimeIntegrator_RungeKuttaDS object

  subroutine Init_TimeIntegrator_RungeKuttaDS(this, problem, flow_op)
    class(TimeIntegrator_RungeKuttaDS), intent(inout) :: this
    class(FlowProblem),   target, intent(in) :: problem !< flow problem
    class(FlowOperators), target, intent(in) :: flow_op !< flow operators

    call this % Init_TimeIntegrator(problem, flow_op)

  end subroutine Init_TimeIntegrator_RungeKuttaDS

  
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
    real(RNP), dimension(:,:,:,:,:,:), allocatable       :: F_s
    real(RNP), dimension(:)          , allocatable, save :: ts         ! node times
    real(RNP), dimension(:,:,:,:,:),   allocatable, save :: u_i, u_0

    ! auxiliary
    integer   :: ns ! number of stages
    integer   :: i, k  ! i = 1,..,s
    real(RNP) :: tau    
    associate(problem => this % problem           &
             , flow_op => this % flow_op          &
             , mesh    => this % flow_op % mesh   &
             , x       => this % flow_op % x      &
             , p       => u  (:,:,:,:,4)          &
             , a_im    => this % imex_rk % a_im   &
             , a_ex    => this % imex_rk % a_ex   &
             , c       => this % imex_rk % c      &
             , b       => this % imex_rk % b      &
             , eop     => this % flow_op % eop_u  )

      

      ! number of stages s........................................................
      !imex = this % imex_rk
      ns   = this % imex_rk % n_stage                           !?
      if ( t < dt ) then                                        !? hier
        call this % imex_rk % Init_IMEX_RK_Method(ns)           !? nicht
      end if                                                    

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
        ts(k) = t + this % imex_rk % c(k) * dt
      end do

      
      ! initialization   .........................................................
      
      F_s  = ZERO
      F_d  = ZERO
      F_d1 = ZERO
      F_d3 = ZERO
      K_im = ZERO
      K_ex = ZERO

      ! compute result............................................................
      ! u(t+dt) = u(t) + dt * Σ b(i) * ( K_im(:,:ns) + K_ex(:,:ns) )

      ! u₀ = u(x,t₀)
      call SetArray(u_0, u, multi=.true.)
      !print*, u(1,2,2,2,4)  !>> test
      ! stage 1
        call TimeDerivative( problem, flow_op, ts(1)    &
                            , u_c  = u                  &
                            , u_d  = u                  &
                            , p    = p                  &
                            , nu   = nu                 &
                            , F    = K_im(:,:,:,:,:,1)  &
                            , F_c  = K_ex(:,:,:,:,:,1)  &
                            , F_d  = F_d (:,:,:,:,:,1)  &
                            , F_d1 = F_d1(:,:,:,:,:,1)  &
                            , F_d3 = F_d3(:,:,:,:,:,1)  &
                            , F_s  = F_s (:,:,:,:,:,1)  &
                           )
    
        ! K_ex(:,:,:,:,:,i) = F_c   

        ! K_im(:,:,:,:,:,i) = F - F_c 

        call MergeArrays(ONE , K_im(:,:,:,:,:,1),-ONE &
                       , K_ex(:,:,:,:,:,1), multi = .true.)

     if ( this % imex_rk % b(1) /= ZERO )  then
        !u = u +  dt * imex % b(1) * (K_im(:,:,:,:,:,1) + K_ex(:,:,:,:,:,1))
        tau = dt * this % imex_rk % b(1)
        call MergeArrays(ONE,u, tau, K_im(:,:,:,:,:,1), multi = .true.)
        call MergeArrays(ONE,u, tau, K_ex(:,:,:,:,:,1), multi = .true.)
      end if
      

      ! stage 2 - ns
      ! determine different stage values
      
      do i=2, ns
        

        call RungeKuttaStage(this , i, ts(i), dt, u_0 , &
                            nu, K_im, K_ex, F_s, F_d, F_d1, F_d3, u_i)
     
        !print*,'u_i', u_i(1,2,2,2,4)  !>> test

        !u = u +  dt * b(i) * (K_ex(*,i) + K_im(*,i))
        tau = dt * this % imex_rk % b(i)

        call MergeArrays(ONE, u, tau, K_ex(:,:,:,:,:,i), multi = .true. )

        call MergeArrays(ONE, u, tau, K_im(:,:,:,:,:,i), multi = .true. )
        if (i == ns) then 
          call SetArray(u(:,:,:,:,4), u_i(:,:,:,:,4), multi = .true.)
        end if
        !print*,'u', u(1,2,2,2,4)  !>> test
      end do
      
      ! prepare next timestep
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

    print*,'p=', p(1,2,2,2)
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
    !print*, p(1,2,2,2)
    !print*, u_i(1,2,2,2,4)
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
    
    call TimeDerivative( problem, flow_op, t        &
                        , u_c  = u_i                &
                        , u_d  = u_i                &
                        , p    = p                  &
                        , nu   = nu                 &
                        , F    = K_im (:,:,:,:,:,i) &
                        , F_c  = K_ex (:,:,:,:,:,i) &
                        , F_d  = F_d  (:,:,:,:,:,i) &
                        , F_d1 = F_d1(:,:,:,:,:,i)  &
                        , F_d3 = F_d3(:,:,:,:,:,i)  &
                        , F_s  = F_s (:,:,:,:,:,i)  &
                        )
    
    ! K_im
    call MergeArrays(ONE, K_im(:,:,:,:,:,i), -ONE, K_ex(:,:,:,:,:,i), multi=.true.)
    
    ! just a print test ...................................................
    !print*,'K_ex=', K_ex(1,1,1,1,1,2)       
    !print*,'K_im=', K_im(1,1,1,1,1,2)   
    !print*,'F_d=', F_d (1,1,1,1,1,2) 
    !print*,'F_d1=', F_d1(1,1,1,1,1,2)
    !print*,'F_d3=', F_d3( 1,1,1,1,1,2)   
    !print*,'F_s=', F_s( 1,1,1,1,1,2) 

    !clean up .............................................................
    !$omp barrier
    !$omp master
    if (allocated(dp  )) deallocate(dp )
    if (allocated(w   )) deallocate(w  )
    if (allocated(f   )) deallocate(f  )
    !$omp end master
    end associate
 end subroutine RungeKuttaStage

  !=============================================================================

end module CART__ISP_Flow__Time_Integrator__Runge_Kutta_DS

