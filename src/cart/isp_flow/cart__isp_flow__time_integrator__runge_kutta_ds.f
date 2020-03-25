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

    type(IMEX_RK_Method)                                 :: imex          
    real(RNP), dimension(:,:,:,:,:,:), allocatable, save :: K_ex, K_im !impl. & expl. K
    real(RNP), dimension(:,:,:,:,:,:), allocatable, save :: F_d1, F_d3
    real(RNP), dimension(:,:,:,:,:,:), allocatable, save :: F_s
    real(RNP), dimension(:,:,:,:,:)  , allocatable, save :: F, F_c 
    real(RNP), dimension(:)          , allocatable, save :: ts         ! node times
    real(RNP), dimension(:,:,:,:,:),   allocatable, save :: u_i, u_0  
 
    ! auxiliary
    integer   :: ns ! number of stages
    integer   :: i  ! i = 1,..,s
    
    associate(problem => this % problem           &
             , flow_op => this % flow_op          &
             , mesh    => this % flow_op % mesh   &
             , x       => this % flow_op % x      &
             , eop     => this % flow_op % eop_u  &
             , p       => u(:,:,:,:,4) )
      
      ! initialization   .........................................................
      
      ! number of stages s 
      imex = this % imex_rk
      ns   = imex % n_stage

      ! workspace ................................................................
      !$omp single
      allocate(K_im ( size(u,1), size(u,2), size(u,3), &
                    & size(u,4), size(u,5), ns )) 

      allocate(K_ex ( size(u,1), size(u,2), size(u,3), &
                    & size(u,4), size(u,5), ns ))

      allocate(F_d1 ( size(u,1), size(u,2), size(u,3), &
                    & size(u,4), size(u,5), ns ))

      allocate(F_d3 ( size(u,1), size(u,2), size(u,3), &
                    & size(u,4), size(u,5), ns ))
      allocate(F_s  ( size(u,1), size(u,2), size(u,3), &
                    & size(u,4), size(u,5), ns ))
      allocate(F   , mold = u)
      allocate(F_c , mold = u)
      allocate(u_i , mold = u)
      allocate(u_0 , mold = u)
      allocate (ts(ns))
      !$omp end single 
      
      ! node times ...............................................................
      do i=1, ns
        ts(i) = t + imex % c(i) * dt 
      end do  

      ! boundary conditions ......................................................
      call GetBoundaryValues(problem, mesh, flow_op%bv_x, t, flow_op%bv_u)
 
      ! steps ....................................................................
      
      ! u' = u₀ ≡ u(x,t₀)
      call SetArray(u_0,u,multi=.true.)

      ! stage 1 
      if (any(imex % a_im(:,1) /= 0).or. imex % b(1) /= 0) then  
        call TimeDerivative( problem, flow_op, ts(1) &
                            , u_c  = u_0             &
                            , u_d  = u_0             &
                            , nu   = nu              &
                            , F    = F               &
                            , F_c  = F_c             &
                            ) 
      
        K_im(:,:,:,:,:,1) = F - F_c
        K_ex(:,:,:,:,:,1) = F_c
        !u = u +  dt * imex % b(1) * (K_im(:,:,:,:,:,1) + K_ex(:,:,:,:,:,1))
        call MergeArrays(ONE,u,dt * imex % b(1), K_im(:,:,:,:,:,1),multi=.true.)
        call MergeArrays(ONE,u,dt * imex % b(1), K_ex(:,:,:,:,:,1),multi=.true.)
      end if 

      ! stage 2 - ns
      ! determine different stage values
      ! compute result
      ! u(t+dt) = u(t) + dt * Σ(i=1,..,ns) b(i) * (K_im(:,i)+K_ex(:,i))
      
      do i=2, ns
         
        call RungeKuttaStage(this, imex, i, dt, u_0 , nu,    &
                             K_im, K_ex, F_s, F_d1, F_d3, u_i)
        ! stage values 
        call TimeDerivative( problem, flow_op, ts(i) &
                          , u_c  = u_i               &
                          , u_d  = u_i               &
                          , nu   = nu                &
                          , F    = F                 &
                          , F_c  = F_c               &
                          , F_d1 = F_d1(:,:,:,:,:,i) &
                          , F_d3 = F_d3(:,:,:,:,:,i) &
                          , F_s  = F_s (:,:,:,:,:,i) &
                          ) 
        ! explicit and implicit values K_im, K_ex
        K_im(:,:,:,:,:,i) = F - F_c
        K_ex(:,:,:,:,:,i) = F_c
        u = u +  dt * imex % b(i) * (K_im(:,:,:,:,:,i) + K_ex(:,:,:,:,:,i))

      end do 
      
      ! prepare next timestep 
      t = t + dt

      ! clean up .................................................................
    
      !$omp barrier
      !$omp master
      if(allocated(K_im ))   deallocate(K_im) 
      if(allocated(K_ex ))   deallocate(K_ex) 
      if(allocated(F_d1 ))   deallocate(F_d1) 
      if(allocated(F_d3 ))   deallocate(F_d3)
      if(allocated(F_s  ))   deallocate(F_s ) 
      if(allocated(u_i  ))   deallocate(u_i )
      if(allocated(u_0  ))   deallocate(u_0 )
      if(allocated(F    ))   deallocate(F   )
      if(allocated(F_c  ))   deallocate(F_c )
      if(allocated(ts   ))   deallocate(ts  )
      !$omp end master
    end associate
    end subroutine TimeStep

  !-----------------------------------------------------------------------------
  !> Executes one IMEX Runge-Kutta stage
  
  ! TODO : generic Runge-Kutta stage, derived from dual-split IMEX Euler

  subroutine RungeKuttaStage(this, imex, i, dt,u_n, nu, &
                             K_im, K_ex, F_s, F_d1, F_d3, u_i)

    class(TimeIntegrator_RungeKuttaDS), intent(inout) :: this
    
    type(IMEX_RK_Method), intent(in) :: imex            !< IMEX butcher tableau
    integer,              intent(in) :: i               !< stage i=2,..,s  
    real(RNP),            intent(in) :: dt              !< time step width 
    real(RNP), optional,  intent(in) :: nu(:,:,:,:,:)   !< variable ν(x,t₀) 
    real(RNP),         intent(in   ) :: u_n (:,:,:,:,:) !< u(x,t₀) 
    real(RNP),         intent(inout) :: u_i (:,:,:,:,:) !< u(x,t₀)
    
    real(RNP), dimension(:,:,:,:,:,:), intent(inout) :: K_im, K_ex, F_d1, F_d3, F_s
 
    ! local variables...................................................
    integer                :: j              ! j=1,..,i
    real(RNP), allocatable :: w(:,:,:,:,:)  ! ....
    real(RNP), allocatable :: u_s(:,:,:,:,:)  ! ....
    real(RNP), allocatable :: p(:,:,:,:)
    ! allocate workspaces ..............................................
    allocate(w , mold = u_i) 
    allocate(u_s , mold = u_i) 
    allocate(p(size(u_i,1), size(u_i,2), size(u_i,3), &
                    & size(u_i,4)))
    
    !step 1: exptrapolation of u'(i) ...................................
    call SetArray(u_i, u_n, multi=.true.)

    do j=1, i-1 

      call MergeArrays(ONE, u_i, dt * imex % a_ex(i,j),K_ex(:,:,:,:,:,j), multi=.true.)

      call MergeArrays(ONE, u_i, dt * imex % a_ex(i,j),K_im(:,:,:,:,:,j), multi=.true.)

      call MergeArrays(ONE, u_i, dt * (imex % a_ex(i,j) - imex % a_im(i,j)) &
                                ,K_ex(:,:,:,:,:,j), multi=.true.)
 
    end do
    ! + for p = i 
    !u_i = u_i + dt*imex % a_im(i,i) * F_s(:,:,:,:,:,i)
    call MergeArrays(ONE, u_i,dt * imex % a_im(i,i), F_s(:,:,:,:,:,i), multi=.true. )

    !step 2: projection => p''(i), ∇ p''(i) and u''(i) = u'(i) - dt * ∇ p''(i)
    call PressureSolver(this % problem, this % flow_op, dt, u_i, p, w)
    call ProjectionStep(this % problem, this % flow_op, dt, p, u_i, w)


    !step 3: diffusion: => u'''(i) = u''(i) + dt * a_im(i,i)*∇.ν(i)∇u'''(i) &
    !                                + dt * Σ(j=1,..,i-1) a_im(i,j)*F_d1(j)  
    !                                - dt * Σ(j=1,..,i-1) a_ex(i,j)*(Fd_1(j)+Fd_3(j))
    do j=1, i-1
      call MergeArrays(ONE, u_i, -dt * imex % a_ex(i,j),&
                          F_d1(:,:,:,:,:,i), multi=.true.)

      call MergeArrays(ONE, u_i, -dt * imex % a_ex(i,j),&
                          F_d3(:,:,:,:,:,i), multi=.true.)

      call MergeArrays(ONE, u_i,  dt * imex % a_im(i,j),& 
                          F_d1(:,:,:,:,:,i), multi=.true.)
    end do 
    call DiffusionStep(this % problem, this % flow_op, dt* imex % a_im(i,i),&
                       f=u_i, u=u_s, w=w, nu=nu)
    call SetArray(u_i,u_s, multi=.true.)
    call SetArray(u_i(:,:,:,:,4),p, multi=.true.)
    !step 4: in case the viscosity is variable execute final projection ...
     
    !clean up ............................................................
      !$omp barrier
      !$omp master
      if (allocated(p   )) deallocate(p   )
      if (allocated(w   )) deallocate(w   )
      if (allocated(u_s )) deallocate(u_s )
      !$omp end master
  end subroutine RungeKuttaStage

  !=============================================================================

end module CART__ISP_Flow__Time_Integrator__Runge_Kutta_DS
