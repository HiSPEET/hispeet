!> summary:  Euler method for incompressible flows with dual splitting
!> author:   Joerg Stiller
!> date:     2020/03/05
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module CART__ISP_Flow__Time_Integrator__Euler
  use Kind_Parameters, only: RNP
  use Constants,       only: ZERO, ONE
  use Array_Assignments

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

  public :: TimeIntegrator_Euler
  public :: TimeIntegrator_Euler_Options

  !-----------------------------------------------------------------------------
  !> IMEX Euler method for incompressible flows with dual splitting

  type, extends(TimeIntegrator) :: TimeIntegrator_Euler
  contains
    procedure :: Init_TimeIntegrator_Euler
    procedure :: Show => Show_TimeIntegrator_Euler
    procedure :: TimeStep
  end type TimeIntegrator_Euler

  ! overloading the constructor
  interface TimeIntegrator_Euler
    module procedure New_TimeIntegrator_Euler
  end interface

  !-----------------------------------------------------------------------------
  !> Type for providing Euler time-integrator options (none, so far)
  !>
  !> Options inherited from base class
  !>
  !>   * `project`                                                            \n
  !>        0: no additional projection step                                  \n
  !>        1: additional projection at the end of the time step
  !>
  !>   * `pressure`
  !>        0: return pressure as computed                                    \n
  !>        1: recompute pressure at the end of the time step

  type, extends(TimeIntegratorOptions) :: TimeIntegrator_Euler_Options
  end type TimeIntegrator_Euler_Options

contains

  !-----------------------------------------------------------------------------
  !> Constructor for objects of type TimeIntegrator_Euler with options

  function New_TimeIntegrator_Euler(problem, flow_op, opt) result(this)
    class(FlowProblem),                            intent(in) :: problem
    class(FlowOperators),                          intent(in) :: flow_op
    class(TimeIntegrator_Euler_Options), optional, intent(in) :: opt
    type(TimeIntegrator_Euler) :: this

    call Init_TimeIntegrator_Euler(this, problem, flow_op, opt)

  end function New_TimeIntegrator_Euler

  !-----------------------------------------------------------------------------
  !> Initialization of a Init_TimeIntegrator_Euler object

  subroutine Init_TimeIntegrator_Euler(this, problem, flow_op, opt)
    class(TimeIntegrator_Euler),                   intent(inout) :: this
    class(FlowProblem),                            intent(in)    :: problem
    class(FlowOperators),                          intent(in)    :: flow_op
    class(TimeIntegrator_Euler_Options), optional, intent(in)    :: opt

    ! intialize parent type
    call this % Init_TimeIntegrator(problem, flow_op, opt)
    this % name = 'IMEX Euler method'

  end subroutine Init_TimeIntegrator_Euler

  !-----------------------------------------------------------------------------
  !> Output of TimeIntegrator_Euler settings

  subroutine Show_TimeIntegrator_Euler(this, unit)
    class(TimeIntegrator_Euler), intent(in) :: this
    integer,           optional, intent(in) :: unit  !< output unit

    ! show parent settings
    call this % Show_TimeIntegrator(unit)

  end subroutine Show_TimeIntegrator_Euler

  !-----------------------------------------------------------------------------
  !> Performs a single IMEX Euler step

  subroutine TimeStep(this, t, dt, u, standby)
    class(TimeIntegrator_Euler), intent(inout) :: this
    real(RNP),         intent(inout) :: t             !< time t₀ → t
    real(RNP),         intent(in)    :: dt            !< step size ∆t = t-t₀
    real(RNP),         intent(inout) :: u (:,:,:,:,:) !< u(x,t₀) → u(x,t)
    logical, optional, intent(in)    :: standby       !< reuse workspace [F]

    ! local variables  .........................................................

    real(RNP), allocatable, save :: u_i  (:,:,:,:,:) ! intermediate solution
    real(RNP), allocatable, save :: nu   (:,:,:,:,:) ! variable diffusivity
    real(RNP), allocatable, save :: F_d1 (:,:,:,:,:) ! ∇·ν∇u
    real(RNP), allocatable, save :: F_d3 (:,:,:,:,:) ! -∇ν(∇·v)
    real(RNP), allocatable, save :: w    (:,:,:,:,:) ! workspace for u
    real(RNP), allocatable, save :: dp   (:,:,:,:)   ! pressure correction

    real(RNP), parameter :: chi = -ONE

    real(RNP), save :: t_0
    integer :: i

    associate( problem => this % problem          &
             , flow_op => this % flow_op          &
             , mesh    => this % flow_op % mesh   &
             , x       => this % flow_op % x      &
             , eop     => this % flow_op % eop_u  &
             , p       => u(:,:,:,:,4)            )

      ! initialization .........................................................

#ifdef CHECK
  !$omp master
  print '(A,G0)', "*** TimeStep *** #1  rank = ", mesh%part
  !$omp end master
#endif

      !$omp single

      ! check for reusable workspace
      if (allocated(u_i)) then
        if (any(shape(u_i) /= shape(u))) then
          deallocate(u_i, F_d1, F_d3, w, dp)
          if (allocated(nu)) deallocate(nu)
        end if
      end if

      ! allocate new workspace, if necessary
      if (.not. allocated(u_i)) then
        allocate(u_i, F_d1, F_d3, w, mold = u)
        allocate(dp, mold = p)
        if (problem % HasVariableProperties()) then
          allocate(nu, mold = u)
        end if
      end if

      t_0 = t
      t   = t + dt

      !$omp end single

      ! nu = ν₀ = ν(x,t₀,u₀)
      if (problem % HasVariableProperties()) then
        call problem % GetDiffusivity(flow_op%x, t_0, u, nu)
      end if

      ! explicit + extrapolated diffusive parts ................................

#ifdef CHECK
  !$omp master
  print '(A,G0)', "*** TimeStep *** #2  rank = ", mesh%part
  !$omp end master
#endif

      ! u' = u₀ ≡ u(x,t₀)
      call SetArray(u_i, u, multi=.true.)

      ! update boundary conditions
      call GetBoundaryValues(problem, mesh, flow_op%bv_x, t_0, flow_op%bv_u)

      ! w = -∇⋅v₀v₀ + ∇⋅ν₀(∇v₀)ᵀ - 2∇(ν₀∇⋅v₀)    for v
      ! w = -∇⋅v₀u₀                              for u \ (v,p)
      call TimeDerivative( problem, flow_op, t_0  &
                         , u      = u             &
                         , chi    = chi           &
                         , nu     = nu            &
                         , F      = w             &
                         , F_d1   = F_d1          &
                         , F_d3   = F_d3          &
                         , source = .false.       &
                         )

      ! u' += ∆t w
      do i = 1, problem % nc
        if (i == 4) cycle
        call MergeArrays(ONE, u_i(:,:,:,:,i), dt, w   (:,:,:,:,i))
        call MergeArrays(ONE, u_i(:,:,:,:,i), dt, F_d3(:,:,:,:,i))
      end do

      ! implicit sources .......................................................

#ifdef CHECK
  !$omp master
  print '(A,G0)', "*** TimeStep *** #3  rank = ", mesh%part
  !$omp end master
#endif

      call problem % GetExternalSources(flow_op%x, t, w)

      do i = 1, problem % nc
        if (i == 4) cycle
        call MergeArrays(ONE, u_i(:,:,:,:,i), dt, w(:,:,:,:,i))
      end do

      ! pressure, continuity and diffusion .....................................

#ifdef CHECK
  !$omp master
  print '(A,G0)', "*** TimeStep *** #4  rank = ", mesh%part
  !$omp end master
#endif

      ! solve for p = p"
      call PressureSolver(problem, flow_op, dt, u_i, p, w)

#ifdef CHECK
  !$omp master
  print '(A,G0)', "*** TimeStep *** #5  rank = ", mesh%part
  !$omp end master
#endif

      ! v" = v' - 1/∆t ∇p" - J(v")
      call ProjectionStep(problem, flow_op, dt, p, u_i, w)

      ! starting values and RHS for diffusion
      do i = 1, problem % nc
        if (i == 4) cycle
        call MergeArrays(ONE, u_i(:,:,:,:,i), -dt, F_d1(:,:,:,:,i))
        if (i > 3) cycle
        call MergeArrays(ONE, u_i(:,:,:,:,i), -dt, F_d3(:,:,:,:,i))
      end do

#ifdef CHECK
  !$omp master
  print '(A,G0)', "*** TimeStep *** #6  rank = ", mesh%part
  !$omp end master
#endif

      ! solve implicit diffusive part for u'''
      call DiffusionStep(problem, flow_op, dt, f=u_i, u=u, w=w, nu=nu)

      ! final projection .......................................................

      if (this % project > 0) then

        ! solve for p = p" + dp
        call SetArray(dp, ZERO)
        call PressureSolver(problem, flow_op, dt, u, dp, w) ! u=u(x,t₀+∆t)
        call MergeArrays(ONE, p, ONE, dp)
        ! v = v''' - 1/∆t ∇p - J(v)
        call ProjectionStep(problem, flow_op, dt, dp, u, w)

      end if

      ! pressure ...............................................................

      if (this % pressure > 0) then
        if (problem % HasVariableProperties()) then
          call problem % GetDiffusivity(flow_op%x, t, u, nu)
        end if
        ! recompute pressure using F_d1 as workspace for F_v
        call TimeDerivative(problem, flow_op, t, u, nu=nu, chi=chi, F=F_d1)
        call ComputePressure(problem, flow_op, t, F_d1, p, w)
      end if

    end associate

    ! clean-up .................................................................

#ifdef CHECK
  !$omp master
  print '(A,G0)', "*** TimeStep *** #7  rank = ", mesh%part
  !$omp end master
#endif

    ! keep workspace in case of standby
    if (present(standby)) then
      if (standby) return
    end if

    !$omp single
    if (allocated(u_i )) deallocate(u_i )
    if (allocated(nu  )) deallocate(nu  )
    if (allocated(F_d1)) deallocate(F_d1)
    if (allocated(F_d3)) deallocate(F_d3)
    if (allocated(w   )) deallocate(w   )
    if (allocated(dp  )) deallocate(dp  )
    !$omp end single

  end subroutine TimeStep

  !=============================================================================

end module CART__ISP_Flow__Time_Integrator__Euler
