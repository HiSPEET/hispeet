!> summary:  Euler method for incompressible flows with dual splitting
!> author:   Joerg Stiller
!> date:     2020/03/05
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> @todo
!>   * handling of variable viscosity
!===============================================================================

module CART__ISP_Flow__Time_Integrator__Euler
  use Kind_Parameters, only: RNP
  use Constants,       only: ONE, ZERO
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
    procedure :: TimeStep
  end type TimeIntegrator_Euler

  ! overloading the constructor
  interface TimeIntegrator_Euler
    module procedure New_TimeIntegrator_Euler
  end interface

  !-----------------------------------------------------------------------------
  !> Type for providing Euler- time-integrator options (none, so far)

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

    call this % Init_TimeIntegrator(problem, flow_op)

    ! no options, so far
    if (present(opt)) return

  end subroutine Init_TimeIntegrator_Euler

  !-----------------------------------------------------------------------------
  !> Performs a single IMEX Euler step

  subroutine TimeStep(this, t, dt, u)
    class(TimeIntegrator_Euler), intent(inout) :: this
    real(RNP), intent(inout) :: t             !< time t₀ → t
    real(RNP), intent(in)    :: dt            !< step size ∆t = t-t₀
    real(RNP), intent(inout) :: u (:,:,:,:,:) !< u(x,t₀) → u(x,t)

    ! local variables  .........................................................

    real(RNP), allocatable, save :: u_i  (:,:,:,:,:) ! intermediate solution
    real(RNP), allocatable, save :: nu   (:,:,:,:,:) ! variable diffusivity
    real(RNP), allocatable, save :: F_d1 (:,:,:,:,:) ! ∇·ν∇u
    real(RNP), allocatable, save :: F_d3 (:,:,:,:,:) ! -χ∇ν(∇·v)
    real(RNP), allocatable, save :: w    (:,:,:,:,:) ! workspace for u
    real(RNP), allocatable, save :: dp   (:,:,:,:)   ! pressure correction

    associate( problem => this % problem          &
             , flow_op => this % flow_op          &
             , mesh    => this % flow_op % mesh   &
             , x       => this % flow_op % x      &
             , eop     => this % flow_op % eop_u  &
             , p       => u(:,:,:,:,4)            )

      ! initialization .........................................................

      ! workspace
      !$omp single
      allocate(u_i , mold = u)
      allocate(F_d1, mold = u)
      allocate(F_d3, mold = u)
      allocate(w   , mold = u)
      allocate(dp  , mold = p)
      if (problem % HasVariableProperties()) then
        allocate(nu, mold = u)
      end if
      !$omp end single

      ! nu = ν₀ = ν(x,t₀,u₀)
      if (problem % HasVariableProperties()) then
        call problem % GetDiffusivity(flow_op%x, t, u, nu)
      end if

      t = t + dt

      ! boundary conditions ....................................................

      call GetBoundaryValues(problem, mesh, flow_op%bv_x, t, flow_op%bv_u)

      ! explicit + extrapolated diffusive parts ................................

      ! u' = u₀ ≡ u(x,t₀)
      call SetArray(u_i, u, multi=.true.)

      ! w = -∇⋅v₀v₀ + ∇⋅ν₀(∇v₀)ᵀ - χ∇(ν₀∇⋅v₀)     for v
      ! w = -∇⋅v₀u₀                               for u \ (v,p)
      call TimeDerivative( problem, flow_op, t  &
                         , u_c  = u             &
                         , u_d  = u             &
                         , nu   = nu            &
                         , F    = w             &
                         , F_d1 = F_d1          &
                         , F_d3 = F_d3          &
                         )

      ! u' += ∆t w
      call MergeArrays(ONE, u_i, dt, w, multi=.true.)

      ! pressure, continuity and diffusion .....................................

      ! solve for p = p"
      call PressureSolver(problem, flow_op, dt, u_i, p, w)

      ! v" = v' - 1/∆t ∇p" - J(v")
      call ProjectionStep(problem, flow_op, dt, p, u_i, w)

      ! solve implicit diffusive part for u'''
      call MergeArrays(ONE, u_i, -dt, F_d1, multi=.true.)
      call MergeArrays(ONE, u_i, -dt, F_d3, multi=.true.)
      call DiffusionStep(problem, flow_op, dt, f=u_i, u=u, w=w, nu=nu)

      ! final projection .......................................................

      if (flow_op % control % div_final) then

        ! solve for p = p" + dp
        call SetArray(dp, ZERO)
        call PressureSolver(problem, flow_op, dt, u, dp, w) ! u=u(x,t₀+∆t)
        call MergeArrays(ONE, p, ONE, dp)
        ! v = v''' - 1/∆t ∇p - J(v)
        call ProjectionStep(problem, flow_op, dt, dp, u, w)

      end if

      ! clean-up ...............................................................

      !$omp barrier
      !$omp master
      if (allocated(u_i )) deallocate(u_i )
      if (allocated(nu  )) deallocate(nu  )
      if (allocated(F_d1)) deallocate(F_d1)
      if (allocated(F_d3)) deallocate(F_d3)
      if (allocated(w   )) deallocate(w   )
      if (allocated(dp  )) deallocate(dp  )
      !$omp end master

    end associate

  end subroutine TimeStep

  !=============================================================================

end module CART__ISP_Flow__Time_Integrator__Euler
