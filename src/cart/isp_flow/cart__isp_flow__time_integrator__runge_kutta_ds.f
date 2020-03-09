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
    type(IMEX_RK_Method) :: imex_rk
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

  subroutine TimeStep(this, t, dt, u, nu)
    class(TimeIntegrator_RungeKuttaDS), intent(inout) :: this
    real(RNP),           intent(inout) :: t             !< time t₀ → t
    real(RNP),           intent(in)    :: dt            !< step size ∆t = t-t₀
    real(RNP), optional, intent(in)    :: nu(:,:,:,:,:) !< variable ν(x,t₀)
    real(RNP),           intent(inout) :: u (:,:,:,:,:) !< u(x,t₀) → u(x,t)

    ! local variables  .........................................................


  end subroutine TimeStep

  !-----------------------------------------------------------------------------
  !> Executes one IMEX Runge-Kutta stage

  subroutine RungeKuttaStage()

  end subroutine RungeKuttaStage

  !=============================================================================

end module CART__ISP_Flow__Time_Integrator__Runge_Kutta_DS
