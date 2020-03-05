!> summary:  Base type of one-step time integrators for incompressible flows
!> author:   Joerg Stiller
!> date:     2020/03/05
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module CART__ISP_Flow__Time_Integrator
  use Kind_Parameters, only: RNP
  implicit none
  private

  public :: TimeIntegrator

  !-----------------------------------------------------------------------------
  !> Abstract type for defining a one-step time integrator incompressible flow

  type, abstract :: TimeIntegrator
    class(FlowProblem),   pointer :: problem => null() !< flow problem
    class(FlowOperators), pointer :: flow_op => null() !< flow operators
  contains
    procedure(TimeStep), deferred :: TimeStep
  end type TimeIntegrator

  abstract interface

    !---------------------------------------------------------------------------
    !> Execution of a single time step with optional variable viscosity ν

    subroutine TimeStep(this, t, dt, u, nu)
      import
      class(TimeIntegrator), intent(inout) :: this
      real(RNP),             intent(inout) :: t             !< time t₀ → t
      real(RNP),             intent(in)    :: dt            !< step size ∆t = t-t₀
      real(RNP), optional,   intent(in)    :: nu(:,:,:,:,:) !< ν(x,t₀)
      real(RNP),             intent(inout) :: u (:,:,:,:,:) !< u(x,t₀) → u(x,t)
    end subroutine TimeStep

  end interface

  !=============================================================================

end module CART__ISP_Flow__Time_Integrator
