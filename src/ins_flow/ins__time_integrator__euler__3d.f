!> summary:  Euler method for incompressible flows
!> author:   Joerg Stiller
!> date:     2022/09/09
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module INS__Time_Integrator__Euler__3D
  use Kind_Parameters
  use XMPI
  use INS__Time_Integrator__3D
  use INS__Problem__3D
  use INS__Operator__3D
  implicit none
  private

  public :: INS_TimeIntegrator_Euler_3D
  public :: INS_TimeIntegrator_Euler_Options_3D

  !-----------------------------------------------------------------------------
  !> Euler method for incompressible flows

  type, extends(INS_TimeIntegrator_3D) :: INS_TimeIntegrator_Euler_3D
  contains
    procedure, non_overridable :: Init_INS_TimeIntegrator_Euler_3D
    procedure :: TimeStep
  end type INS_TimeIntegrator_Euler_3D

  ! constructor
  interface INS_TimeIntegrator_Euler_3D
    module procedure New_INS_TimeIntegrator_Euler_3D
  end interface

  !-----------------------------------------------------------------------------
  !> Type for providing Euler time-integrator options (none, so far)

  type, extends(INS_TimeIntegratorOptions_3D) :: &
    INS_TimeIntegrator_Euler_Options_3D
  end type INS_TimeIntegrator_Euler_Options_3D

contains

  !-----------------------------------------------------------------------------
  !> Constructor for objects of type INS_TimeIntegrator_Euler_3D

  function New_INS_TimeIntegrator_Euler_3D(problem, flow_op, opt) result(this)
    class(INS_Problem_3D),  intent(in) :: problem
    class(INS_Operator_3D), intent(in) :: flow_op
    class(INS_TimeIntegrator_Euler_Options_3D), optional, intent(in) :: opt
    type(INS_TimeIntegrator_Euler_3D) :: this

    call Init_INS_TimeIntegrator_Euler_3D(this, problem, flow_op, opt)

  end function New_TimeIntegrator_Euler

  !-----------------------------------------------------------------------------
  !> Initialization of a INS_TimeIntegrator_Euler_3D object

  subroutine Init_INS_TimeIntegrator_Euler_3D(this, problem, flow_op, opt)
    class(INS_TimeIntegrator_Euler_3D), intent(inout) :: this
    class(INS_Problem_3D),              intent(in)    :: problem
    class(INS_Operator_3D),             intent(in)    :: flow_op
    class(INS_TimeIntegrator_Euler_Options_3D), optional, intent(in) :: opt

    ! intialize parent type
    call this % Init_INS_TimeIntegrator_3D(problem, flow_op, opt)
    this % name = 'Euler method'

  end subroutine Init_INS_TimeIntegrator_Euler_3D

  !---------------------------------------------------------------------------
  !> Execution of an Euler time step

  subroutine TimeStep(this, t, dt, u, standby)
    class(INS_TimeIntegrator_Euler_3D), intent(inout) :: this
    real(RNP),         intent(inout) :: t            !< time t₀ → t
    real(RNP),         intent(in)    :: dt           !< step size ∆t = t-t₀
    real(RNP),         intent(inout) :: u(:,:,:,:,:) !< u(x,t₀) → u(x,t)
    logical, optional, intent(in)    :: standby      !< reuse workspace
  end subroutine TimeStep

  !=============================================================================

end module INS__Time_Integrator__Euler__3D
