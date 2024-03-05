!> summary:  Base class of time integrators for 1D conservation laws
!> author:   Robin Fränzel, Joerg Stiller
!> date:     2023/05/04
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module CL__Time_Integrator__1D

  use, intrinsic :: ISO_Fortran_Env, only: OUTPUT_UNIT

  use Kind_Parameters, only: RNP
  use CL__Operator__1D
  use CL__Problem__1D

  implicit none
  private

  public :: CL_TimeIntegrator_1D
  public :: CL_TimeIntegrator_Options_1D

  !-----------------------------------------------------------------------------
  !> Abstract type of a one-step time integrator for 1D conservation laws

  type, abstract :: CL_TimeIntegrator_1D

    character(len=80) :: name = '' !< time-integrator name
    integer   :: imex_mode         !< IMEX approach, 0/1/2: EX/IMEX/EX→IM
    integer   :: diffusion_method  !< implicit diffusion method
    integer   :: diffusion_i_max   !< max num iterations
    real(RNP) :: diffusion_r_red   !< residual reduction
    real(RNP) :: diffusion_r_max   !< max residual

  contains

    procedure, non_overridable :: Init_CL_TimeIntegrator_1D
    procedure, non_overridable :: Show_CL_TimeIntegrator_1D
    procedure :: Show => Show_CL_TimeIntegrator_1D
    procedure(TimeStep), deferred :: TimeStep
  end type CL_TimeIntegrator_1D

  abstract interface

    !---------------------------------------------------------------------------
    !> Execution of a single time step

    subroutine TimeStep(this, cl_problem, cl_operator, dt, t_0, u_0, u)
      import
      class(CL_TimeIntegrator_1D), intent(in) :: this
      class(CL_Problem_1D),  intent(in)    :: cl_problem
      class(CL_Operator_1D), intent(in)    :: cl_operator
      real(RNP),             intent(in)    :: dt          !< step size ∆t
      real(RNP),             intent(in)    :: t_0         !< initial time
      real(RNP), contiguous, intent(in)    :: u_0(0:,:,:) !< u(t₀)
      real(RNP), contiguous, intent(inout) :: u  (0:,:,:) !< u(t₀+∆t)
    end subroutine TimeStep

  end interface

  !-----------------------------------------------------------------------------
  !> Base type for providing time integrator options
  !>
  !> Use `imex_mode` to select the implicity of the time-integration method:
  !>
  !>   - `0`  fully explicit (EX)
  !>   - `1`  conventional implicit-explicit (IMEX)
  !>   - `2`  fully explicit step followed by implicit diffusion
  !>
  !> `diffusion_method` specifies the method for solving the algebraic systems
  !> The method for solving the algebraic system that result from the implicit
  !> treatment of diffusion. Common choices are:
  !>
  !>   - `1`  direct hybrid solver (hybrid DG, constant coefficients)
  !>   - `2`  Conjugate gradient method
  !>   - `3`  Schwarz method (mainly as smoother)
  !>   - `4`  IPCG or FGMRES (as a solver or smoother)
  !>
  !> Especially the latter two require proper initialization of `cl_operator`.
  !> See there and the problem specific implementations for more details.

  type CL_TimeIntegrator_Options_1D
    integer   :: imex_mode        = 0     !< IMEX approach
    integer   :: diffusion_method = 1     !< implicit diffusion method
    integer   :: diffusion_i_max  = 10    !< max num iterations
    real(RNP) :: diffusion_r_red  = 1e-10 !< residual reduction
    real(RNP) :: diffusion_r_max  = 1e-12 !< max residual
  end type CL_TimeIntegrator_Options_1D

contains

  !=============================================================================
  ! CL_TimeIntegrator_1D: type-bound procedures

  !-----------------------------------------------------------------------------
  !> Initialization of CL_TimeIntegrator_1D object

  subroutine Init_CL_TimeIntegrator_1D(this, opt)
    class(CL_TimeIntegrator_1D),         intent(inout) :: this
    class(CL_TimeIntegrator_Options_1D), intent(in)    :: opt

    this % imex_mode        = opt % imex_mode
    this % diffusion_method = opt % diffusion_method
    this % diffusion_i_max  = opt % diffusion_i_max
    this % diffusion_r_red  = opt % diffusion_r_red
    this % diffusion_r_max  = opt % diffusion_r_max

  end subroutine Init_CL_TimeIntegrator_1D

  !-----------------------------------------------------------------------------
  !> Output of CL_TimeIntegrator_1D settings

  subroutine Show_CL_TimeIntegrator_1D(this, unit)
    class(CL_TimeIntegrator_1D), intent(in) :: this
    integer, optional, intent(in) :: unit  !< output unit

    integer :: io

    if (present(unit)) then
      io = unit
    else
      io = OUTPUT_UNIT
    end if

    write(io,'(/,A)') 'CL_TimeIntegrator_1D settings'
    write(io,'(A,/)') repeat('=',80)

    write(io,'(2X,2A)') 'name: ', trim(this % name)

    write(io,'(/,A)') 'Solver'
    write(io,'(A,/)') repeat('-',80)
    write(io,'(2X,A,T22,I0)')     'imex_mode:'        , this % imex_mode
    write(io,'(2X,A,T22,I0)')     'diffusion_method:' , this % diffusion_method
    write(io,'(2X,A,T22,I0)')     'diffusion_i_max:'  , this % diffusion_i_max
    write(io,'(2X,A,T21,ES12.5)') 'diffusion_r_red:'  , this % diffusion_r_red
    write(io,'(2X,A,T21,ES12.5)') 'diffusion_r_max:'  , this % diffusion_r_max

  end subroutine Show_CL_TimeIntegrator_1D

  !=============================================================================

end module CL__Time_Integrator__1D
