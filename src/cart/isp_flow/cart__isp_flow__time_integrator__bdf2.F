!> summary:  IMEX BDF2 method for incompressible flows with dual splitting
!> author:   Joerg Stiller
!> date:     2020/07/26
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> @note
!>   * uses internal workspace for storing results dating from t₀-∆t
!===============================================================================

module CART__ISP_Flow__Time_Integrator__BDF2
  use Kind_Parameters, only: RNP
  use Constants,       only: ZERO, ONE, HALF
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

  public :: TimeIntegrator_BDF2
  public :: TimeIntegrator_BDF2_Options

  !-----------------------------------------------------------------------------
  !> IMEX BDF2 method for incompressible flows with dual splitting

  type, extends(TimeIntegrator) :: TimeIntegrator_BDF2
  contains
    procedure :: Init_TimeIntegrator_BDF2
    procedure :: Show => Show_TimeIntegrator_BDF2
    procedure :: TimeStep
  end type TimeIntegrator_BDF2

  ! overloading the constructor
  interface TimeIntegrator_BDF2
    module procedure New_TimeIntegrator_BDF2
  end interface

  !-----------------------------------------------------------------------------
  !> Type for providing BDF2 time-integrator options (none, so far)
  !>
  !> Options inherited from base class
  !>
  !>   * `splitting` -- defines the splitting scheme used in the stages:      \n
  !>        1: standard velocity correction with χ = -1                       \n
  !>        2: rotational velocity correction with χ = -2                     \n
  !>           and F_d3/2 removed after extrapolation                         \n
  !>        3: "native" velocity correction with χ = 0                        \n
  !>        4: velocity correction with χ = -1                                \n
  !>           and F_d3 removed after extrapolation                           \n
  !>        5: velocity correction with χ = -2                                \n
  !>           and F_d3 removed after extrapolation
  !>
  !>   * `project`                                                            \n
  !>        0: no additional projection step                                  \n
  !>        1: additional projection at the end of the time step
  !>
  !>   * `pressure`
  !>        0: return pressure as computed                                    \n
  !>        1: recompute pressure at the end of the time step

  type, extends(TimeIntegratorOptions) :: TimeIntegrator_BDF2_Options
  end type TimeIntegrator_BDF2_Options

contains

  !-----------------------------------------------------------------------------
  !> Constructor for objects of type TimeIntegrator_BDF2 with options

  function New_TimeIntegrator_BDF2(problem, flow_op, opt) result(this)
    class(FlowProblem),                           intent(in) :: problem
    class(FlowOperators),                         intent(in) :: flow_op
    class(TimeIntegrator_BDF2_Options), optional, intent(in) :: opt
    type(TimeIntegrator_BDF2) :: this

    call Init_TimeIntegrator_BDF2(this, problem, flow_op, opt)

  end function New_TimeIntegrator_BDF2

  !-----------------------------------------------------------------------------
  !> Initialization of a Init_TimeIntegrator_BDF2 object

  subroutine Init_TimeIntegrator_BDF2(this, problem, flow_op, opt)
    class(TimeIntegrator_BDF2),                   intent(inout) :: this
    class(FlowProblem),                            intent(in)    :: problem
    class(FlowOperators),                          intent(in)    :: flow_op
    class(TimeIntegrator_BDF2_Options), optional, intent(in)    :: opt

    ! intialize parent type
    call this % Init_TimeIntegrator(problem, flow_op, opt)
    this % name = 'IMEX BDF2 method'

  end subroutine Init_TimeIntegrator_BDF2

  !-----------------------------------------------------------------------------
  !> Output of TimeIntegrator_BDF2 settings

  subroutine Show_TimeIntegrator_BDF2(this, unit)
    class(TimeIntegrator_BDF2), intent(in) :: this
    integer,           optional, intent(in) :: unit  !< output unit

    ! show parent settings
    call this % Show_TimeIntegrator(unit)

  end subroutine Show_TimeIntegrator_BDF2

  !-----------------------------------------------------------------------------
  !> Performs a single IMEX BDF2 step

  subroutine TimeStep(this, t, dt, u, standby)
    class(TimeIntegrator_BDF2), intent(inout) :: this
    real(RNP),         intent(inout) :: t             !< time t₀ → t
    real(RNP),         intent(in)    :: dt            !< step size ∆t = t-t₀
    real(RNP),         intent(inout) :: u (:,:,:,:,:) !< u(x,t₀) → u(x,t)
    logical, optional, intent(in)    :: standby       !< reuse workspace [T]

    ! local variables  .........................................................

    ! ∆F is the part of extrapolated RHS to be removed after projection
    real(RNP), allocatable, save :: u_1  (:,:,:,:,:) ! solution at t₀-∆t
    real(RNP), allocatable, save :: u_i  (:,:,:,:,:) ! intermediate solution
    real(RNP), allocatable, save :: nu   (:,:,:,:,:) ! var. diffusivity (t₀   )
    real(RNP), allocatable, save :: nu_1 (:,:,:,:,:) ! var. diffusivity (t₀-∆t)
    real(RNP), allocatable, save :: nu_x (:,:,:,:,:) ! var. diffusivity (t₀+∆t)
    real(RNP), allocatable, save :: F_c  (:,:,:,:,:) ! ∇·vu
    real(RNP), allocatable, save :: F_d1 (:,:,:,:,:) ! ∇·ν∇u
    real(RNP), allocatable, save :: F_d2 (:,:,:,:,:) ! ∇·ν(∇u)ᵀ
    real(RNP), allocatable, save :: F_d3 (:,:,:,:,:) ! χ∇ν(∇·v)
    real(RNP), allocatable, save :: dF   (:,:,:,:,:) ! ∆F
    real(RNP), allocatable, save :: w    (:,:,:,:,:) ! workspace for u

    ! IMEX BDF2 coefficients
    real(RNP), parameter :: gamma_0 =  3 * HALF
    real(RNP), parameter :: alpha_0 =  2
    real(RNP), parameter :: alpha_1 = -HALF
    real(RNP), parameter :: beta_0  =  2
    real(RNP), parameter :: beta_1  = -1

    ! control
    real(RNP), save :: t_0, t_1 = -huge(ONE)
    logical,   save :: euler

    ! auxiliary
    real(RNP) :: a0, a1, tau, tau_cd3, tau_s
    integer   :: np, ne, nc
    integer   :: c, e, i, j, k

    np = size(u,1)
    ne = size(u,4)
    nc = size(u,5)

    associate( problem => this % problem          &
             , flow_op => this % flow_op          &
             , mesh    => this % flow_op % mesh   &
             , x       => this % flow_op % x      &
             , eop     => this % flow_op % eop_u  &
             , chi     => this % chi              &
             , cd3     => this % cd3              &
             , p       => u(:,:,:,:,4)            )

      ! initialization .........................................................


      ! workspace  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

      !$omp single

      if (allocated(u_1)) then
        if (any(shape(u_1) /= shape(u))) then
          deallocate(u_1, u_i, F_c, F_d1, F_d2, F_d3, dF, w)
          if (allocated(nu  )) deallocate(nu)
          if (allocated(nu_1)) deallocate(nu_1)
          if (allocated(nu_x)) deallocate(nu_x)
        end if
      end if

      if (allocated(u_1)) then

        euler = abs(t_1 + dt - t) > epsilon(ONE)

      else

        euler = .true.

        allocate( u_1  , mold = u )
        allocate( u_i  , mold = u )
        allocate( F_c  , mold = u )
        allocate( F_d1 , mold = u )
        allocate( F_d2 , mold = u )
        allocate( F_d3 , mold = u )
        allocate( dF   , mold = u )
        allocate( w    , mold = u )

        if (problem % HasVariableProperties()) then
          allocate( nu   , mold = u )
          allocate( nu_1 , mold = u )
          allocate( nu_x , mold = u )
        end if

      end if

      t_0 = t
      t   = t + dt

      !$omp end single

      ! variable viscosity
      if (problem % HasVariableProperties()) then
        ! nu = ν(*,t₀)
        call problem % GetDiffusivity(flow_op%x, t_0, u, nu)
        ! nu_x ≈ ν(*,t₀+∆t)
        if (euler) then
          call SetArray(nu_x, nu, multi = .true.)
        else
          call SetArray(nu_x, nu, multi = .true.)
          call MergeArrays(beta_0, nu_x, beta_1, nu_1, multi = .true.)
        end if
        call SetArray(nu_1, nu, multi = .true.)
      end if

      ! boundary conditions ....................................................

      call GetBoundaryValues(problem, mesh, flow_op%bv_x, t, flow_op%bv_u)

      ! explicit + extrapolated diffusive parts ................................

      if (euler) then

        call SetArray(u_i, u   , multi = .true.)
        call SetArray(dF , ZERO, multi = .true.)

      else ! BDF2 with u_1, F_c, F_d1, F_d2, F_d3 storing values from t₀-∆t

        a0       =  alpha_0 / gamma_0
        a1       =  alpha_1 / gamma_0
        tau      =  dt / gamma_0 * beta_1
        tau_cd3  =  tau*cd3

        !$omp do collapse(2)
        do c = 1, 3
        do e = 1, ne
          do k = 1, np
          do j = 1, np
          do i = 1, np

            u_i(i,j,k,e,c)  =  a0  *   u    (i,j,k,e,c) &
                            +  a1  *   u_1  (i,j,k,e,c) &
                            +  tau * ( F_c  (i,j,k,e,c) &
                                     + F_d1 (i,j,k,e,c) &
                                     + F_d2 (i,j,k,e,c) &
                                     + F_d3 (i,j,k,e,c) )

            dF(i,j,k,e,c)  =  tau     * F_d1(i,j,k,e,c)  &
                           +  tau_cd3 * F_d3(i,j,k,e,c)

          end do
          end do
          end do
        end do
        end do

        !$omp do collapse(2)
        do c = 5, nc
        do e = 1, ne
          do k = 1, np
          do j = 1, np
          do i = 1, np

            u_i(i,j,k,e,c)  =  a0  *   u    (i,j,k,e,c) &
                            +  a1  *   u_1  (i,j,k,e,c) &
                            +  tau * ( F_c  (i,j,k,e,c) &
                                     + F_d1 (i,j,k,e,c) )

            dF(i,j,k,e,c)  =  tau * F_d1(i,j,k,e,c)

          end do
          end do
          end do
        end do
        end do

      end if

      call TimeDerivative( problem, flow_op, t  &
                         , u    = u             &
                         , nu   = nu            &
                         , chi  = chi           &
                         , F_c  = F_c           &
                         , F_d1 = F_d1          &
                         , F_d2 = F_d2          &
                         , F_d3 = F_d3          &
                         , F_s  = w             &
                         )

      ! save t₀ and u(x,t₀) for next step
      !$omp single
      t_1 = t_0
      !$omp end single
      call SetArray(u_1, u, multi = .true.)

      if (euler) then
        tau   = dt
        tau_s = dt
      else
        tau   = dt / gamma_0 * beta_0
        tau_s = dt / gamma_0
      end if
      tau_cd3 = tau*cd3

      !$omp do collapse(2)
      do c = 1, 3
      do e = 1, ne
        do k = 1, np
        do j = 1, np
        do i = 1, np

          u_i(i,j,k,e,c)  =          u_i  (i,j,k,e,c)   &
                          +  tau * ( F_c  (i,j,k,e,c)   &
                                   + F_d1 (i,j,k,e,c)   &
                                   + F_d2 (i,j,k,e,c)   &
                                   + F_d3 (i,j,k,e,c) ) &
                          +  tau_s * w    (i,j,k,e,c)

          dF(i,j,k,e,c)  =            dF   (i,j,k,e,c)  &
                         +  tau     * F_d1 (i,j,k,e,c)  &
                         +  tau_cd3 * F_d3 (i,j,k,e,c)

        end do
        end do
        end do
      end do
      end do

      !$omp do collapse(2)
      do c = 5, nc
      do e = 1, ne
        do k = 1, np
        do j = 1, np
        do i = 1, np

          u_i(i,j,k,e,c)  =          u_i  (i,j,k,e,c)   &
                          +  tau * ( F_c  (i,j,k,e,c)   &
                                   + F_d1 (i,j,k,e,c) ) &
                          +  tau_s * w    (i,j,k,e,c)

          dF(i,j,k,e,c)  =        dF   (i,j,k,e,c)  &
                         +  tau * F_d1 (i,j,k,e,c)

        end do
        end do
        end do
      end do
      end do

      ! pressure, continuity and diffusion .....................................

      if (euler) then
        tau = dt
      else
        tau = dt / gamma_0
      end if

      ! solve for p = p"
      call PressureSolver(problem, flow_op, tau, u_i, p, w)

      ! v" = v' - 1/τ ∇p" - J(v")
      call ProjectionStep(problem, flow_op, tau, p, u_i, w)

      ! starting values and RHS for diffusion
      do i = 1, problem % nc
        if (i == 4) cycle
        call SetArray(u(:,:,:,:,i), u_i(:,:,:,:,i))
        call MergeArrays(ONE, u_i(:,:,:,:,i), -ONE, dF(:,:,:,:,i))
      end do

      ! solve implicit diffusive part for u'''
      call DiffusionStep(problem, flow_op, tau, f=u_i, u=u, w=w, nu=nu_x)

      ! final projection .......................................................

      if (this % project > 0) then
        associate(dp => dF(:,:,:,:,4))

          ! solve for p = p" + dp
          call SetArray(dp, ZERO)
          call PressureSolver(problem, flow_op, tau, u, dp, w)
          call MergeArrays(ONE, p, ONE, dp)
          ! v = v''' - 1/τ ∇p - J(v)
          call ProjectionStep(problem, flow_op, tau, dp, u, w)

        end associate
      end if

      ! pressure ...............................................................

      if (this % pressure > 0) then
        if (problem % HasVariableProperties()) then
          call problem % GetDiffusivity(flow_op%x, t, u, nu)
        end if
        ! recompute pressure using dF as workspace for F_v
        call TimeDerivative(problem, flow_op, t, u, nu=nu, chi=chi, F=dF)
        call ComputePressure(problem, flow_op, t, dF, p, w)
      end if

    end associate

    ! clean-up .................................................................

    if (.not. present(standby)) then
      return
    else if (standby) then
      return
    end if

    !$omp single
    if (allocated( u_1  )) deallocate( u_1  )
    if (allocated( u_i  )) deallocate( u_i  )
    if (allocated( nu   )) deallocate( nu   )
    if (allocated( nu_1 )) deallocate( nu_1 )
    if (allocated( nu_x )) deallocate( nu_x )
    if (allocated( F_c  )) deallocate( F_c  )
    if (allocated( F_d1 )) deallocate( F_d1 )
    if (allocated( F_d2 )) deallocate( F_d2 )
    if (allocated( F_d3 )) deallocate( F_d3 )
    if (allocated( dF   )) deallocate( dF   )
    if (allocated( w    )) deallocate( w    )
    !$omp end single

  end subroutine TimeStep

  !=============================================================================

end module CART__ISP_Flow__Time_Integrator__BDF2
