!------------------------------------------------------------------------------!
! This file is part of HiSPEET: High-order Spectral Element Techniques         !
!                                                                              !
! Copyright (C) 2026 by the HiSPEET authors and the                            !
! Chair of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany                 !
!                                                                              !
! HiSPEET is free software: you can redistribute it and/or modify              !
! it under the terms of the GNU General Public License as published by         !
! the Free Software Foundation, either version 3 of the License, or            !
! (at your option) any later version.                                          !
!                                                                              !
! HiSPEET is distributed in the hope that it will be useful,                   !
! but WITHOUT ANY WARRANTY; without even the implied warranty of               !
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                         !
! See the GNU General Public License for more details.                         !
!                                                                              !
! You should have received a copy of the GNU General Public License            !
! along with HiSPEET. If not, see <http://www.gnu.org/licenses/>.              !
!------------------------------------------------------------------------------!

!> summary:  IMEX trapezoidal rule for incompressible flows
!> author:   Joerg Stiller
!> date:     2020/08/31
!===============================================================================

module CART__ISP_Flow__Time_Integrator__TR
  use, intrinsic :: ISO_Fortran_Env, only: OUTPUT_UNIT
  use Kind_Parameters, only: RNP
  use Constants,       only: ZERO, ONE
  use XMPI
  use Array_Assignments

  use ISP_Flow_Problem

  use CART__ISP_Flow__Boundary_Values
  use CART__ISP_Flow__Diffusion
  use CART__ISP_Flow__Operators
  use CART__ISP_Flow__Pressure
  use CART__ISP_Flow__Projection
  use CART__ISP_Flow__SGS_Diffusivity
  use CART__ISP_Flow__Time_Derivative
  use CART__ISP_Flow__Time_Integrator

  implicit none
  private

  public :: TimeIntegrator_TR
  public :: TimeIntegrator_TR_Options

  !-----------------------------------------------------------------------------
  !> IMEX trapezoidal rule for incompressible flows

  type, extends(TimeIntegrator) :: TimeIntegrator_TR
    integer, allocatable :: i_max_p1 !< max pressure  iterations in substep 1
    integer, allocatable :: i_max_u1 !< max diffusion iterations in substep 1
    integer, allocatable :: i_max_p2 !< max pressure  iterations in substep 2
    integer, allocatable :: i_max_u2 !< max diffusion iterations in substep 2
  contains
    procedure :: Init_TimeIntegrator_TR
    procedure :: Show => Show_TimeIntegrator_TR
    procedure :: TimeStep
  end type TimeIntegrator_TR

  ! overloading the constructor
  interface TimeIntegrator_TR
    module procedure New_TimeIntegrator_TR
  end interface

  !-----------------------------------------------------------------------------
  !> Type for providing TR time-integrator options (none, so far)
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

  type, extends(TimeIntegratorOptions) :: TimeIntegrator_TR_Options
    integer :: i_max_p1 = -1 !< max pressure  iterations in substep 1, if > 0
    integer :: i_max_u1 = -1 !< max diffusion iterations in substep 1, if > 0
    integer :: i_max_p2 = -1 !< max pressure  iterations in substep 2, if > 0
    integer :: i_max_u2 = -1 !< max diffusion iterations in substep 2, if > 0
  contains
    procedure :: Bcast => Bcast_TimeIntegrator_TR_Options
  end type TimeIntegrator_TR_Options

contains

  !-----------------------------------------------------------------------------
  !> Constructor for objects of type TimeIntegrator_TR with options

  function New_TimeIntegrator_TR(problem, flow_op, opt) result(this)
    class(FlowProblem),                         intent(in) :: problem
    class(FlowOperators),                       intent(in) :: flow_op
    class(TimeIntegrator_TR_Options), optional, intent(in) :: opt
    type(TimeIntegrator_TR) :: this

    call Init_TimeIntegrator_TR(this, problem, flow_op, opt)

  end function New_TimeIntegrator_TR

  !-----------------------------------------------------------------------------
  !> Initialization of a Init_TimeIntegrator_TR object

  subroutine Init_TimeIntegrator_TR(this, problem, flow_op, opt)
    class(TimeIntegrator_TR),                   intent(inout) :: this
    class(FlowProblem),                         intent(in)    :: problem
    class(FlowOperators),                       intent(in)    :: flow_op
    class(TimeIntegrator_TR_Options), optional, intent(in)    :: opt

    ! intialize parent type
    call this % Init_TimeIntegrator(problem, flow_op, opt)
    this % name = 'IMEX trapezoidal rule'

    ! specific settings
    if (opt % i_max_p1 > 0) this % i_max_p1 = opt % i_max_p1
    if (opt % i_max_u1 > 0) this % i_max_u1 = opt % i_max_u1
    if (opt % i_max_p2 > 0) this % i_max_p2 = opt % i_max_p2
    if (opt % i_max_u2 > 0) this % i_max_u2 = opt % i_max_u2

  end subroutine Init_TimeIntegrator_TR

  !-----------------------------------------------------------------------------
  !> Output of TimeIntegrator_TR settings

  subroutine Show_TimeIntegrator_TR(this, unit)
    class(TimeIntegrator_TR), intent(in) :: this
    integer,        optional, intent(in) :: unit  !< output unit

    integer :: io

    if (present(unit)) then
      io = unit
    else
      io = OUTPUT_UNIT
    end if

    ! show parent settings
    call this % Show_TimeIntegrator(unit)

    write(io,'(A)')           'TimeIntegrator_TR settings'
    write(io,'(A,/)')         repeat('-',80)

    call ShowAllocatableOption('i_max_p1', this % i_max_p1)
    call ShowAllocatableOption('i_max_u1', this % i_max_u1)
    call ShowAllocatableOption('i_max_p2', this % i_max_p2)
    call ShowAllocatableOption('i_max_u2', this % i_max_u2)

    write(io,*)

  contains

    subroutine ShowAllocatableOption(name, opt)
      character(len=*),  intent(in) :: name
      integer, optional, intent(in) :: opt

      if (present(opt)) then
        write(io,'(2X,2A,T15,G0)') name, ':', opt
      else
        write(io,'(2X,2A,T15,G0)') name, ':', 'not set'
      end if
    end subroutine ShowAllocatableOption

  end subroutine Show_TimeIntegrator_TR

  !-----------------------------------------------------------------------------
  !> Performs a single IMEX TR step

  subroutine TimeStep(this, t, dt, u, standby)
    class(TimeIntegrator_TR), intent(inout) :: this
    real(RNP),         intent(inout) :: t             !< time t₀ → t
    real(RNP),         intent(in)    :: dt            !< step size ∆t = t-t₀
    real(RNP),         intent(inout) :: u (:,:,:,:,:) !< u(x,t₀) → u(x,t)
    logical, optional, intent(in)    :: standby       !< reuse workspace [F]

    ! local variables  .........................................................

    real(RNP), allocatable, save :: u_i (:,:,:,:,:) ! intermediate solution
    real(RNP), allocatable, save :: u_1 (:,:,:,:,:) ! solution after substep 1
    real(RNP), allocatable, save :: F_0 (:,:,:,:,:) ! F - F_p at t₀
    real(RNP), allocatable, save :: F_1 (:,:,:,:,:) ! F - F_p at t₁ = t₀ + ∆t
    real(RNP), allocatable, save :: dF  (:,:,:,:,:) ! RHS correction for diffusion
    real(RNP), allocatable, save :: F_d3(:,:,:,:,:) ! -∇(ν∇⋅v)
    real(RNP), allocatable, save :: nu  (:,:,:,:,:) ! variable diffusivity
    real(RNP), allocatable, save :: w   (:,:,:,:,:) ! workspace for u

    real(RNP), parameter :: chi = -ONE

    real(RNP), save :: t_0
    integer :: i

    associate( problem => this % problem          &
             , flow_op => this % flow_op          &
             , mesh    => this % flow_op % mesh   &
             , x       => this % flow_op % x      &
             , eop     => this % flow_op % eop_u  &
             , p       => u(:,:,:,:,4)            )

      !-------------------------------------------------------------------------
      ! intialization

      !$omp single

      ! check for reusable workspace
      if (allocated(u_i)) then
        if (any(shape(u_i) /= shape(u))) then
          deallocate(u_i, u_1, F_0, F_1, dF, F_d3, w)
          if (allocated(nu  )) deallocate(nu)
        end if
      end if

      ! allocate new workspace, if necessary
      if (.not. allocated(u_i)) then
        allocate(u_i, u_1, F_0, F_1, dF, F_d3, w, mold = u)
        if ( problem % HasVariableProperties() .or.   &
             flow_op % control % sgs > 0            ) &
        then
          allocate(nu, mold = u)
        end if
      end if

      t_0 = t
      t   = t + dt

      !$omp end single

      ! boundary conditions at time t ..........................................

      call GetBoundaryValues(problem, mesh, flow_op%bv_x, t, flow_op%bv_u)

      !-------------------------------------------------------------------------
      ! substep 1

      ! FW Euler skipping pressure .............................................

      if (allocated(nu)) then
        ! molecular viscosity,  nu = ν(*,t₀)
        call problem % GetDiffusivity(flow_op%x, t_0, u, nu)
        ! add SGS viscosity
        if (flow_op % control % sgs > 0) then
          call Add_SGS_Diffusivity(flow_op, u, nu)
        end if
      end if

      ! F_0 = [F_c + F_d + F_s](x,t₀,ν₀,u₀), u_i = F_d3(ν₀,u₀)
      call TimeDerivative( problem, flow_op, t_0  &
                         , u    = u               &
                         , nu   = nu              &
                         , chi  = chi             &
                         , F    = F_0             &
                         , F_d1 = dF              &
                         , F_d3 = F_d3            &
                         )

      ! u₁' = u₀ + ∆t F_0
      call SetArray(u_i, u, multi=.true.)
      do i = 1, problem % nc
        if (i == 4) cycle
        call MergeArrays(ONE, u_i(:,:,:,:,i), dt, F_0 (:,:,:,:,i))
        call MergeArrays(ONE, u_i(:,:,:,:,i), dt, F_d3(:,:,:,:,i))
      end do

      ! pressure and continuity ................................................

      ! solve for p ≈ (p₀+p₁)
      call PressureSolver( problem, flow_op, dt, u_i, p, w &
                         , i_max = this % i_max_p1         )

      ! v₁ = v₁' - 1/∆t ∇p - J(v₁)
      call ProjectionStep(problem, flow_op, dt, p, u_i, w)

      ! BW Euler diffusion .....................................................

      do i = 1, problem % nc
        if (i == 4) cycle
        call SetArray(u_1(:,:,:,:,i), u_i(:,:,:,:,i))
        call MergeArrays(ONE, u_i(:,:,:,:,i), -dt, dF  (:,:,:,:,i))
        if (i > 3) cycle
        call MergeArrays(ONE, u_i(:,:,:,:,i), -dt, F_d3(:,:,:,:,i))
      end do

      ! solve implicit diffusive part for u₂'''
      call DiffusionStep( problem, flow_op, dt, f=u_i, u=u_1, w=w, nu=nu &
                        , i_max = this % i_max_u1                        )

      !-------------------------------------------------------------------------
      ! substep 2

      ! explicit trapezoidal rule skipping pressure  ...........................

      ! recompute ν₁ = ν(x,t₁,u₁) if required
      if ( problem % HasVariableProperties() .or.   &! ν_phy is variable
           flow_op % control % sgs > 1            ) &! ν_sgs must be updated
      then
        call problem % GetDiffusivity(flow_op%x, t, u_1, nu)
        ! add SGS viscosity, required even if no ν_sgs update requested
        if (flow_op % control % sgs > 0) then
          call Add_SGS_Diffusivity(flow_op, u_1, nu)
        end if
      end if

      ! F_1 = [F_c + F_d + F_s](x,t₁,ν₁,u₁)
      call TimeDerivative( problem, flow_op, t  &
                         , u    = u_1           &
                         , nu   = nu            &
                         , chi  = chi           &
                         , F    = F_1           &
                         , F_d1 = dF            &
                         , F_d3 = w             &
                         )

      ! save dF = F_d1(ν₁,u₁) + F_d3(ν₀,u₀) + F_d3(ν₁,u₁)
      do i = 1, 3
        call MergeArrays(ONE, dF(:,:,:,:,i), ONE, F_d3(:,:,:,:,i))
        call MergeArrays(ONE, dF(:,:,:,:,i), ONE, w   (:,:,:,:,i))
      end do

      ! implicit TR diffusion ..................................................

      ! u₂' = u₀ + ∆t/2 [(F + F_d3)₀ + (F + F_d3)₁]
      do i = 1, problem % nc
        if (i == 4) cycle
        call SetArray(u_i(:,:,:,:,i), u(:,:,:,:,i))
        call MergeArrays(ONE, u_i(:,:,:,:,i), dt/2, F_0 (:,:,:,:,i)) ! F₀
        call MergeArrays(ONE, u_i(:,:,:,:,i), dt/2, F_d3(:,:,:,:,i)) ! F_d3₀
        call MergeArrays(ONE, u_i(:,:,:,:,i), dt/2, F_1 (:,:,:,:,i)) ! F₁
        call MergeArrays(ONE, u_i(:,:,:,:,i), dt/2, w   (:,:,:,:,i)) ! F_d3₁
      end do

      ! pressure and continuity ................................................

      ! solve for p ≈ (p₀+p₂)
      call PressureSolver( problem, flow_op, dt, u_i, p, w &
                         , i_max = this % i_max_p2         )

      ! v₂" = v₂' - 1/∆t ∇p - J(v₂")
      call ProjectionStep(problem, flow_op, dt, p, u_i, w)

      ! implicit TR diffusion ..................................................

      ! starting values and RHS
      do i = 1, problem % nc
        if (i == 4) cycle
        call SetArray(u(:,:,:,:,i), u_i(:,:,:,:,i))
        call MergeArrays(ONE, u_i(:,:,:,:,i), -dt/2, dF(:,:,:,:,i))
      end do

      ! nu = ν₂ = ν(x,t₂,u₂")
      if (problem % HasVariableProperties()) then
        call problem % GetDiffusivity(flow_op%x, t, u, nu)
      end if

      ! solve implicit diffusive part for u₂'''
      call DiffusionStep( problem, flow_op, dt/2, f=u_i, u=u, w=w, nu=nu &
                        , i_max = this % i_max_u2                       )

      ! final projection .......................................................

      if (this % project > 0) then
        call SetArray(p, ZERO)
        call PressureSolver(problem, flow_op, dt, u, p, w)
        call ProjectionStep(problem, flow_op, dt, p, u, w)
      end if

      ! pressure ...............................................................

      if (this % pressure > 0) then
        ! molecular + SGS viscosity
        if (allocated(nu)) then
          call problem % GetDiffusivity(flow_op%x, t, u, nu)
          if (flow_op % control % sgs > 0) then
            call Add_SGS_Diffusivity(flow_op, u, nu)
          end if
        end if
        ! recompute pressure using F_1 as workspace for F_v
        call TimeDerivative(problem, flow_op, t, u, nu=nu, chi=chi, F=F_1)
        call ComputePressure(problem, flow_op, t, F_1, p, w)
      end if

    end associate

    ! clean-up .................................................................

    if (present(standby)) then
      if (standby) return
    end if

    !$omp single
    if (allocated(u_i )) deallocate(u_i )
    if (allocated(u_1 )) deallocate(u_1 )
    if (allocated(F_0 )) deallocate(F_0 )
    if (allocated(F_1 )) deallocate(F_1 )
    if (allocated(dF  )) deallocate(dF  )
    if (allocated(F_d3)) deallocate(F_d3)
    if (allocated(w   )) deallocate(w   )
    if (allocated(nu  )) deallocate(nu  )
    !$omp end single

  end subroutine TimeStep

  !=============================================================================
  ! TimeIntegrator_TR_Options: type-bound procedures

  !-----------------------------------------------------------------------------
  !> MPI broadcasting of time-integrator options

  subroutine Bcast_TimeIntegrator_TR_Options(this, root, comm)
    class(TimeIntegrator_TR_Options), intent(inout) :: this
    integer,        intent(in) :: root !< rank of broadcast root
    type(MPI_Comm), intent(in) :: comm !< MPI communicator

    type(MPI_Request) :: request(4)
    type(MPI_Status)  :: stat(size(request))
    integer :: n

    ! broadcast options of parent class
    call this % TimeIntegratorOptions % Bcast(root, comm)

    n = 1
    call XMPI_Ibcast( this % i_max_p1, root, comm, request(n) );  n = n + 1
    call XMPI_Ibcast( this % i_max_u1, root, comm, request(n) );  n = n + 1
    call XMPI_Ibcast( this % i_max_p2, root, comm, request(n) );  n = n + 1
    call XMPI_Ibcast( this % i_max_u2, root, comm, request(n) )

    call MPI_Waitall(n, request, stat)

  end subroutine Bcast_TimeIntegrator_TR_Options

  !=============================================================================

end module CART__ISP_Flow__Time_Integrator__TR
