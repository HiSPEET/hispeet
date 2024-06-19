!> summary:  Explicit TVD Runge-Kutta method of order 3
!> author:   Joerg Stiller
!> date:     2024/02/04
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module CL__Time_Integrator__TVD_RK3__1D

  use, intrinsic :: ISO_Fortran_Env, only: OUTPUT_UNIT
  use, intrinsic :: IEEE_Arithmetic

  use Kind_Parameters, only: RNP
  use Constants
  use Array_Assignments

  use CL__Problem__1D
  use CL__Operator__1D
  use CL__Time_Integrator__1D

  implicit none
  private

  public :: CL_TimeIntegrator_TVD_RK3_1D
  public :: CL_TimeIntegrator_Options_TVD_RK3_1D

  !-----------------------------------------------------------------------------
  !> IMEX Euler method for Dahlquist equation

  type, extends(CL_TimeIntegrator_1D) :: CL_TimeIntegrator_TVD_RK3_1D
  contains
    procedure :: Init_CL_TimeIntegrator_TVD_RK3_1D
    procedure :: Show => Show_CL_TimeIntegrator_TVD_RK3_1D
    procedure :: TimeStep
  end type CL_TimeIntegrator_TVD_RK3_1D

  ! overloading the constructor
  interface CL_TimeIntegrator_TVD_RK3_1D
    module procedure New_CL_TimeIntegrator_TVD_RK3_1D
  end interface

  !-----------------------------------------------------------------------------
  !> Type for providing Euler time-integrator options (none, so far)

  type, extends(CL_TimeIntegrator_Options_1D) :: &
    CL_TimeIntegrator_Options_TVD_RK3_1D
  end type CL_TimeIntegrator_Options_TVD_RK3_1D

contains

  !-----------------------------------------------------------------------------
  !> Constructor for objects of type CL_TimeIntegrator_TVD_RK3_1D with options

  function New_CL_TimeIntegrator_TVD_RK3_1D(opt) result(this)
    class(CL_TimeIntegrator_Options_TVD_RK3_1D), intent(in) :: opt
    type(CL_TimeIntegrator_TVD_RK3_1D) :: this

    call Init_CL_TimeIntegrator_TVD_RK3_1D(this, opt)

  end function New_CL_TimeIntegrator_TVD_RK3_1D

  !-----------------------------------------------------------------------------
  !> Initialization of a Init_CL_TimeIntegrator_TVD_RK3_1D object

  subroutine Init_CL_TimeIntegrator_TVD_RK3_1D(this, opt)
    class(CL_TimeIntegrator_TVD_RK3_1D),         intent(inout) :: this
    class(CL_TimeIntegrator_Options_TVD_RK3_1D), intent(in)    :: opt

    ! intialize parent type
    call this % Init_CL_TimeIntegrator_1D(opt)
    this % name = 'TVD Runge-Kutta method of order 3'

  end subroutine Init_CL_TimeIntegrator_TVD_RK3_1D

  !-----------------------------------------------------------------------------
  !> Output of CL_TimeIntegrator_TVD_RK3_1D settings

  subroutine Show_CL_TimeIntegrator_TVD_RK3_1D(this, unit)
    class(CL_TimeIntegrator_TVD_RK3_1D), intent(in) :: this
    integer, optional, intent(in) :: unit  !< output unit

    integer :: io

    if (present(unit)) then
      io = unit
    else
      io = OUTPUT_UNIT
    end if

    ! show parent settings
    call this % Show_CL_TimeIntegrator_1D(unit)

  end subroutine Show_CL_TimeIntegrator_TVD_RK3_1D

  !-----------------------------------------------------------------------------
  !> Performs an IMEX Euler step

  subroutine TimeStep(this, cl_problem, cl_operator, dt, t_0, u_0, u)
    class(CL_TimeIntegrator_TVD_RK3_1D), intent(in) :: this
    class(CL_Problem_1D),  intent(in)    :: cl_problem
    class(CL_Operator_1D), intent(in)    :: cl_operator
    real(RNP),             intent(in)    :: dt          !< step size ∆t
    real(RNP),             intent(in)    :: t_0         !< initial time
    real(RNP), contiguous, intent(in)    :: u_0(0:,:,:) !< u(t₀)
    real(RNP), contiguous, intent(inout) :: u  (0:,:,:) !< u(t₀+∆t)

    real(RNP), allocatable, save :: r_c(:,:,:)
    real(RNP), allocatable, save :: r_d(:,:,:)
    real(RNP), allocatable, save :: f_s(:,:,:)
    real(RNP), allocatable, save :: u_1(:,:,:)
    real(RNP), allocatable, save :: u_2(:,:,:)
    real(RNP), allocatable, save :: bv(:,:)

    real(RNP), allocatable :: Me_inv(:)
    real(RNP) :: t_1, t_2, dt_1, dt_2, dt_3
    real(RNP) :: a, b
    integer   :: e, k

    associate( nc       => cl_problem  % nc       &
             , po       => cl_operator % eop % po &
             , ne       => cl_operator % ne       &
             , Me       => cl_operator % Me       &
             , activity => cl_operator % activity )

      !$omp master

      ! initialization .........................................................

      allocate(r_c, mold = u)
      allocate(r_d, mold = u)
      allocate(f_s, mold = u)
      allocate(u_1, mold = u)
      allocate(u_2, mold = u)
      allocate(bv(nc,2))

      allocate(Me_inv(0:po), source = 1/Me)

      t_1 = t_0 + dt
      t_2 = t_0 + dt/2

      dt_1 = dt
      dt_2 = dt/4
      dt_3 = 2*dt/3

      ! stage 1 ................................................................

      call cl_problem % GetBoundaryValues(t_0, bv)
      call cl_problem % GetSources(cl_operator, t_0, u_0, f_s)
      call cl_problem % GetConvectionTerm(cl_operator, bv, u_0, r_c)

      if (cl_problem % HasDiffusion()) then
        call cl_problem % GetDiffusionTerm(cl_operator, bv, u_0, r_d)
      else
        call SetArray(r_d, ZERO, multi=.true.)
      end if

      do k = 1, nc
      do e = 1, ne
        if (activity(e) < 1) cycle
        u_1(:,e,k) = u_0(:,e,k) &
                   + dt_1 * (Me_inv * (r_c(:,e,k) + r_d(:,e,k)) + f_s(:,e,k))
      end do
      end do

      if (cl_problem % limiting_scope > 1) then
        if (cl_problem % limiting_method == 1) then
          call cl_problem % MomentLimiter(cl_operator, u_1)
        end if
      end if

      ! stage 2 ................................................................

      call cl_problem % GetBoundaryValues(t_1, bv)
      call cl_problem % GetSources(cl_operator, t_1, u_1, f_s)
      call cl_problem % GetConvectionTerm(cl_operator, bv, u_1, r_c)

      if (cl_problem % HasDiffusion()) then
        call cl_problem % GetDiffusionTerm(cl_operator, bv, u_1, r_d)
      end if

      a = real(THREE / 4, RNP)
      b = real(ONE   / 4, RNP)

      do k = 1, nc
      do e = 1, ne
        if (activity(e) < 1) cycle
        u_2(:,e,k) = a * u_0(:,e,k) &
                   + b * u_1(:,e,k) &
                   + dt_2 * (Me_inv * (r_c(:,e,k) + r_d(:,e,k)) + f_s(:,e,k))
      end do
      end do

      if (cl_problem % limiting_scope > 1) then
        if (cl_problem % limiting_method == 1) then
          call cl_problem % MomentLimiter(cl_operator, u_2)
        end if
      end if

      ! stage 3 ................................................................

      call cl_problem % GetBoundaryValues(t_2, bv)
      call cl_problem % GetSources(cl_operator, t_2, u_2, f_s)
      call cl_problem % GetConvectionTerm(cl_operator, bv, u_2, r_c)

      if (cl_problem % HasDiffusion()) then
        call cl_problem % GetDiffusionTerm(cl_operator, bv, u_2, r_d)
      end if

      a = 1 * THIRD
      b = 2 * THIRD

      do k = 1, nc
      do e = 1, ne
        if (activity(e) < 1) cycle
        u(:,e,k) = a * u_0(:,e,k) &
                 + b * u_2(:,e,k) &
                 + dt_3 * (Me_inv * (r_c(:,e,k) + r_d(:,e,k)) + f_s(:,e,k))
      end do
      end do

      if (cl_problem % limiting_scope > 0) then
        if (cl_problem % limiting_method == 1) then
          call cl_problem % MomentLimiter(cl_operator, u)
        end if
      end if

      ! finalization ...........................................................

      deallocate(r_c, r_d, f_s, u_1, u_2, bv)

    end associate

  end subroutine TimeStep

  !=============================================================================

end module CL__Time_Integrator__TVD_RK3__1D
