!> summary:  Streamline-diffusion method of order 1 for 1D conservation laws
!> author:   Joerg Stiller
!> date:     2023/05/05
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module CL__Time_Integrator__ISD1__1D

  use, intrinsic :: ISO_Fortran_Env, only: OUTPUT_UNIT

  use Kind_Parameters, only: RNP
  use Constants,       only: ZERO, ONE
  use Array_Assignments

  use CL__Problem__1D
  use CL__Operator__1D
  use CL__Time_Integrator__1D

  implicit none
  private

  public :: CL_TimeIntegrator_ISD1_1D
  public :: CL_TimeIntegrator_Options_ISD1_1D

  !-----------------------------------------------------------------------------
  !> IMEX ISD1 method for Dahlquist equation

  type, extends(CL_TimeIntegrator_1D) :: CL_TimeIntegrator_ISD1_1D
  contains
    procedure :: Init_CL_TimeIntegrator_ISD1_1D
    procedure :: Show => Show_CL_TimeIntegrator_ISD1_1D
    procedure :: TimeStep
  end type CL_TimeIntegrator_ISD1_1D

  ! overloading the constructor
  interface CL_TimeIntegrator_ISD1_1D
    module procedure New_CL_TimeIntegrator_ISD1_1D
  end interface

  !-----------------------------------------------------------------------------
  !> Type for providing ISD1 time-integrator options (none, so far)

  type, extends(CL_TimeIntegrator_Options_1D) :: &
    CL_TimeIntegrator_Options_ISD1_1D
  end type CL_TimeIntegrator_Options_ISD1_1D

contains

  !-----------------------------------------------------------------------------
  !> Constructor for objects of type CL_TimeIntegrator_ISD1_1D with options

  function New_CL_TimeIntegrator_ISD1_1D(opt) result(this)
    class(CL_TimeIntegrator_Options_ISD1_1D), optional, intent(in) :: opt
    type(CL_TimeIntegrator_ISD1_1D) :: this

    call Init_CL_TimeIntegrator_ISD1_1D(this, opt)

  end function New_CL_TimeIntegrator_ISD1_1D

  !-----------------------------------------------------------------------------
  !> Initialization of a Init_CL_TimeIntegrator_ISD1_1D object

  subroutine Init_CL_TimeIntegrator_ISD1_1D(this, opt)
    class(CL_TimeIntegrator_ISD1_1D),                   intent(inout) :: this
    class(CL_TimeIntegrator_Options_ISD1_1D), optional, intent(in)    :: opt

    ! intialize parent type
    call this % Init_CL_TimeIntegrator_1D(opt)
    this % name = 'Streamline-diffusion method of order 1'

  end subroutine Init_CL_TimeIntegrator_ISD1_1D

  !-----------------------------------------------------------------------------
  !> Output of CL_TimeIntegrator_ISD1_1D settings

  subroutine Show_CL_TimeIntegrator_ISD1_1D(this, unit)
    class(CL_TimeIntegrator_ISD1_1D), intent(in) :: this
    integer, optional, intent(in) :: unit  !< output unit

    integer :: io

    if (present(unit)) then
      io = unit
    else
      io = OUTPUT_UNIT
    end if

    ! show parent settings
    call this % Show_CL_TimeIntegrator_1D(unit)

    write(io,'(2X,A,T15,G0)')  'impl:', this % impl

  end subroutine Show_CL_TimeIntegrator_ISD1_1D

  !-----------------------------------------------------------------------------
  !> Performs an IMEX ISD1 step: u₁ = u₀ + ∆t (iλᵢ u₀ + λᵣ u₁)

  subroutine TimeStep(this, cl_problem, cl_operator, dt, t_0, u_0, u)
    class(CL_TimeIntegrator_ISD1_1D), intent(inout) :: this
    class(CL_Problem_1D),  intent(in) :: cl_problem
    class(CL_Operator_1D), intent(in) :: cl_operator
    real(RNP), intent(in)    :: dt          !< step size ∆t
    real(RNP), intent(in)    :: t_0         !< initial time
    real(RNP), intent(in)    :: u_0(0:,:,:) !< u(t₀)
    real(RNP), intent(inout) :: u  (0:,:,:) !< u(t₀+∆t)

    real(RNP), allocatable, save :: r_c(:,:,:)
    real(RNP), allocatable, save :: r_d(:,:,:)
    real(RNP), allocatable, save :: r_sd(:,:,:)
    real(RNP), allocatable, save :: f_s(:,:,:)
    real(RNP), allocatable, save :: u_i(:,:,:)
    real(RNP), allocatable, save :: bv(:,:)

    real(RNP), save :: t

    real(RNP), allocatable :: Me_inv(:)
    integer :: e, k

    associate( nc   => cl_problem  % nc       &
             , bc   => cl_problem  % bc       &
             , po   => cl_operator % eop % po &
             , ne   => cl_operator % ne       &
             , Me   => cl_operator % Me       &
             , mask => cl_operator % mask     )

      !$omp master

      ! initialization .........................................................

      t = t_0 + dt

      allocate(r_c , mold = u)
      allocate(r_d , mold = u)
      allocate(r_sd, mold = u)
      allocate(f_s , mold = u)
      allocate(u_i , mold = u)
      allocate(bv(nc,2))

      allocate(Me_inv(0:po), source = 1/Me)

      select case(this%impl)

      case(0)

        ! explicit ISD1 step ..................................................

        call cl_problem % GetBoundaryValues(t_0, bv)
        call cl_problem % GetConvectionTerm(cl_operator, bv, u_0, r_c)
        call cl_problem % GetDiffusionTerm(cl_operator, bv, u_0, r_d)
        call cl_problem % GetSDTerm(cl_operator, dt, bv, u_0, r_sd)
        call cl_problem % GetSources(cl_operator, t_0, u_0, f_s)

        do k = 1, nc
        do e = 1, ne
          if (mask(e)) then
            u(:,e,k) = u_0(:,e,k) &
                     + dt * ( Me_inv * (r_c(:,e,k) + r_d(:,e,k) + r_sd(:,e,k)) &
                            + f_s(:,e,k))
          else
            u(:,e,k) = u_0(:,e,k)
          end if
        end do
        end do

      case(1)

        ! IMEX ISD1 step ......................................................

        call cl_problem % GetBoundaryValues(t_0, bv)
        call cl_problem % GetConvectionTerm(cl_operator, bv, u_0, r_c)

        call cl_problem % GetBoundaryValues(t, bv)
        call cl_problem % GetSources(cl_operator, t, u_0, f_s)

        ! intermediate solution
        do k = 1, nc
        do e = 1, ne
          if (mask(e)) then
            u_i(:,e,k) = u_0(:,e,k) + dt * (Me_inv * r_c(:,e,k) + f_s(:,e,k))
          else
            u_i(:,e,k) = u_0(:,e,k)
          end if
        end do
        end do

        ! implicit diffusion step
        call cl_problem % DiffusionSolver( cl_operator, dt, dt, bv, u_i, u &
                                         , method = 4, i_max = 20          )

      end select

      ! finalization ...........................................................

      deallocate(r_c, r_d, r_sd, f_s, u_i, bv)

    end associate

  end subroutine TimeStep

  !=============================================================================

end module CL__Time_Integrator__ISD1__1D
