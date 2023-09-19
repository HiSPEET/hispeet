!> summary:  Euler time integrators for 1D conservation laws
!> author:   Robin Fränzel, Joerg Stiller
!> date:     2023/05/04
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module CL__Time_Integrator__Euler__1D

  use, intrinsic :: ISO_Fortran_Env, only: OUTPUT_UNIT

  use Kind_Parameters, only: RNP
  use Constants,       only: ZERO, ONE
  use Array_Assignments

  use CL__Problem__1D
  use CL__Operator__1D
  use CL__Time_Integrator__1D

  implicit none
  private

  public :: CL_TimeIntegrator_Euler_1D
  public :: CL_TimeIntegrator_Options_Euler_1D

  !-----------------------------------------------------------------------------
  !> IMEX Euler method for Dahlquist equation

  type, extends(CL_TimeIntegrator_1D) :: CL_TimeIntegrator_Euler_1D
  contains
    procedure :: Init_CL_TimeIntegrator_Euler_1D
    procedure :: Show => Show_CL_TimeIntegrator_Euler_1D
    procedure :: TimeStep
  end type CL_TimeIntegrator_Euler_1D

  ! overloading the constructor
  interface CL_TimeIntegrator_Euler_1D
    module procedure New_CL_TimeIntegrator_Euler_1D
  end interface

  !-----------------------------------------------------------------------------
  !> Type for providing Euler time-integrator options (none, so far)

  type, extends(CL_TimeIntegrator_Options_1D) :: &
    CL_TimeIntegrator_Options_Euler_1D
  end type CL_TimeIntegrator_Options_Euler_1D

contains

  !-----------------------------------------------------------------------------
  !> Constructor for objects of type CL_TimeIntegrator_Euler_1D with options

  function New_CL_TimeIntegrator_Euler_1D(opt) result(this)
    class(CL_TimeIntegrator_Options_Euler_1D), intent(in) :: opt
    type(CL_TimeIntegrator_Euler_1D) :: this

    call Init_CL_TimeIntegrator_Euler_1D(this, opt)

  end function New_CL_TimeIntegrator_Euler_1D

  !-----------------------------------------------------------------------------
  !> Initialization of a Init_CL_TimeIntegrator_Euler_1D object

  subroutine Init_CL_TimeIntegrator_Euler_1D(this, opt)
    class(CL_TimeIntegrator_Euler_1D),         intent(inout) :: this
    class(CL_TimeIntegrator_Options_Euler_1D), intent(in)    :: opt

    ! intialize parent type
    call this % Init_CL_TimeIntegrator_1D(opt)
    this % name = 'Euler method'

  end subroutine Init_CL_TimeIntegrator_Euler_1D

  !-----------------------------------------------------------------------------
  !> Output of CL_TimeIntegrator_Euler_1D settings

  subroutine Show_CL_TimeIntegrator_Euler_1D(this, unit)
    class(CL_TimeIntegrator_Euler_1D), intent(in) :: this
    integer, optional, intent(in) :: unit  !< output unit

    integer :: io

    if (present(unit)) then
      io = unit
    else
      io = OUTPUT_UNIT
    end if

    ! show parent settings
    call this % Show_CL_TimeIntegrator_1D(unit)

  end subroutine Show_CL_TimeIntegrator_Euler_1D

  !-----------------------------------------------------------------------------
  !> Performs an IMEX Euler step

  subroutine TimeStep(this, cl_problem, cl_operator, dt, t_0, u_0, u)
    class(CL_TimeIntegrator_Euler_1D), intent(in) :: this
    class(CL_Problem_1D),  intent(in)    :: cl_problem
    class(CL_Operator_1D), intent(in)    :: cl_operator
    real(RNP),             intent(in)    :: dt          !< step size ∆t
    real(RNP),             intent(in)    :: t_0         !< initial time
    real(RNP), contiguous, intent(in)    :: u_0(0:,:,:) !< u(t₀)
    real(RNP), contiguous, intent(inout) :: u  (0:,:,:) !< u(t₀+∆t)

    real(RNP), allocatable, save :: r_c(:,:,:)
    real(RNP), allocatable, save :: r_d(:,:,:)
    real(RNP), allocatable, save :: f_s(:,:,:)
    real(RNP), allocatable, save :: u_i(:,:,:)
    real(RNP), allocatable, save :: bv(:,:)

    real(RNP), allocatable :: Me_inv(:)
    real(RNP) :: t
    integer   :: e, k

    associate( nc       => cl_problem  % nc       &
             , po       => cl_operator % eop % po &
             , ne       => cl_operator % ne       &
             , Me       => cl_operator % Me       &
             , activity => cl_operator % activity )

      !$omp master

      ! initialization .........................................................

      t = t_0 + dt

      allocate(r_c, mold = u)
      allocate(r_d, mold = u)
      allocate(f_s, mold = u)
      allocate(u_i, mold = u)
      allocate(bv(nc,2))

      allocate(Me_inv(0:po), source = 1/Me)

      select case(this%impl)

      case(0)

        ! explicit Euler step ..................................................

        call cl_problem % GetBoundaryValues(t_0, bv)
        call cl_problem % GetConvectionTerm(cl_operator, bv, u_0, r_c)
        call cl_problem % GetDiffusionTerm(cl_operator, bv, u_0, r_d)
        call cl_problem % GetSources(cl_operator, t_0, u_0, f_s)

        do k = 1, nc
        do e = 1, ne
          if (activity(e) > 0) then
            u(:,e,k) = u_0(:,e,k) &
                     + dt * (Me_inv * (r_c(:,e,k) + r_d(:,e,k)) + f_s(:,e,k))
          else
            u(:,e,k) = u_0(:,e,k)
          end if
        end do
        end do

      case(1)

        ! IMEX Euler step ......................................................

        call cl_problem % GetBoundaryValues(t_0, bv)
        call cl_problem % GetConvectionTerm(cl_operator, bv, u_0, r_c)

        call cl_problem % GetBoundaryValues(t, bv)
        call cl_problem % GetSources(cl_operator, t, u_0, f_s)

        ! intermediate solution
        do k = 1, nc
        do e = 1, ne
          if (activity(e) > 0) then
            u_i(:,e,k) = u_0(:,e,k) + dt * (Me_inv * r_c(:,e,k) + f_s(:,e,k))
          else
            u_i(:,e,k) = u_0(:,e,k)
          end if
        end do
        end do

        ! implicit diffusion step
        call cl_problem % DiffusionSolver( cl_operator, dt, ZERO, bv         &
                                         , u_i, u_0, u                       &
                                         , method = this % diffusion_method  &
                                         , i_max  = this % diffusion_i_max   &
                                         , r_red  = this % diffusion_r_red   &
                                         , r_max  = this % diffusion_r_max   )

      end select

      ! finalization ...........................................................

      deallocate(r_c, r_d, f_s, u_i, bv)

    end associate

  end subroutine TimeStep

  !=============================================================================

end module CL__Time_Integrator__Euler__1D
