module DQ__Time_Integrator__RK

  use, intrinsic :: ISO_Fortran_Env, only: OUTPUT_UNIT

  use Kind_Parameters, only: RNP
  use Constants,       only: ZERO, ONE
  use IMEX_Runge_Kutta_Method
  use DQ__Time_Integrator

  implicit none
  private

  public :: DQ_TimeIntegrator_RK
  public :: DQ_TimeIntegrator_Options_RK

  !-----------------------------------------------------------------------------
  !> IMEX Runge-Kutta method for Dahlquist equation

  type, extends(DQ_TimeIntegrator) :: DQ_TimeIntegrator_RK
    type(IMEX_RK_Method) :: imex_rk  !< IMEX Runge-Kutta method
  contains
    procedure :: Init_DQ_TimeIntegrator_RK
    procedure :: Show => Show_DQ_TimeIntegrator_RK
    procedure :: TimeStep
  end type DQ_TimeIntegrator_RK

  ! overloading the constructor
  interface DQ_TimeIntegrator_RK
    module procedure New_DQ_TimeIntegrator_RK
  end interface

  !-----------------------------------------------------------------------------
  !> Type for providing RK time-integrator options

  type, extends(DQ_TimeIntegratorOptions) :: DQ_TimeIntegrator_Options_RK
    integer :: n_stage = 4  !< number of stages
    integer :: method  = 1  !< RK method selector, if more than one exist
  end type DQ_TimeIntegrator_Options_RK

  !-----------------------------------------------------------------------------
  !> Imaginary unit

  complex(RNP), parameter :: im1 = (ZERO, ONE)

contains

  !-----------------------------------------------------------------------------
  !> Constructor for objects of type DQ_TimeIntegrator_RK with options

  function New_DQ_TimeIntegrator_RK(opt) result(this)
    class(DQ_TimeIntegrator_Options_RK), optional, intent(in) :: opt
    type(DQ_TimeIntegrator_RK) :: this

    call Init_DQ_TimeIntegrator_RK(this, opt)

  end function New_DQ_TimeIntegrator_RK

  !-----------------------------------------------------------------------------
  !> Initialization of a Init_DQ_TimeIntegrator_RK object

  subroutine Init_DQ_TimeIntegrator_RK(this, opt)
    class(DQ_TimeIntegrator_RK),                   intent(inout) :: this
    class(DQ_TimeIntegrator_Options_RK), optional, intent(in)    :: opt

    ! intialize parent type
    call this % Init_DQ_TimeIntegrator(opt)
    this % name = 'IMEX Runge-Kutta method'

    ! initialize RK method
    call this % imex_rk % Init_IMEX_RK_Method(opt % n_stage, opt % method)

  end subroutine Init_DQ_TimeIntegrator_RK

  !-----------------------------------------------------------------------------
  !> Output of DQ_TimeIntegrator_RK settings

  subroutine Show_DQ_TimeIntegrator_RK(this, unit)
    class(DQ_TimeIntegrator_RK), intent(in) :: this
    integer, optional, intent(in) :: unit  !< output unit

    integer :: io

    if (present(unit)) then
      io = unit
    else
      io = OUTPUT_UNIT
    end if

    ! show parent settings
    call this % Show_DQ_TimeIntegrator(unit)

    write(io,'(A)')           'DQ_TimeIntegrator_RungeKutta settings'
    write(io,'(A,/)')         repeat('-',80)

    ! show IMEX RK settings
    call this % imex_rk % Show(unit)

  end subroutine Show_DQ_TimeIntegrator_RK

  !-----------------------------------------------------------------------------
  !> Performs an RK time step

  subroutine TimeStep(this, lambda, dt, u)
    class(DQ_TimeIntegrator_RK), intent(inout) :: this
    real(RNP),    intent(in)    :: dt      !< step size ∆t
    complex(RNP), intent(in)    :: lambda  !< ...
    complex(RNP), intent(inout) :: u       !< u(t) → u(t+ ∆t)

    select case(this%impl)
    case(0)
      call TimeStep_EX(this, lambda, dt, u)
    case(1)
      call TimeStep_IM(this, lambda, dt, u)
    case default
      call TimeStep_IMEX(this, lambda, dt, u)
    end select

  end subroutine TimeStep

  !-----------------------------------------------------------------------------
  !> Performs an IMEX RK step

  subroutine TimeStep_IMEX(this, lambda, dt, u)
    class(DQ_TimeIntegrator_RK), intent(inout) :: this
    real(RNP),    intent(in)    :: dt      !< step size ∆t
    complex(RNP), intent(in)    :: lambda  !< ...
    complex(RNP), intent(inout) :: u       !< u(t) → u(t+ ∆t)

    complex(RNP), allocatable :: u_s(:)  ! stage solutions
    integer :: i, j

    associate( a_im    => this % imex_rk % a_im    &
             , a_ex    => this % imex_rk % a_ex    &
             , b_im    => this % imex_rk % b_im    &
             , b_ex    => this % imex_rk % b_ex    &
             , ns      => this % imex_rk % n_stage )

      ! initialization .........................................................

      allocate(u_s(ns))

      ! stage 1 ................................................................

      u_s(1) = u

      ! stages 2:ns ............................................................

      do i = 2, ns
        u_s(i) = u
        do j = 1, i-1
          u_s(i) = u_s(i)                                             &
                 + dt * a_ex(i,j) * (ZERO, ONE) * lambda%im * u_s(j)  &
                 + dt * a_im(i,j) *               lambda%re * u_s(j)
        end do
        u_s(i) = u_s(i) / (ONE - dt * a_im(i,i) * lambda%re)
      end do

      ! assembly ...............................................................

      do i = 1, ns
        u = u + dt * b_ex(i) * (ZERO, ONE) * lambda%im * u_s(i)  &
              + dt * b_im(i) *               lambda%re * u_s(i)
      end do

    end associate

  end subroutine TimeStep_IMEX

  !-----------------------------------------------------------------------------
  !> Performs an implicit RK step

  subroutine TimeStep_IM(this, lambda, dt, u)
    class(DQ_TimeIntegrator_RK), intent(inout) :: this
    real(RNP),    intent(in)    :: dt      !< step size ∆t
    complex(RNP), intent(in)    :: lambda  !< ...
    complex(RNP), intent(inout) :: u       !< u(t) → u(t+ ∆t)

    complex(RNP), allocatable :: u_s(:)  ! stage solutions
    integer :: i, j

    associate( a_im  => this % imex_rk % a_im    &
             , b_im  => this % imex_rk % b_im    &
             , ns    => this % imex_rk % n_stage )

      ! initialization .........................................................

      allocate(u_s(ns))

      ! stage 1 ................................................................

      u_s(1) = u

      ! stages 2:ns ............................................................

      do i = 2, ns
        u_s(i) = u
        do j = 1, i-1
          u_s(i) = u_s(i) + dt * a_im(i,j) * lambda * u_s(j)
        end do
        u_s(i) = u_s(i) / (ONE - dt * a_im(i,i) * lambda)
      end do

      ! assembly ...............................................................

      do i = 1, ns
        u = u + dt * b_im(i) * lambda * u_s(i)
      end do

    end associate

  end subroutine TimeStep_IM

  !-----------------------------------------------------------------------------
  !> Performs a single explicit RK step

  subroutine TimeStep_EX(this, lambda, dt, u)
    class(DQ_TimeIntegrator_RK), intent(inout) :: this
    real(RNP),    intent(in)    :: dt      !< step size ∆t
    complex(RNP), intent(in)    :: lambda  !< ...
    complex(RNP), intent(inout) :: u       !< u(t) → u(t+ ∆t)

    complex(RNP), allocatable :: u_s(:)  ! stage solutions
    integer :: i, j

    associate( a_ex  => this % imex_rk % a_ex    &
             , b_ex  => this % imex_rk % b_ex    &
             , ns    => this % imex_rk % n_stage )

      ! initialization .........................................................

      allocate(u_s(ns))

      ! stage 1 ................................................................

      u_s(1) = u

      ! stages 2:ns ............................................................

      do i = 2, ns
        u_s(i) = u
        do j = 1, i-1
          u_s(i) = u_s(i) + dt * a_ex(i,j) * lambda * u_s(j)
        end do
      end do

      ! assembly ...............................................................

      do i = 1, ns
        u = u + dt * b_ex(i) * lambda * u_s(i)
      end do

    end associate

  end subroutine TimeStep_EX

  !=============================================================================

end module DQ__Time_Integrator__RK
