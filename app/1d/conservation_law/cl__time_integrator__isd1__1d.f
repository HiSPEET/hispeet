!> summary:  Streamline-diffusion method of order 1 for 1D conservation laws
!> author:   Joerg Stiller
!> date:     2023/05/05
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module CL__Time_Integrator__ISD1__1D

  use, intrinsic :: ISO_Fortran_Env, only: OUTPUT_UNIT

  use Kind_Parameters, only: RNP
  use Constants,       only: ZERO, ONE, HALF
  use Array_Assignments

  use CL__Problem__1D
  use CL__Operator__1D
  use CL__Time_Integrator__1D

  implicit none
  private

  public :: CL_TimeIntegrator_ISD1_1D
  public :: CL_TimeIntegrator_Options_ISD1_1D

  !-----------------------------------------------------------------------------
  !> IMEX ISD1 method for 1D conservation laws

  type, extends(CL_TimeIntegrator_1D) :: CL_TimeIntegrator_ISD1_1D
    integer   :: n_stage           !< number of stages, ignored with `imex_mode=0`
    integer   :: diffusion_i_max_1 !< max number of iterations in stage 1
    integer   :: diffusion_i_max_2 !< max number of iterations in stage 2
    real(RNP) :: sigma             !< ISD scaling factor
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
    integer   :: n_stage = 2            !< number of stages
    integer   :: diffusion_i_max_1 = -1 !< max number of iterations in stage 1
    integer   :: diffusion_i_max_2 = -1 !< max number of iterations in stage 2
    real(RNP) :: sigma = HALF           !< ISD scaling factor
  end type CL_TimeIntegrator_Options_ISD1_1D

contains

  !-----------------------------------------------------------------------------
  !> Constructor for objects of type CL_TimeIntegrator_ISD1_1D with options

  function New_CL_TimeIntegrator_ISD1_1D(opt) result(this)
    class(CL_TimeIntegrator_Options_ISD1_1D), intent(in) :: opt
    type(CL_TimeIntegrator_ISD1_1D) :: this

    call Init_CL_TimeIntegrator_ISD1_1D(this, opt)

  end function New_CL_TimeIntegrator_ISD1_1D

  !-----------------------------------------------------------------------------
  !> Initialization of a Init_CL_TimeIntegrator_ISD1_1D object

  subroutine Init_CL_TimeIntegrator_ISD1_1D(this, opt)
    class(CL_TimeIntegrator_ISD1_1D),         intent(inout) :: this
    class(CL_TimeIntegrator_Options_ISD1_1D), intent(in)    :: opt

    ! intialize parent type
    call this % Init_CL_TimeIntegrator_1D(opt)

    if (this % imex_mode > 0) then
      this % n_stage = opt % n_stage
    else
      this % n_stage = 1
    end if

    if (opt % diffusion_i_max_1 >= 0) then
      this % diffusion_i_max_1 = opt % diffusion_i_max_1
    else
      this % diffusion_i_max_1 = opt % diffusion_i_max
    end if

    if (opt % diffusion_i_max_2 >= 0) then
      this % diffusion_i_max_2 = opt % diffusion_i_max_2
    else
      this % diffusion_i_max_2 = opt % diffusion_i_max
    end if

    this % sigma = opt % sigma

    if (this%n_stage == 2) then
      this % name = 'Streamline-diffusion method of order 1 with two stages'
    else
      this % name = 'Streamline-diffusion method of order 1 with one stage'
    end if

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

    write(io,*)
    write(io,'(2X,A,T22,I0)') 'diffusion_i_max_1:', this % diffusion_i_max_1
    write(io,'(2X,A,T22,I0)') 'diffusion_i_max_2:', this % diffusion_i_max_2
    write(io,'(2X,A,T21,ES12.5)') 'sigma:', this % sigma

  end subroutine Show_CL_TimeIntegrator_ISD1_1D

  !-----------------------------------------------------------------------------
  !> Performs an IMEX ISD1 step: u₁ = u₀ + ∆t (iλᵢ u₀ + λᵣ u₁)

  subroutine TimeStep(this, cl_problem, cl_operator, dt, t_0, u_0, u)
    class(CL_TimeIntegrator_ISD1_1D), intent(in) :: this
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
    real(RNP) :: t, theta
    integer   :: e, k
    integer   :: i_max_1, i_max_2

    associate( nc       => cl_problem  % nc       &
             , eop      => cl_operator % eop      &
             , dx       => cl_operator % dx       &
             , po       => cl_operator % eop % po &
             , ne       => cl_operator % ne       &
             , activity => cl_operator % activity )

      !$omp master

      ! initialization .........................................................

      t = t_0 + dt
      theta = 2 * this%sigma * dt

      allocate(r_c , mold = u)
      allocate(r_d , mold = u)
      allocate(f_s , mold = u)
      allocate(u_i , mold = u)
      allocate(bv(nc,2))

      if (.not. cl_problem % HasDiffusion()) then
        call SetArray(r_d, ZERO)
      end if

      allocate(Me_inv(0:po), source = ONE/(dx/2 * eop%w))

      select case(this%imex_mode)

      case(0)

        ! explicit ISD1 step ..................................................

        call cl_problem % GetBoundaryValues(t_0, bv)
        call cl_problem % GetConvectionTerm(cl_operator, bv, u_0, r_c)
        call cl_problem % GetHybridDiffusionTerm &
                              (cl_operator, 'T', theta, bv, u_0, u_0, r_d)
        call cl_problem % GetSources(cl_operator, t_0, u_0, f_s)

        do k = 1, nc
        do e = 1, ne
          if (activity(e) > 0) then
            u(:,e,k) = u_0(:,e,k) &
                     + dt * (Me_inv * (r_c(:,e,k) + r_d(:,e,k)) + f_s(:,e,k))
          end if
        end do
        end do

      case(1:2)

        ! semi-implicit ISD1, stage 1 ..........................................

        call cl_problem % GetBoundaryValues(t_0, bv)
        call cl_problem % GetConvectionTerm(cl_operator, bv, u_0, r_c)

        if (this % imex_mode == 2) then
          call cl_problem % GetHybridDiffusionTerm &
                                (cl_operator, 'T', theta, bv, u_0, u_0, r_d)
        else
          call SetArray(r_d, ZERO, multi=.true.)
        end if

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

        ! start values
        if (this%imex_mode == 2) then
          do k = 1, nc
          do e = 1, ne
            if (activity(e) > 0) then
              u(:,e,k) = u_i(:,e,k) + dt *  Me_inv * r_d(:,e,k)
            end if
          end do
          end do
        else
          call SetArray(u, u_0, multi=.true.)
        end if

        ! implicit diffusion step
        call cl_problem % DiffusionSolver( cl_operator, dt, theta, bv        &
                                         , f      = u_i                      &
                                         , u_0    = u_0                      &
                                         , u      = u                        &
                                         , method = this % diffusion_method  &
                                         , i_max  = this % diffusion_i_max_1 &
                                         , r_red  = this % diffusion_r_red   &
                                         , r_max  = this % diffusion_r_max   )

        if (this % n_stage == 2) then

          ! semi-implicit ISD1, stage 2 ........................................

          if (cl_problem % limiting_method == 1 .and. &
              cl_problem % limiting_scope  == 2       ) then
            call cl_problem % MomentLimiter(cl_operator, u)
          end if

          call cl_problem % GetConvectionTerm(cl_operator, bv, u, r_c)

          if (this % imex_mode == 2) then
            call cl_problem % GetHybridDiffusionTerm &
                                  (cl_operator, 'T', theta, bv, u_0, u, r_d)
          end if

          ! intermediate solution
          do k = 1, nc
          do e = 1, ne
            if (activity(e) > 0) then
              u_i(:,e,k) = u_0(:,e,k) + dt * (Me_inv * r_c(:,e,k) + f_s(:,e,k))
            end if
          end do
          end do

          ! start values
          if (this%imex_mode == 2) then
            do k = 1, nc
            do e = 1, ne
              if (activity(e) > 0) then
                u(:,e,k) = u_i(:,e,k) + dt *  Me_inv * r_d(:,e,k)
              end if
            end do
            end do
          end if

          ! implicit diffusion step
          call cl_problem % DiffusionSolver( cl_operator, dt, theta, bv        &
                                           , f      = u_i                      &
                                           , u_0    = u_0                      &
                                           , u      = u                        &
                                           , method = this % diffusion_method  &
                                           , i_max  = this % diffusion_i_max_2 &
                                           , r_red  = this % diffusion_r_red   &
                                           , r_max  = this % diffusion_r_max   )

        end if

      end select

      ! limiting ...............................................................

      if (cl_problem % limiting_method == 1) then
        call cl_problem % MomentLimiter(cl_operator, u)
      end if

      ! finalization ...........................................................

      deallocate(r_c, r_d, f_s, u_i, bv)

    end associate

  end subroutine TimeStep

  !=============================================================================

end module CL__Time_Integrator__ISD1__1D
