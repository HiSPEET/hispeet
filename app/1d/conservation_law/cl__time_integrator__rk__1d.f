!> summary:  Runge-Kutta time integrators for 1D conservation laws
!> author:   Robin Fränzel, Joerg Stiller
!> date:     2023/06/06
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> @note
!> This method is essentially an implementation of the IMEX Runge-Kutta method,
!> but can also be run in fully explicit mode.
!===============================================================================

module CL__Time_Integrator__RK__1D

  use, intrinsic :: ISO_Fortran_Env, only: OUTPUT_UNIT
  use, intrinsic :: IEEE_Arithmetic

  use Kind_Parameters, only: RNP
  use Constants,       only: ZERO, ONE
  use Logging_Levels
  use Array_Assignments
  use IMEX_Runge_Kutta_Method

  use CL__Problem__1D
  use CL__Operator__1D
  use CL__Time_Integrator__1D

  implicit none
  private

  public :: CL_TimeIntegrator_RK_1D
  public :: CL_TimeIntegrator_Options_RK_1D

  !-----------------------------------------------------------------------------
  !> IMEX Euler method for Dahlquist equation

  type, extends(CL_TimeIntegrator_1D) :: CL_TimeIntegrator_RK_1D
    type(IMEX_RK_Method) :: imex_rk  !< IMEX Runge-Kutta method
  contains
    procedure :: Init_CL_TimeIntegrator_RK_1D
    procedure :: Show => Show_CL_TimeIntegrator_RK_1D
    procedure :: TimeStep
  end type CL_TimeIntegrator_RK_1D

  ! overloading the constructor
  interface CL_TimeIntegrator_RK_1D
    module procedure New_CL_TimeIntegrator_RK_1D
  end interface

  !-----------------------------------------------------------------------------
  !> Type for providing Euler time-integrator options (none, so far)

  type, extends(CL_TimeIntegrator_Options_1D) :: CL_TimeIntegrator_Options_RK_1D
    integer :: n_stage = 4  !< number of stages
    integer :: method  = 1  !< RK method selector, if more than one exist
  end type CL_TimeIntegrator_Options_RK_1D

contains

  !-----------------------------------------------------------------------------
  !> Constructor for objects of type CL_TimeIntegrator_RK_1D with options

  function New_CL_TimeIntegrator_RK_1D(opt) result(this)
    class(CL_TimeIntegrator_Options_RK_1D), intent(in) :: opt
    type(CL_TimeIntegrator_RK_1D) :: this

    call Init_CL_TimeIntegrator_RK_1D(this, opt)

  end function New_CL_TimeIntegrator_RK_1D

  !-----------------------------------------------------------------------------
  !> Initialization of a Init_CL_TimeIntegrator_RK_1D object

  subroutine Init_CL_TimeIntegrator_RK_1D(this, opt)
    class(CL_TimeIntegrator_RK_1D),         intent(inout) :: this
    class(CL_TimeIntegrator_Options_RK_1D), intent(in)    :: opt

    ! intialize parent type
    call this % Init_CL_TimeIntegrator_1D(opt)
    this % name = 'IMEX Runge-Kutta method'

    ! initialize RK method
    call this % imex_rk % Init_IMEX_RK_Method(opt % n_stage, opt % method)

  end subroutine Init_CL_TimeIntegrator_RK_1D

  !-----------------------------------------------------------------------------
  !> Output of CL_TimeIntegrator_RK_1D settings

  subroutine Show_CL_TimeIntegrator_RK_1D(this, unit)
    class(CL_TimeIntegrator_RK_1D), intent(in) :: this
    integer, optional, intent(in) :: unit  !< output unit

    integer :: io

    if (present(unit)) then
      io = unit
    else
      io = OUTPUT_UNIT
    end if

    ! show parent settings
    call this % Show_CL_TimeIntegrator_1D(unit)

    ! show IMEX RK settings
    call this % imex_rk % Show(unit)

  end subroutine Show_CL_TimeIntegrator_RK_1D

  !-----------------------------------------------------------------------------
  !> Performs an IMEX Runge-Kutta step

  subroutine TimeStep(this, cl_problem, cl_operator, dt, t_0, u_0, u)
    class(CL_TimeIntegrator_RK_1D), intent(in) :: this
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

    real(RNP), allocatable, save :: f_ex(:,:,:,:) ! explicit RHS
    real(RNP), allocatable, save :: f_im(:,:,:,:) ! implicit RHS

    real(RNP), allocatable :: Me_inv(:)
    real(RNP) :: d_ex, d_im, dt_a_ii, t_i
    integer   :: e, i, j, k

!### CHECK
if (log_level > 0) then
print '(99(G0,X))','TI RK: #0'
end if
!### CHECK END
    associate( a_ex     => this % imex_rk % a_ex    &
             , a_im     => this % imex_rk % a_im    &
             , b_ex     => this % imex_rk % b_ex    &
             , b_im     => this % imex_rk % b_im    &
             , c        => this % imex_rk % c       &
             , ns       => this % imex_rk % n_stage &
             , nc       => cl_problem  % nc         &
             , po       => cl_operator % eop % po   &
             , ne       => cl_operator % ne         &
             , Me       => cl_operator % Me         &
             , activity => cl_operator % activity   )

      !$omp master

      ! initialization .........................................................
!### CHECK
if (log_level > 1) then
print '(99(G0,X))','TI RK: #0.1 '
print '(99(G0,X))','TI RK: #0.1  minval(u_0(:,:,1)) =',minval(u_0(:,:,1))
print '(99(G0,X))','TI RK: #0.1  minval(u_0(:,:,2)) =',minval(u_0(:,:,2))
print '(99(G0,X))','TI RK: #0.1  minval(u_0(:,:,3)) =',minval(u_0(:,:,3))
end if
!### CHECK END

      allocate(r_c, mold = u)
      allocate(r_d, mold = u)
      allocate(f_s, mold = u)
      allocate(u_i, mold = u)
      allocate(bv(nc,2))
!### CHECK
if (log_level > 1) then
print '(99(G0,X))','TI RK: #0.2'
end if
!### CHECK END

      allocate(f_ex(0:po,ne,nc,ns))
      allocate(f_im, mold = f_ex)
!### CHECK
if (log_level > 1) then
print '(99(G0,X))','TI RK: #0.3'
end if
!### CHECK END

      allocate(Me_inv(0:po), source = 1/Me)
!### CHECK
if (log_level > 1) then
print '(99(G0,X))','TI RK: #0.4'
end if
!### CHECK END

      ! stage 1 ................................................................

      t_i = t_0

      call SetArray(u_i, u_0, multi = .true.)
!### CHECK
if (log_level > 1) then
print '(99(G0,X))','TI RK: #0.5'
print '(99(G0,X))','TI RK: #0.5  minval(u_i(:,:,1)) =',minval(u_i(:,:,1))
print '(99(G0,X))','TI RK: #0.5  minval(u_i(:,:,2)) =',minval(u_i(:,:,2))
print '(99(G0,X))','TI RK: #0.5  minval(u_i(:,:,3)) =',minval(u_i(:,:,3))
end if
!### CHECK END

      call cl_problem % GetBoundaryValues(t_i, bv)
!### CHECK
if (log_level > 1) then
print '(99(G0,X))','TI RK: #0.6'
end if
!### CHECK END
      call cl_problem % GetConvectionTerm(cl_operator, bv, u_i, r_c)
!### CHECK
if (log_level > 1) then
print '(99(G0,X))','TI RK: #0.7'
end if
!### CHECK END
      call cl_problem % GetHybridDiffusionTerm &
                            (cl_operator, 'T', ZERO, bv, u_0, u_i, r_d)
!### CHECK
if (log_level > 1) then
print '(99(G0,X))','TI RK: #0.8'
end if
!### CHECK END
      call cl_problem % GetSources(cl_operator, t_i, u_i, f_s)
!### CHECK
if (log_level > 0) then
print '(99(G0,X))','TI RK: #1  any(ieee_is_nan(r_c))  =',any(ieee_is_nan(r_c))
print '(99(G0,X))','TI RK: #1  any(ieee_is_nan(r_d))  =',any(ieee_is_nan(r_d))
print '(99(G0,X))','TI RK: #1  any(ieee_is_nan(f_s))  =',any(ieee_is_nan(f_s))
end if
!### CHECK END

      do k = 1, nc
      do e = 1, ne
        if (activity(e) > 0) then
          f_ex(:,e,k,1) = Me_inv * r_c(:,e,k)
          f_im(:,e,k,1) = Me_inv * r_d(:,e,k) + f_s(:,e,k)
        else
          f_ex(:,e,k,1) = 0
          f_im(:,e,k,1) = 0
        end if
      end do
      end do

      if (this % impl == 0) then
        ! explicit method: f_ex = f_ex + f_im, f_im = 0
        call MergeArrays(ONE, f_ex(:,:,:,1), ONE, f_im(:,:,:,1), multi=.true.)
        call SetArray(f_im(:,:,:,1), ZERO)
      end if
!### CHECK
if (log_level > 0) then
print '(99(G0,X))','TI RK: #1  any(ieee_is_nan(f_ex)) =',any(ieee_is_nan(f_ex))
print '(99(G0,X))','TI RK: #1  any(ieee_is_nan(f_im)) =',any(ieee_is_nan(f_im))
end if
!### CHECK END

      ! stages 2:ns ............................................................

      Stages: do i = 2, ns

        ! initialization . . . . . . . . . . . . . . . . . . . . . . . . . . . .

        t_i      =  t_0 + dt * c(i)
        dt_a_ii  =  dt * a_im(i,i)

        do k = 1, nc
        do e = 1, ne
          u_i(:,e,k) = u_0(:,e,k)  ! stage solution
          r_d(:,e,k) = u_0(:,e,k)  ! diffusion RHS
        end do
        end do

        ! lower stage contributions  . . . . . . . . . . . . . . . . . . . . . .

        do j = 1, i-1
        do k = 1, nc
        do e = 1, ne
          if (activity(e) > 0) then

            ! extrapolation of stage solution
            u_i(:,e,k) = u_i(:,e,k) &
                       + dt * a_ex(i,j) * (f_ex(:,e,k,j) + f_im(:,e,k,j))

            ! composition of diffusion RHS, used only in IMEX mode
            r_d(:,e,k) = r_d(:,e,k)                       &
                       + dt * ( a_ex(i,j) * f_ex(:,e,k,j) &
                              + a_im(i,j) * f_im(:,e,k,j) )
          end if
        end do
        end do
        end do

        ! implicit part  . . . . . . . . . . . . . . . . . . . . . . . . . . . .

        if (this % impl == 1) then

          call cl_problem % GetBoundaryValues(t_i, bv)
          call cl_problem % GetSources(cl_operator, t_i, u_i, f_s)
          call MergeArrays(ONE, r_d, dt_a_ii, f_s, multi=.true.)

          call cl_problem % DiffusionSolver( cl_operator, dt_a_ii, ZERO, bv    &
                                           , f      = r_d                      &
                                           , u_0    = u_0                      &
                                           , u      = u_i                      &
                                           , method = this % diffusion_method  &
                                           , i_max  = this % diffusion_i_max   &
                                           , r_red  = this % diffusion_r_red   &
                                           , r_max  = this % diffusion_r_max   )

        end if

        ! limiting . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

        if (cl_problem % limiting_method == 1 .and. &
            cl_problem % limiting_scope  == 2       ) then
          call cl_problem % MomentLimiter(cl_operator, u)
        end if

        ! RHS contributions  . . . . . . . . . . . . . . . . . . . . . . . . . .

        call cl_problem % GetConvectionTerm(cl_operator, bv, u_i, r_c)
        call cl_problem % GetHybridDiffusionTerm &
                              (cl_operator, 'T', ZERO, bv, u_0, u_i, r_d)

        do k = 1, nc
        do e = 1, ne
          if (activity(e) > 0) then
            f_ex(:,e,k,i) = Me_inv * r_c(:,e,k)
            f_im(:,e,k,i) = Me_inv * r_d(:,e,k) + f_s(:,e,k)
          else
            f_ex(:,e,k,i) = 0
            f_im(:,e,k,i) = 0
          end if
        end do
        end do

        if (this % impl == 0) then
          ! explicit method: f_ex = f_ex + f_im, f_im = 0
          call MergeArrays(ONE, f_ex(:,:,:,i), ONE, f_im(:,:,:,i), multi=.true.)
          call SetArray(f_im(:,:,:,i), ZERO)
        end if
!### CHECK
if (log_level > 0) then
print '(99(G0,X))','TI RK: #',i,' minval(u_i(:,:,1))     =',minval(u_i(:,:,1))
print '(99(G0,X))','TI RK: #',i,' minval(u_i(:,:,2))     =',minval(u_i(:,:,2))
print '(99(G0,X))','TI RK: #',i,' minval(u_i(:,:,3))     =',minval(u_i(:,:,3))
print '(99(G0,X))','TI RK: #',i,' any(ieee_is_nan(r_c))  =',any(ieee_is_nan(r_c))
print '(99(G0,X))','TI RK: #',i,' any(ieee_is_nan(r_d))  =',any(ieee_is_nan(r_d))
print '(99(G0,X))','TI RK: #',i,' any(ieee_is_nan(f_s))  =',any(ieee_is_nan(f_s))
print '(99(G0,X))','TI RK: #',i,' any(ieee_is_nan(f_ex)) =',any(ieee_is_nan(f_ex))
print '(99(G0,X))','TI RK: #',i,' any(ieee_is_nan(f_im)) =',any(ieee_is_nan(f_im))
end if
!### CHECK END

      end do Stages

      ! assembly ...............................................................

      do i = 1, ns
        d_ex = dt * (b_ex(i) - a_ex(ns,i))
        d_im = dt * (b_im(i) - a_im(ns,i))
        if (d_ex /= ZERO .or. d_im /= ZERO) cycle
        do k = 1, nc
        do e = 1, ne
          if (activity(e) > 0) then
            u_i(:,e,k) = u_i(:,e,k)            &
                       + d_ex * f_ex(:,e,k,i)  &
                       + d_im * f_im(:,e,k,i)
          end if
        end do
        end do
      end do

      call SetArray(u, u_i)

      ! limiting ...............................................................

      if (cl_problem % limiting_method == 1) then
        call cl_problem % MomentLimiter(cl_operator, u)
      end if

      ! finalization ...........................................................

      deallocate(r_c, r_d, f_s, u_i, bv, f_ex, f_im)

    end associate

  end subroutine TimeStep

  !=============================================================================

end module CL__Time_Integrator__RK__1D
