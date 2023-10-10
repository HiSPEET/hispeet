program Conservation_Law_ML
  use Kind_Parameters
  use Constants

  use CL__Operator__1D
  use CL__Problem__1D
  use CL__Problem__Burgers__Moving_Front__1D
  use CL__Time_Integrator__1D
  use CL__Time_Integrator__Euler__1D
  use CL__Time_Integrator__ISD1__1D
  use CL__SDC__Method__1D
  use CL__SDC__Method__ISD1__1D
  use CL__MLSDC__1D
  use CL__MLSDC__Level__1D
  use CL__MLSDC__Variable__1D

  implicit none

  class(CL_TimeIntegrator_Options_1D), allocatable :: opt_pre
  class(CL_SDC_Options_1D)           , allocatable :: opt_sdc
  class(CL_Problem_1D)               , allocatable :: cl_problem

  type(CL_MLSDC_1D)          :: mlsdc
  type(CL_MLSDC_Options_1D)  :: opt
  type(CL_MLSDC_Variable_1D) :: u

  real(RNP) :: t_0     = 0.25  ! slab start time
  real(RNP) :: t_1     = 0.75  ! slab end time
  integer   :: n_level = 2     ! number of space-time levels
  integer   :: l

  ! options
  opt = CL_MLSDC_Options_1D(n_level)
  do l = 1, n_level
    opt % n_space(l) = 2 ** (l-1)
    opt % p_space(l) = 5
    opt % q_conv(l)  = (opt % p_space(l) * 3 + 1) / 2
    opt % n_time(l)  = 1
    opt % p_time(l)  = 2 ** (l-1)
  end do
  print '(A,99I5)', 'n_level = ', opt % n_level
  print '(A,99I5)', 'n_space = ', opt % n_space
  print '(A,99I5)', 'p_space = ', opt % p_space
  print '(A,99I5)', 'q_conv  = ', opt % q_conv
  print '(A,99I5)', 'n_time  = ', opt % n_time
  print '(A,99I5)', 'p_time  = ', opt % p_time

  ! problem
  allocate(CL_Problem_Burgers_MovingFront_1D :: cl_problem)
  call cl_problem % SetProblem()

  ! predictor options (ISD1, so far)
  allocate(CL_TimeIntegrator_Options_ISD1_1D :: opt_pre)
  opt_pre % impl             = 1
  opt_pre % diffusion_method = 4
  opt_pre % diffusion_i_max  = 100

  ! SDC options (ISD1, so far)
  allocate(CL_SDC_Options_ISD1_1D :: opt_sdc)
  opt_sdc % point_set        = 2
  opt_sdc % diffusion_method = 4
  opt_sdc % diffusion_i_max  = 100

  ! MLSDC data structure
  mlsdc = CL_MLSDC_1D(opt, opt_pre, opt_sdc, cl_problem)

  ! MLSDC variable
  u = CL_MLSDC_Variable_1D(mlsdc)

  ! set variable
  do l = 1, n_level
    call GetExactSolution(mlsdc%level(l), t_0, t_1, u%level(l)%val)
  end do

  ! evaluate, plot etc.

contains

  !---------------------------------------------------------------------------
  !> Gets the exact solution of for a given mesh level

  subroutine GetExactSolution(level, t_0, t_1, u)
    class(CL_MLSDC_Level_1D), intent(in)  :: level !< space-time level
    real(RNP),                intent(in)  :: t_0   !< initial time
    real(RNP),                intent(in)  :: t_1   !< final time
    real(RNP), contiguous,    intent(out) :: u(0:,:,:,0:,:)

    real(RNP), allocatable :: t(:)
    real(RNP) :: dt
    integer   :: nt, pt
    integer   :: i, j

    pt = ubound(u,4)
    nt = ubound(u,5)
    dt = (t_1 - t_0) / nt

    allocate(t(0:pt))

    associate( cl_problem  => level % cl_problem  &
             , cl_operator => level % cl_operator &
             , cl_sdc      => level % cl_sdc      )

      do j = 1, nt
        t(0:) = cl_sdc % IntermediateTimes(t_0 + (j-1)*dt, dt)
        do i = 0, pt
          call cl_problem % GetExactSolution(cl_operator, t(i), u(:,:,:,i,j))
        end do
      end do

    end associate

  end subroutine GetExactSolution

  !=============================================================================

end program Conservation_Law_ML
