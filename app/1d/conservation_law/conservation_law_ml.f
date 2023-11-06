program Conservation_Law_ML
  use Kind_Parameters
  use Constants
  use Execution_Control

  use CL__Operator__1D
  use CL__Problem__1D
  use CL__Problem__Convection_Diffusion__Wave_Package__1D
  use CL__Problem__Burgers__Wave_Package__1D
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

  ! declarations: control ......................................................

  character(len=80) :: case_name = 'conservation_law_ml'
  character(len=80) :: case_file

  ! declarations: problem ......................................................

  class(CL_Problem_1D), allocatable :: cl_problem

  character(len=80) :: problem_name = ''
  ! available problems, so far:
  !   - 'convection_diffusion__wave_package'
  !   - 'burgers__wave_package'
  !   - 'burgers__moving_front'

  real(RNP) :: t_start = 0.4  ! start time
  real(RNP) :: t_end   = 0.6  ! end time

  namelist/problem_prm/ problem_name, t_start, t_end

  ! declarations: discretization ...............................................

  integer, parameter :: max_n_level = 20 ! upper bound for number of levels

  real(RNP) :: dt_slab = 0.2  ! thickness of one time slab
  integer   :: n_level = 2    ! number of space-time levels

  namelist/discretization_prm/ dt_slab, n_level

  ! time integration -- still needs to be configured
  class(CL_TimeIntegrator_Options_1D), allocatable :: opt_pre
  class(CL_SDC_Options_1D)           , allocatable :: opt_sdc

  ! MLSDC
  type(CL_MLSDC_1D)          :: mlsdc
  type(CL_MLSDC_Options_1D)  :: mlsdc_opt

  integer :: n_space (max_n_level) = -1 ! number of elements in space
  integer :: p_space (max_n_level) = -1 ! polynomial degree of elements space
  integer :: n_time  (max_n_level) = -1 ! number of time steps in one slice
  integer :: p_time  (max_n_level) = -1 ! polynomial degree of time step

  namelist/discretization_prm/ n_space, p_space, n_time, p_time


  type(CL_MLSDC_Variable_1D) :: u_h, u_x

  real(RNP) :: t_0, t_1
  real(RNP) :: err_max
  logical   :: exists
  integer   :: io, stat
  integer   :: l

  ! initialization .............................................................

  ! greeting
  write(*,'(/,A)') 'MLSDC DG-SEM for 1D Conservation laws'

  ! identify case
  call get_command_argument(1, case_name, status=stat)
  if (stat /= 0 .or. len_trim(case_name) == 0) then
    call get_command_argument(0, case_name, status=stat)
  end if

  case_file = trim(case_name) // '.prm'
  inquire(file=case_file, exist=exists)
  if (exists) then
    open(newunit=io, file=case_file)
    read(io, nml = problem_prm)
    read(io, nml = discretization_prm)
    close(io)
  end if

  ! MLSDC options
  mlsdc_opt = CL_MLSDC_Options_1D(n_level)
  do l = 1, n_level
    mlsdc_opt % n_space (l) = n_space(l)
    mlsdc_opt % p_space (l) = p_space(l)
    mlsdc_opt % q_conv  (l) = (p_space(l) * 3 + 1) / 2
    mlsdc_opt % n_time  (l) = n_time(l)
    mlsdc_opt % p_time  (l) = p_time(l)
  end do
  write(*,*)
  write(*,'(A,99I5)') 'n_level = ', mlsdc_opt % n_level
  write(*,'(A,99I5)') 'n_space = ', mlsdc_opt % n_space
  write(*,'(A,99I5)') 'p_space = ', mlsdc_opt % p_space
  write(*,'(A,99I5)') 'q_conv  = ', mlsdc_opt % q_conv
  write(*,'(A,99I5)') 'n_time  = ', mlsdc_opt % n_time
  write(*,'(A,99I5)') 'p_time  = ', mlsdc_opt % p_time
  write(*,*)

  ! problem
  select case(problem_name)
  case('convection_diffusion__wave_package')
    write(*,'(A)') 'Initializing Convection-Diffusion Wave Package problem'
    allocate(CL_Problem_ConvectionDiffusion_WavePackage_1D :: cl_problem)
  case('burgers__wave_package')
    write(*,'(A)') 'Initializing Burgers Wave Package problem'
    allocate(CL_Problem_Burgers_WavePackage_1D :: cl_problem)
  case('burgers__moving_front')
    write(*,'(A)') 'Initializing Burgers Moving Front problem'
    allocate(CL_Problem_Burgers_MovingFront_1D :: cl_problem)
  case default
    call Error('Conservation_Law', 'Invalid problem name')
  end select
  call cl_problem % SetProblem() ! may not work with wave package
! use the following to configure the problem (later)
! call cl_problem % SetProblem(case_name)

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
  mlsdc = CL_MLSDC_1D(mlsdc_opt, opt_pre, opt_sdc, cl_problem)

  ! MLSDC variables
  u_h = CL_MLSDC_Variable_1D(mlsdc)
  u_x = CL_MLSDC_Variable_1D(mlsdc)

  ! fine-to-coarse projection test .............................................

  write(*,'(/,A)') 'fine-to-coarse projection test'

  t_0 = t_start
  t_1 = t_start + dt_slab

  ! set to exact solution
  do l = 1, n_level
    call GetExactSolution(mlsdc%level(l), t_0, t_1, u_h % level(l)%val)
    call GetExactSolution(mlsdc%level(l), t_0, t_1, u_x % level(l)%val)
  end do

  do l = n_level, 2, -1
    associate( u_hf => u_h % level(l  ) % val &
             , u_hc => u_h % level(l-1) % val &
             , u_xc => u_x % level(l-1) % val )

      call mlsdc % level(l) % Project_FC(u_hf, u_hc)

      err_max = maxval(abs(u_hc - u_xc))
      write(*,'(2X,2(A,I3),A,ES10.3)') 'level',l,' to',l-1,': err_max =',err_max

    end associate
  end do

  ! coarse-to-fine interpolation test ..........................................

  write(*,'(/,A)') 'coarse-to-fine interpolation test'

  call GetExactSolution(mlsdc%level(1), t_0, t_1, u_h % level(1)%val)

  do l = 1, n_level-1
    associate( u_hc => u_h % level(l  ) % val &
             , u_hf => u_h % level(l+1) % val &
             , u_xf => u_x % level(l+1) % val )

      call mlsdc % level(l) % Interpolate_CF(u_hc, u_hf, complete=.true.)

      err_max = maxval(abs(u_hf - u_xf))
      write(*,'(2X,2(A,I3),A,ES10.3)') 'level',l,' to',l+1,': err_max =',err_max

    end associate
  end do

!### CHECK
block
  real(RNP), allocatable :: t(:)
  real(RNP) :: dt
  character(len=80) plot_file
  integer :: io
  integer :: c, i, j, k, m, n
  c = 1
  do l = 1, n_level
    write(plot_file,'(9G0)') 'result_l', l, '.dat'
    open(newunit=io, file=plot_file)
    write(io,'(A)') '# x, u_h, u_x, err'

    ! spatial
    m = mlsdc % level(l) % p_time
    n = mlsdc % level(l) % n_time
    !m = 0
    !n = 1
    do j = 1, mlsdc % level(l) % n_space
    do i = 0, mlsdc % level(l) % p_space
      write(io,'(99(ES17.10,1X))')                &
        mlsdc % level(l) % cl_operator % x(i,j) , &
        u_h % level(l) % val(i,j,c,m,n)         , &
        u_x % level(l) % val(i,j,c,m,n)         , &
        u_h % level(l) % val(i,j,c,m,n) -         &
        u_x % level(l) % val(i,j,c,m,n)
    end do
    end do

!!     ! temporal
!!     allocate(t(0:mlsdc % level(l) % p_time))
!!     dt = (t_1 - t_0) / mlsdc % level(l) % n_time
!!     i = 0
!!     j = 1
!!     do n = 1, mlsdc % level(l) % n_time
!!       t(0:) = mlsdc % level(l) % cl_sdc % IntermediateTimes(t_0 + (n-1)*dt, dt)
!!       do m = 0, mlsdc % level(l) % p_time
!!         write(io,'(99(ES17.10,1X))')                &
!!           t(m)                                    , &
!!           u_h % level(l) % val(i,j,c,m,n)         , &
!!           u_x % level(l) % val(i,j,c,m,n)         , &
!!           u_h % level(l) % val(i,j,c,m,n) -         &
!!           u_x % level(l) % val(i,j,c,m,n)
!!       end do
!!     end do
!!     deallocate(t)

    close(io)
  end do
end block
!### CHECK END

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
