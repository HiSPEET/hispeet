program Conservation_Law_ML
  use Kind_Parameters
  use Constants
  use Execution_Control

  use CL__Operator__1D
  use CL__Problem__1D
  use CL__Problem__CNS__Acoustic_Wave__1D
  use CL__Problem__CNS__Shu_Osher__1D
  use CL__Problem__Convection_Diffusion__Wave_Package__1D
  use CL__Problem__Burgers__Wave_Package__1D
  use CL__Problem__Burgers__Moving_Front__1D

  use CL__Time_Integrator__1D
  use CL__Time_Integrator__Euler__1D
  use CL__Time_Integrator__ISD1__1D
  use CL__SDC__Method__1D
  use CL__SDC__Method__Euler__1D
  use CL__SDC__Method__ISD1__1D
  use CL__MLSDC__1D
  use CL__MLSDC__Level__1D
  use CL__MLSDC__Variable__1D
  use CL__MLSDC__Upward_Leg__1D
  use CL__MLSDC__V_Cycle__1D

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

  real(RNP) :: t_start = 0.00   ! start time
  real(RNP) :: t_end   = 1e-4   ! end time
  integer   :: nt_max  = -1     ! maximal time steps

  namelist/problem_prm/ problem_name, t_start, t_end, nt_max

  ! declarations: discretization ...............................................

  integer, parameter :: max_n_level = 20 ! upper bound for number of levels

  real(RNP) :: dt_slab  = 0.1  ! thickness of one time slab
  integer   :: n_level  = 2    ! number of space-time levels
  integer   :: n_cycle  = 2    ! number of v-cycles
  integer   :: n_coarse = 2    ! number of coarse sweeps

  namelist/discretization_prm/ dt_slab, n_level, n_cycle, n_coarse

  integer :: p_space (max_n_level) = -1 ! polynomial degree of elements space
  integer :: n_space (max_n_level) = -1 ! number of elements in space
  integer :: p_time  (max_n_level) = -1 ! polynomial degree of time step
  integer :: n_time  (max_n_level) = -1 ! number of time steps in one slice

  namelist/discretization_prm/ p_space, n_space, p_time, n_time

  ! MLSDC
  type(CL_MLSDC_1D)          :: mlsdc
  type(CL_MLSDC_Options_1D)  :: mlsdc_opt

  ! time integration -- still needs to be configured
  class(CL_TimeIntegrator_Options_1D), allocatable :: opt_pre
  class(CL_SDC_Options_1D)           , allocatable :: opt_sdc

  integer :: time_method = 2
  integer :: sdc_method  = 2

  namelist/time_integration_prm/ time_method, sdc_method

  type(CL_TimeIntegrator_Options_ISD1_1D) :: opt_pre_isd1
  type(CL_SDC_Options_Euler_1D)           :: opt_sdc_euler
  type(CL_SDC_Options_ISD1_1D)            :: opt_sdc_isd1

  namelist/time_integration_prm/ opt_pre_isd1, opt_sdc_euler, opt_sdc_isd1

  ! auxiliary variables
  type(CL_MLSDC_Variable_1D) :: u_h, u_x
  real(RNP)                  :: t_0, t_1, t
  !real(RNP), allocatable     :: err_2(:), err(:)
  real(RNP)                  :: err_2, err
  real(RNP)                  :: err_max, t_run, t_run_0
  logical                    :: exists
  integer                    :: io, stat
  integer                    :: l, nt, i, k, c

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
    read(io, nml = time_integration_prm)
    close(io)
  end if

  nt = nint((t_end-t_start) / dt_slab)
  if (nt_max >= 0) then
    nt = min(nt, nt_max)
  end if

  ! MLSDC options
  mlsdc_opt = CL_MLSDC_Options_1D(n_level)
  do l = 1, n_level
    mlsdc_opt % p_space (l) = p_space(l)
    mlsdc_opt % n_space (l) = n_space(l)
    mlsdc_opt % p_time  (l) = p_time(l)
    mlsdc_opt % n_time  (l) = n_time(l)
    if (index(problem_name, 'burgers') == 1) then
      mlsdc_opt % q_conv  (l) = (p_space(l) * 3 + 1) / 2
    elseif (index(problem_name, 'cns') == 1) then
      mlsdc_opt % q_conv  (l) = p_space(l) * 2
    else
      mlsdc_opt % q_conv  (l) = p_space(l)
    end if
  end do
  write(*,*)
  write(*,'(A,99I5)') 'n_level = ', mlsdc_opt % n_level
  write(*,'(A,99I5)') 'p_space = ', mlsdc_opt % p_space
  write(*,'(A,99I5)') 'q_conv  = ', mlsdc_opt % q_conv
  write(*,'(A,99I5)') 'n_space = ', mlsdc_opt % n_space
  write(*,'(A,99I5)') 'p_time  = ', mlsdc_opt % p_time
  write(*,'(A,99I5)') 'n_time  = ', mlsdc_opt % n_time
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
  case('cns__acoustic_wave')
    write(*,'(A)') 'Initializing CNS acoustic wave problem'
    allocate(CL_Problem_CNS_AcousticWave_1D :: cl_problem)
  case('cns__shu_osher')
    write(*,'(A)') 'Initializing CNS Shu-Osher problem'
    allocate(CL_Problem_CNS_ShuOsher_1D :: cl_problem)
  case default
    call Error('Conservation_Law', 'Invalid problem name')
  end select

  ! convdiff
  call cl_problem % SetProblem(case_name)
  ! burgers
  !call cl_problem % SetProblem()

  ! predictor options (ISD1, so far)
  select case(time_method)
  case(2)
    allocate(CL_TimeIntegrator_Options_ISD1_1D :: opt_pre)
    opt_pre = opt_pre_isd1
  end select 

  ! SDC options (ISD1 and Euler so far)
  select case(sdc_method)
  case(1)
    allocate(CL_SDC_Options_Euler_1D :: opt_sdc)
    opt_sdc = opt_sdc_euler
  case(2)
    allocate(CL_SDC_Options_ISD1_1D :: opt_sdc)
    opt_sdc = opt_sdc_isd1
  end select

  ! MLSDC data structure
  mlsdc = CL_MLSDC_1D(mlsdc_opt, opt_pre, opt_sdc, cl_problem)

  ! MLSDC variables
  u_h = CL_MLSDC_Variable_1D(mlsdc)
  u_x = CL_MLSDC_Variable_1D(mlsdc)

  ! initialiuation
  t_0 = t_start
  t_1 = t_start + dt_slab

  ! time integration ...........................................................

  write(*,*)
  write(*,'(A)') 'Time integration'
  write(*,'(A)') repeat('=',80)
  write(*,*)

  call cl_problem % GetInitialValues( mlsdc % level(1) % cl_operator     &
                                    , u_h   % level(1) % val(:,:,:,0,1)  )

  print*, u_h%level(1)%val(1,1,1,0,1)

  call cpu_time(t_run_0)

  t = t_0
  k = 1

  do i = 1, nt

    ! upward leg
    call CL_MLSDC_Upward_Leg_1D(mlsdc, t, dt_slab, u_h)

    ! check converged
    
    ! enter v cycle
    call CL_MLSDC_V_Cycle_1D(mlsdc, t, dt_slab, 1, 1, n_coarse, n_cycle, u_h, u_x)

  !  if (i<nt) then 
  !    ! set new initial value for the coarse grid
  !    u_h%level(1)%val(:,:,:,0 ,1 ) = u_h%level(1)%val(:,:,:                 &
  !                                                  ,mlsdc%level(1)%m_time &
  !                                                  ,n_time(1)             ) 
  !    u_h%level(1)%val(:,:,:,1:,2:) = 0
  !  end if 

    ! Update t
    t = t + dt_slab
    
    if (10 * (t-t_0) >= k * (t_end-t_start)) then
      write(*,'(2X,I3,"%")') 10*k
      k = k + 1
    end if

  end do 

  call cpu_time(t_run)

  ! evaluation and output of results ...........................................

  associate( u_h => u_h       % level  (n_level) % val              &
           , u_x => u_x       % level  (n_level) % val              &
           , ne  => mlsdc_opt % n_space(n_level)                    &
           , po  => mlsdc_opt % p_space(n_level)                    &
           , mt  => mlsdc     % level  (n_level) % m_time           &
           , nt  => mlsdc_opt % n_time (n_level)                    &
           , Me  => mlsdc     % level  (n_level) % cl_operator % Me &
           , nc  => mlsdc     % level  (n_level) % cl_problem  % nc )

    ! set exact solution for now
    do l = 1, n_level
      call GetExactSolution( mlsdc%level(l), t_start, t_end, u_x)
    end do

    open(newunit=io, file='conservation_law_ml.dat')
    do k = 1, ne
    do i = 0, po
      err   = u_h(i, k, 1, mt, nt) - u_x(i, k, 1, mt, nt)
      err_2 = err_2 + Me(i) * err**2
      !do c = 1,nc
       ! write(io,'(ES17.9E3,1X,ES17.9E3)') mlsdc%level(n_level) % cl_operator % x(i,k), u_h(i,k,nc,mt,nt)
       write(io,'(ES17.9E3,1X,ES17.9E3)') mlsdc%level(n_level) % cl_operator % x(i,k), u_h(i,k,1,mt,nt) 
     !end do
    end do
    end do
    err_max = maxval(abs(u_h(:,:,:,mt,nt) - u_x(:,:,:,mt,nt)))
    close(io)
 
    write(*,'(2X,A,99(ES12.5,1X))') 't_run  =', t_run  - t_run_0
    write(*,'(2X,A,I3,A,ES10.3)') 'Error on level ',n_level,': err_2 =',sqrt(err_2)
    write(*,'(2X,A,I3,A,ES10.3)') 'Error on level ',n_level,': err_max =',err_max

  end associate

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
    integer   :: nt, mt
    integer   :: i, j

    mt = ubound(u,4)
    nt = ubound(u,5)
    dt = (t_1 - t_0) / nt

    allocate(t(0:mt))

    associate( cl_problem  => level % cl_problem  &
             , cl_operator => level % cl_operator &
             , cl_sdc      => level % cl_sdc      )

      do j = 1, nt
        t(0:) = cl_sdc % SubintervalPoints(t_0 + (j-1)*dt, dt)
        do i = 0, mt
          call cl_problem % GetExactSolution(cl_operator, t(i), u(:,:,:,i,j))
        end do
      end do

    end associate

  end subroutine GetExactSolution

  !=============================================================================

end program Conservation_Law_ML
