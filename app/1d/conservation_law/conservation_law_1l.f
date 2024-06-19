! PROVISORIUM
program Conservation_Law_1L
  use Kind_Parameters
  use Constants
  use Array_Assignments
  use Execution_Control
  use Logging_Levels

  use CL__Operator__1D
  use CL__Problem__1D
  use CL__Problem__Convection_Diffusion__Wave_Package__1D
  use CL__Problem__Burgers__Sine_Wave__1D
  use CL__Problem__Burgers__Moving_Front__1D
  use CL__Problem__Burgers__Wave_Package__1D
  use CL__Problem__CNS__Acoustic_Wave__1D
  use CL__Problem__CNS__Contact_Layer__1D
  use CL__Problem__CNS__Shock_Tube__1D
  use CL__Problem__CNS__Shu_Osher__1D

  use CL__Time_Integrator__1D
  use CL__Time_Integrator__Euler__1D
  use CL__Time_Integrator__ISD1__1D
  use CL__Time_Integrator__ISD2__1D
  use CL__Time_Integrator__RK__1D
  use CL__Time_Integrator__TVD_RK3__1D

  use CL__SDC__Method__1D
  use CL__SDC__Method__Euler__1D
  use CL__SDC__Method__ISD1__1D

  use, intrinsic :: IEEE_Arithmetic

  implicit none

  ! declarations: control ......................................................

  character(len=80) :: case_name = 'conservation_law_1l'
  character(len=80) :: case_file

  ! namelist input of log levels imported from module `Logging_Levels`
  namelist/control_prm/ log_level
  namelist/control_prm/ log_level_inner_iteration
  namelist/control_prm/ log_level_outer_iteration

  ! declarations: problem ......................................................

  class(CL_Problem_1D), allocatable :: cl_problem

  character(len=80) :: problem_name = ''
  ! available problems, so far:
  !   - 'convection_diffusion__wave_package'
  !   - 'burgers__wave_package'
  !   - 'burgers__moving_front'
  !   - 'cns__acoustic_wave'
  !   - 'cns__contact_layer'
  !   - 'cns__shock_tube'
  !   - 'cns__shu_osher'

  namelist/problem_prm/ problem_name

  ! declarations: space discretization .........................................

  type(CL_Operator_1D)         :: cl_operator
  type(CL_Operator_Options_1D) :: cl_operator_opt

  namelist/space_discretization_prm/ cl_operator_opt

  ! declarations: time integration .............................................

  real(RNP) :: t_end     =  0.1    ! problem time to reach
  real(RNP) :: dt        =  0.001  ! time step width
  integer   :: nt_max    = -1      ! max number of time steps
  integer   :: nt                  ! number of time steps to execute
  logical   :: adjust_dt = .false. ! adjust time step to dt = t_end / nt

  namelist/time_integration_prm/ t_end, dt, nt_max, adjust_dt

  integer   :: time_method = 1  ! standalone integrator or predictor
  integer   :: sdc_method  = 0  ! SDC method

  namelist/time_integration_prm/ time_method, sdc_method

  class(CL_TimeIntegrator_1D), allocatable   :: cl_tint
  type(CL_TimeIntegrator_Options_Euler_1D)   :: cl_tint_euler_opt
  type(CL_TimeIntegrator_Options_ISD1_1D)    :: cl_tint_isd1_opt
  type(CL_TimeIntegrator_Options_ISD2_1D)    :: cl_tint_isd2_opt
  type(CL_TimeIntegrator_Options_RK_1D)      :: cl_tint_rk_opt
  type(CL_TimeIntegrator_Options_TVD_RK3_1D) :: cl_tint_tvd_rk3_opt

  namelist/time_integration_prm/ cl_tint_euler_opt,  &
                                 cl_tint_isd1_opt,   &
                                 cl_tint_isd2_opt,   &
                                 cl_tint_rk_opt,     &
                                 cl_tint_tvd_rk3_opt

  class(CL_SDC_Method_1D), allocatable :: cl_sdc
  type(CL_SDC_Options_Euler_1D) :: cl_sdc_euler_opt
  type(CL_SDC_Options_ISD1_1D)  :: cl_sdc_isd1_opt

  namelist/time_integration_prm/ cl_sdc_euler_opt, cl_sdc_isd1_opt

  ! declarations: variables ....................................................

  real(RNP), dimension(:,:,:), allocatable :: u_0, u

  ! declarations: auxiliary ....................................................

  real(RNP), allocatable :: err(:), err_2(:)
  real(RNP) :: t, t_run, t_run_0
  real(RNP) :: tau_conv, tau_diff
  logical   :: exists
  integer   :: io, stat
  integer   :: i, k

  ! initialization .............................................................

  ! greeting
  write(*,'(/,A)') 'DG-SEM for 1D Conservation laws'

  ! identify case
  call get_command_argument(1, case_name, status=stat)
  if (stat /= 0 .or. len_trim(case_name) == 0) then
    call get_command_argument(0, case_name, status=stat)
  end if

  case_file = trim(case_name) // '.prm'
  inquire(file=case_file, exist=exists)
  if (exists) then
    open(newunit=io, file=case_file)
    read(io, nml = control_prm)
    read(io, nml = problem_prm)
    read(io, nml = space_discretization_prm)
    read(io, nml = time_integration_prm)
    close(io)
  end if

  nt = nint(t_end / dt)
  if (nt_max >= 0) then
    nt = min(nt, nt_max)
  end if
  if (adjust_dt .and. nt > 0) then
    dt = t_end / nt
  end if

  ! problem
  select case(problem_name)
  case('convection_diffusion__wave_package')
    write(*,'(A)') 'Initializing Convection-Diffusion Wave Package problem'
    allocate(CL_Problem_ConvectionDiffusion_WavePackage_1D :: cl_problem)
  case('burgers__sine_wave')
    write(*,'(A)') 'Initializing Burgers Sine Wave problem'
    allocate(CL_Problem_Burgers_SineWave_1D :: cl_problem)
  case('burgers__moving_front')
    write(*,'(A)') 'Initializing Burgers Moving Front problem'
    allocate(CL_Problem_Burgers_MovingFront_1D :: cl_problem)
  case('burgers__wave_package')
    write(*,'(A)') 'Initializing Burgers Wave Package problem'
    allocate(CL_Problem_Burgers_WavePackage_1D :: cl_problem)
  case('cns__acoustic_wave')
    write(*,'(A)') 'Initializing CNS acoustic wave problem'
    allocate(CL_Problem_CNS_AcousticWave_1D :: cl_problem)
  case('cns__contact_layer')
    write(*,'(A)') 'Initializing CNS contact layer problem'
    allocate(CL_Problem_CNS_ContactLayer_1D :: cl_problem)
  case('cns__shock_tube')
    write(*,'(A)') 'Initializing CNS shock tube problem'
    allocate(CL_Problem_CNS_ShockTube_1D :: cl_problem)
  case('cns__shu_osher')
    write(*,'(A)') 'Initializing CNS Shu-Osher problem'
    allocate(CL_Problem_CNS_ShuOsher_1D :: cl_problem)
  case default
    call Error('Conservation_Law', 'Invalid problem name')
  end select
  call cl_problem % SetProblem(case_name)

  ! operators
  cl_operator_opt % xb1 = cl_problem % xb1
  cl_operator_opt % xb2 = cl_problem % xb2
  cl_operator = CL_Operator_1D(cl_operator_opt)

  ! time integration
  select case(sdc_method)

  case(1)

    ! SDC based on Euler
    select case(time_method)
    case(1)
      cl_sdc = CL_SDC_Method_Euler_1D(cl_tint_euler_opt, cl_sdc_euler_opt)
    case(2)
      cl_sdc = CL_SDC_Method_Euler_1D(cl_tint_isd1_opt, cl_sdc_euler_opt)
    case(3)
      cl_sdc = CL_SDC_Method_Euler_1D(cl_tint_isd2_opt, cl_sdc_euler_opt)
    case(4)
      cl_sdc = CL_SDC_Method_Euler_1D(cl_tint_rk_opt, cl_sdc_euler_opt)
    end select

  case(2)

    ! SDC based on ISD1
    select case(time_method)
    case(1)
      cl_sdc = CL_SDC_Method_ISD1_1D(cl_tint_euler_opt, cl_sdc_isd1_opt)
    case(2)
      cl_sdc = CL_SDC_Method_ISD1_1D(cl_tint_isd1_opt, cl_sdc_isd1_opt)
    case(3)
      cl_sdc = CL_SDC_Method_ISD1_1D(cl_tint_isd2_opt, cl_sdc_isd1_opt)
    case(4)
      cl_sdc = CL_SDC_Method_ISD1_1D(cl_tint_rk_opt, cl_sdc_isd1_opt)
    end select

  case default

    ! standalone integrator
    select case(time_method)
    case(1)
      cl_tint = CL_TimeIntegrator_Euler_1D(cl_tint_euler_opt)
    case(2)
      cl_tint = CL_TimeIntegrator_ISD1_1D(cl_tint_isd1_opt)
    case(3)
      cl_tint = CL_TimeIntegrator_ISD2_1D(cl_tint_isd2_opt)
    case(4)
      cl_tint = CL_TimeIntegrator_RK_1D(cl_tint_rk_opt)
    case(5)
      cl_tint = CL_TimeIntegrator_TVD_RK3_1D(cl_tint_tvd_rk3_opt)
    end select
  end select

  ! show settings
  if (sdc_method > 0) then
    call cl_sdc % Show()
  else
    call cl_tint % Show()
  end if

  ! variables
  allocate(u(0:cl_operator%eop%po, cl_operator%ne, cl_problem%nc))
  allocate(u_0, mold = u)

  ! initial values
  call cl_problem % GetInitialValues(cl_operator, u)
  if (cl_problem % limiting_initial) then
    if (cl_problem % limiting_method == 1) then
      call cl_problem % MomentLimiter(cl_operator, u)
    end if
  end if
  call SetArray(u_0, u)

  ! time integration ...........................................................

  write(*,*)
  write(*,'(A)') 'Time integration'
  write(*,'(A,/)') repeat('=',80)
  write(*,'(2X,A,ES12.5)') 'dt  =' , dt
  write(*,'(2X,A,I0)')     'nt  = ', nt
  write(*,*)

  ! CFL and diffusion numbers
  call cl_problem % GetTimeScales(cl_operator, u, tau_conv, tau_diff)
  if (tau_conv > 0) write(*,'(2X,A,ES12.5)') 'c_conv =', dt/tau_conv
  if (tau_diff > 0) write(*,'(2X,A,ES12.5)') 'c_diff =', dt/tau_diff
  write(*,*)

  t = 0
  k = 1

  call cpu_time(t_run_0)

  do i = 1, nt

    if (sdc_method > 0) then
      call cl_sdc % TimeStep(cl_problem, cl_operator, dt, t, u_0, u)
    else
      call cl_tint % TimeStep(cl_problem, cl_operator, dt, t, u_0, u)
    end if

    if (any(ieee_is_nan(u))) then
      call Error('Conservation_Law_1L', 'detected NaN')
    end if

    call SetArray(u_0, u)
    t = t + dt

    if (10 * t >= k * t_end) then
      write(*,'(2X,I3,"%")') 10*k
      if (k == 5) then
        call cl_problem % GetTimeScales(cl_operator, u, tau_conv, tau_diff)
        if (tau_conv > 0) write(*,'(2X,A,ES12.5)') 'c_conv =', dt/tau_conv
        if (tau_diff > 0) write(*,'(2X,A,ES12.5)') 'c_diff =', dt/tau_diff
        write(*,*)
      end if
      k = k + 1
    end if

  end do

  call cpu_time(t_run)

  call cl_problem % GetTimeScales(cl_operator, u, tau_conv, tau_diff)
  if (tau_conv > 0) write(*,'(2X,A,ES12.5)') 'c_conv =', dt/tau_conv
  if (tau_diff > 0) write(*,'(2X,A,ES12.5)') 'c_diff =', dt/tau_diff
  write(*,*)

  ! evaluation and output of results ...........................................

  associate( nc => cl_problem  % nc        &
           , ne => cl_operator % ne        &
           , po => cl_operator % eop % po  &
           , Me => cl_operator % Me        )

    if (cl_problem % HasExactSolution()) then
      call cl_problem % GetExactSolution(cl_operator, t, u_0)
    else
      u_0 = u + 1
    end if

    allocate(err  (nc), source = ZERO)
    allocate(err_2(nc), source = ZERO)

    open(newunit=io, file=trim(case_name)//'.dat')

    write(io,'(A)') '# x, u, u_ex, err'
    do k = 1, ne
    do i = 0, po
      err   = u(i,k,:) - u_0(i,k,:)
      err_2 = err_2 + Me(i) * err**2
      write(io,'(99(ES17.9E3,1X))') cl_operator % x(i,k), u(i,k,:), u_0(i,k,:), err
    end do
    end do
    close(io)

    t_run = t_run  - t_run_0
    if (cl_problem % HasExactSolution()) then
      err_2 = sqrt(err_2)
    else
      err_2 = -1
    end if

    write(*,*)
    write(*,'(2X,A,99(ES12.5,1X))') 't_end  =', t
    write(*,'(2X,A,99(ES12.5,1X))') 't_run  =', t_run
    write(*,'(2X,A,99(ES12.5,1X))') 'err_2  =', err_2

    open(newunit=io, file=trim(case_name)//'.run')
    write(io,'(A)',advance='NO') '# cfl, dt, t_end, t_run'
    do k = 1, nc
      write(io,'(A,I0)',advance='NO') ', err_',k
    end do
    write(io,'(A)',advance='YES')
    write(io,'(99(ES17.9E3,1X))') dt/tau_conv, dt, t_end, t_run, err_2
    close(io)

  end associate

  !=============================================================================

end program Conservation_Law_1L
