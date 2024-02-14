!> summary:  Testprogram for the MLSDC Components
!> author:   Joerg Stiller, Erik Pfister
!> date:     2024/01/22
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

program Test_Components_MLSDC
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
  use CL__SDC__Method__Euler__1D
  use CL__MLSDC__1D
  use CL__MLSDC__Level__1D
  use CL__MLSDC__Variable__1D
  use CL__MLSDC__V_Cycle__1D

  use Export_VTK_Spacetime_Data__1D

  implicit none

  ! declarations: control ......................................................

  character(len=80) :: case_name = 'test_components_ml'
  character(len=80) :: case_file

  ! declarations: problem ......................................................

  class(CL_Problem_1D), allocatable :: cl_problem

  character(len=80) :: problem_name = ''
  ! available problems, so far:
  !   - 'convection_diffusion__wave_package'
  !   - 'burgers__wave_package'
  !   - 'burgers__moving_front'

  real(RNP) :: t_start = 0.40  ! start time
  real(RNP) :: t_end   = 0.41  ! end time

  namelist/problem_prm/ problem_name, t_start, t_end

  ! declarations: discretization ...............................................

  integer, parameter :: max_n_level = 20 ! upper bound for number of levels

  real(RNP) :: dt_slab  = 0.01 ! thickness of one time slab
  integer   :: n_level  = 2    ! number of space-time levels
  integer   :: n_cycle  = 2    ! number of v-cycles
  integer   :: n_coarse = 2    ! number of coarse sweeps

  namelist/discretization_prm/ dt_slab, n_level, n_cycle, n_coarse

  integer :: p_space (max_n_level) = -1 ! polynomial degree of elements space
  integer :: n_space (max_n_level) = -1 ! number of elements in space
  integer :: p_time  (max_n_level) = -1 ! polynomial degree of time step
  integer :: n_time  (max_n_level) = -1 ! number of time steps in one slice
  integer :: n_sweep (max_n_level) = -1 ! number of SDC sweeps corrector

  namelist/discretization_prm/ p_space, n_space, p_time, n_time, n_sweep

  ! time integration -- still needs to be configured
  class(CL_TimeIntegrator_Options_1D), allocatable :: opt_pre
  class(CL_SDC_Options_1D)           , allocatable :: opt_sdc

  integer :: sdc_method = 2

  namelist/time_integration_prm/ sdc_method

  type(CL_SDC_Options_Euler_1D) :: opt_sdc_euler
  type(CL_SDC_Options_ISD1_1D)  :: opt_sdc_isd1

  namelist/time_integration_prm/ opt_sdc_euler, opt_sdc_isd1

  ! component testing
  integer :: projection_switch    ! switch for projection test
  integer :: interpolation_switch ! switch for interpolation test
  integer :: residual_switch      ! switch for collocation residual test
  integer :: restriction_switch   ! switch for restriction test
  integer :: predictor_switch     ! switch for predictor test
  integer :: corrector_switch     ! switch for corrector test
  integer :: vcycle_switch        ! switch for v cycle

  namelist/component_testing_prm/ projection_switch, interpolation_switch, &
                                  residual_switch, restriction_switch,     &
                                  predictor_switch, corrector_switch,      &
                                  vcycle_switch

  ! MLSDC
  type(CL_MLSDC_1D)          :: mlsdc
  type(CL_MLSDC_Options_1D)  :: mlsdc_opt

  type(CL_MLSDC_Variable_1D)     :: u_h, u_x
  real(RNP)                      :: t_0, t_1, dt
  real(RNP)                      :: err_max, r_max
  logical                        :: exists
  integer                        :: io, stat
  integer                        :: l, i, j, k, c
  character(len=80)              :: filename
  real(RNP), allocatable         :: t_ges(:,:)
  real(RNP), allocatable         :: s(:,:,:,:,:)
  real(RNP), allocatable         :: u_hf(:,:,:,:,:)
  real(RNP), allocatable         :: r_f(:,:,:,:,:)
  real(RNP), allocatable         :: r_c(:,:,:,:,:)
  real(RNP), allocatable         :: r_fc(:,:,:,:,:)
  real(RNP), allocatable         :: v(:,:,:,:,:)
  real(RNP), allocatable         :: g(:,:,:,:,:)
  character(len=20), allocatable :: sname(:)
  real(RNP)                      :: dt_f, int_f, int_c

  ! initialization .............................................................

  ! greeting
  write(*,'(/,A)') 'MLSDC Component Testprogram'

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
    read(io, nml = component_testing_prm)
    close(io)
  end if

  ! MLSDC options
  mlsdc_opt = CL_MLSDC_Options_1D(n_level)
  do l = 1, n_level
    mlsdc_opt % p_space (l) = p_space(l)
    mlsdc_opt % n_space (l) = n_space(l)
    mlsdc_opt % p_time  (l) = p_time(l)
    mlsdc_opt % n_time  (l) = n_time(l)
    mlsdc_opt % q_conv  (l) = (p_space(l) * 3 + 1) / 2
  end do
  write(*,*)
  write(*,'(A,99I5)') 'n_level = ', mlsdc_opt % n_level
  write(*,'(A,99I5)') 'p_space = ', mlsdc_opt % p_space
  write(*,'(A,99I5)') 'q_conv  = ', mlsdc_opt % q_conv
  write(*,'(A,99I5)') 'n_space = ', mlsdc_opt % n_space
  write(*,'(A,99I5)') 'p_time  = ', mlsdc_opt % p_time
  write(*,'(A,99I5)') 'n_time  = ', mlsdc_opt % n_time
  write(*,'(A,99I5)') 'n_sweep = ', n_sweep(1:n_level)
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

  ! convdiff
  call cl_problem % SetProblem(case_name)
  ! burgers
  !call cl_problem % SetProblem()

  ! predictor options (ISD1, so far)
  allocate(CL_TimeIntegrator_Options_ISD1_1D :: opt_pre)
  opt_pre % imex_mode        =   1
  opt_pre % diffusion_method =   4
  opt_pre % diffusion_i_max  = 100

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

  ! testing of the components ..................................................

  t_0 = t_start
  t_1 = t_start + dt_slab

  ! set to exact solution
  do l = 1, n_level
    call GetExactSolution(mlsdc%level(l), t_0, t_1, u_h % level(l)%val)
    call GetExactSolution(mlsdc%level(l), t_0, t_1, u_x % level(l)%val)
  end do

  ! fine-to-coarse projection test .............................................

  if (projection_switch > 0) then

    write(*,'(/,A)') 'fine-to-coarse projection test'

    do l = n_level, 2, -1
      associate( u_hf     => u_h       % level (l  ) % val    &
               , u_hc     => u_h       % level (l-1) % val    &
               , u_xc     => u_x       % level (l-1) % val    &
               , m_time_c => mlsdc     % level (l-1) % m_time &
               , n_time_c => mlsdc_opt % n_time(l-1)          )

        call mlsdc % level(l) % Project_FC(u_hf, u_hc)

        err_max = maxval(abs(u_hc - u_xc))
        write(*,'(2X,2(A,I3),A,ES10.3)') 'level',l,' to',l-1,': err_max =',err_max

        if (projection_switch > 1) then ! visualization

          allocate(t_ges(0:m_time_c, 1:n_time_c))
          call mlsdc % level(l-1) % GetTimeMesh(t_0, t_1, t_ges)
          call Pack_SpacetimeData(u=u_hc, s=s, uname=['u_hc'], sname=sname)
          call Pack_SpacetimeData(u=u_xc, s=s, uname=['u_xc'], sname=sname)
          call Pack_SpacetimeData(u=u_hc-u_xc, s=s, uname=['error'], sname=sname)

          write(filename, '(A,A,I0)') &
              trim(problem_name), '__fc_projection_to_level_',l-1

          call ExportVTK_SpacetimeData( x     = mlsdc%level(l-1)%cl_operator%x &
                                      , t     = t_ges                          &
                                      , s     = s                              &
                                      , sname = sname                          &
                                      , file  = trim(adjustl(filename))        )
          deallocate(t_ges, s, sname)

        end if

      end associate
    end do

  end if

  ! coarse-to-fine interpolation test ..........................................

  if (interpolation_switch > 0) then

    write(*,'(/,A)') 'coarse-to-fine interpolation test'

    call GetExactSolution(mlsdc%level(1), t_0, t_1, u_h % level(1)%val)

    do l = 1, n_level-1
      associate( u_hc     => u_h       % level (l  ) % val &
               , u_hf     => u_h       % level (l+1) % val &
               , u_xf     => u_x       % level (l+1) % val &
               , m_time_f => mlsdc % level(l+1) % m_time   &
               , n_time_f => mlsdc % level(l+1) % n_time   )

        call mlsdc % level(l) % Interpolate_CF(u_hc, u_hf, complete=.true.)

        err_max = maxval(abs(u_hf - u_xf))
        write(*,'(2X,2(A,I3),A,ES10.3)') 'level',l,' to',l+1,': err_max =',err_max

        if (interpolation_switch > 1) then ! visualization

          allocate(t_ges(0:m_time_f, 1:n_time_f))
          call mlsdc % level(l+1) % GetTimeMesh(t_0, t_1, t_ges)
          call Pack_SpacetimeData(u=u_hf, s=s, uname=['u_hf'], sname=sname)
          call Pack_SpacetimeData(u=u_xf, s=s, uname=['u_xf'], sname=sname)
          call Pack_SpacetimeData(u=u_hf-u_xf, s=s, uname=['error'], sname=sname)

          write(filename, '(A,A,I0)') &
              trim(problem_name),'__cf_interpol_to_level_',l+1

          call ExportVTK_SpacetimeData( x     = mlsdc%level(l+1)%cl_operator%x &
                                      , t     = t_ges                          &
                                      , s     = s                              &
                                      , sname = sname                          &
                                      , file  = trim(adjustl(filename))        )
          deallocate(t_ges, s, sname)

        end if

      end associate
    end do

  end if

   ! collocation residual test ..................................................


   if (residual_switch > 0) then

     write(*,'(/,A)') 'collocation residual test'

     do l = 1, n_level
       associate( u      => u_x       % level (l) % val    &
                , r      => u_h       % level (l) % val    &
                , m_time => mlsdc     % level (l) % m_time &
                , n_time => mlsdc_opt % n_time(l)          )

         call mlsdc % level(l) % GetResidual(dt_slab, t_0, u, r)

         r_max = maxval(abs(r))
         write(*,'(2X,A,I3,A,ES10.3)') 'level',l,': r_max =',r_max

         if (residual_switch > 1) then ! visualization

           allocate(t_ges(0:m_time, 1:n_time))
           call mlsdc % level(l) % GetTimeMesh(t_0, t_1, t_ges)
           call Pack_SpacetimeData(u=r, s=s, uname=['r'], sname=sname)
           write(filename,'(A,A,I0)') trim(problem_name),'__residual_level_',l

           call ExportVTK_SpacetimeData( x     = mlsdc%level(l)%cl_operator%x &
                                       , t     = t_ges                        &
                                       , s     = s                            &
                                       , sname = sname                        &
                                       , file  = trim(adjustl(filename))      )
           deallocate(t_ges, s, sname)

         end if

       end associate
     end do

   end if

   ! fine-to-coarse restriction test ..............................................

   if (restriction_switch > 0) then

     write(*,'(/,A)') 'fine-to-coarse restriction test'

     do l = n_level, 2, -1
       associate( u_xf      => u_x       % level  (l  ) % val              &
                , p_space_f => mlsdc_opt % p_space(l  )                    &
                , n_space_f => mlsdc_opt % n_space(l  )                    &
                , m_time_f  => mlsdc     % level  (l  ) % m_time           &
                , n_time_f  => mlsdc_opt % n_time (l  )                    &
                , m_time_c  => mlsdc     % level  (l-1) % m_time           &
                , n_time_c  => mlsdc_opt % n_time (l-1)                    &
                , r_hf      => u_h       % level  (l  ) % val              &
                , r_hc      => u_h       % level  (l-1) % val              &
                , Me_f      => mlsdc     % level  (l  ) % cl_operator % Me &
                , wt_f      => mlsdc     % level  (l  ) % cl_sdc      % w  &
                , nc        => mlsdc     % level  (l  ) % cl_problem  % nc )

         ! timestep fine
         dt_f = (t_1 - t_0) / n_time_f

         ! residual like variable
         do i = 1, n_space_f
           do c = 1, nc
             do j = 0, m_time_f
               do k = 1, n_time_f
                 r_hf(:,i,c,j,k) = Me_f * wt_f(j) * dt_f * u_xf(:,i,c,j,k)
               end do
             end do
           end do
         end do

         ! restrict
         call mlsdc % level(l-1) % Restrict_FC(r_hf, r_hc)

         ! integral
         int_f = sum(r_hf)
         int_c = sum(r_hc)

         err_max = abs(int_f - int_c)
         write(*,'(2X,2(A,I3),A,ES10.3)') 'level',l,' to',l-1,': err_max =',err_max

         if (restriction_switch > 1) then ! visualization

           allocate(t_ges(0:m_time_c, 1:n_time_c))
           call mlsdc%level(l-1)%GetTimeMesh(t_0, t_1, t_ges)
           call Pack_SpacetimeData(u=r_hc, s=s, uname=['r_h'], sname=sname)

           write(filename,'(A,A,I0)') &
               trim(problem_name),'__fc_restrict_to_level_',l-1

           call ExportVTK_SpacetimeData( x     = mlsdc%level(l-1)%cl_operator%x &
                                       , t     = t_ges                          &
                                       , s     = s                              &
                                       , sname = sname                          &
                                       , file  = trim(adjustl(filename))        )
           deallocate(t_ges, s, sname)

         end if

       end associate
     end do

   end if

   ! test of the predictor ......................................................

   if (predictor_switch > 0) then

     write(*,'(/,A)') 'test of the predictor'

     do l = n_level, 1, -1
       associate( u_x    => u_x       % level (l) % val    &
                , u_h    => u_h       % level (l) % val    &
                , m_time => mlsdc     % level (l) % m_time &
                , n_time => mlsdc_opt % n_time(l)          )

         ! set approximate solution to zero and initial condition
         u_h = 0.
         u_h(:,:,:,0,1) = u_x(:,:,:,0,1)
         ! predictor
         call mlsdc % level(l) % ApplyPredictor(dt_slab, t_0, u_h)

         err_max = maxval(abs(u_h - u_x))
         write(*,'(2X,A,I3,A,ES10.3)') 'error on level ',l,': err_max =',err_max

         if (predictor_switch > 1) then ! visualization

           allocate(t_ges(0:m_time, 1:n_time))
           call mlsdc%level(l)%GetTimeMesh(t_0, t_1, t_ges)
           call Pack_SpacetimeData(u=u_h, s=s, uname=['u_h'], sname=sname)
           call Pack_SpacetimeData(u=u_x, s=s, uname=['u_x'], sname=sname)
           call Pack_SpacetimeData(u=u_h-u_x, s=s, uname=['error'], sname=sname)

           write(filename,'(A,A,I0)') trim(problem_name),'__predictor_level_',l

           call ExportVTK_SpacetimeData( x     = mlsdc%level(l)%cl_operator%x &
                                       , t     = t_ges                        &
                                       , s     = s                            &
                                       , sname = sname                        &
                                       , file  = trim(adjustl(filename))      )
           deallocate(t_ges, s, sname)

         end if

       end associate
     end do

   end if

   ! test of the corrector ......................................................

   if (corrector_switch > 0) then

     write(*,'(/,A)') 'test of the corrector'

     do l = n_level, 1, -1
       associate( u_x     => u_x       % level  (l) % val              &
                , u_h     => u_h       % level  (l) % val              &
                , p_space => mlsdc_opt % p_space(l)                    &
                , n_space => mlsdc_opt % n_space(l)                    &
                , m_time  => mlsdc     % level  (l) % m_time           &
                , n_time  => mlsdc_opt % n_time (l)                    &
                , nc      => mlsdc     % level  (l) % cl_problem  % nc )

         ! set approximate solution to zero and initial condition
         u_h = 0.
         u_h(:,:,:,0,1) = u_x(:,:,:,0,1)

         ! predictor
         call mlsdc % level(l) % ApplyPredictor(dt_slab, t_0, u_h)

         ! set FAS defect correction
         allocate(g(0:p_space, 1:n_space, 1:nc, 0:m_time, 1:n_time))
         if (l == n_level) then
           g = 0.
         else
           allocate(v(0:p_space, 1:n_space, 1:nc, 0:m_time, 1:n_time))
           ! fine solution restriction
           call mlsdc % level(l+1) % Project_FC(u_hf, v)
           ! get residual
           allocate(r_c(0:p_space, 1:n_space, 1:nc, 0:m_time, 1:n_time))
           call mlsdc % level(l) % GetResidual(dt_slab, t_0, v, r_c)
           ! restrict fine residual
           allocate(r_fc(0:p_space, 1:n_space, 1:nc, 0:m_time, 1:n_time))
           call mlsdc % level(l) % Restrict_FC(r_f, r_fc)
           ! calculate G
           g = r_fc - r_c
           deallocate(u_hf, r_f, r_c, r_fc, v)
         end if

         ! corrector
         call mlsdc % level(l) % ApplyCorrector(dt_slab, t_0, g, u_h, n_sweep(l))
         deallocate(g)

         err_max = maxval(abs(u_h - u_x))
         write(*,'(2X,A,I3,A,ES10.3)') 'sdc error on level ',l,': err_max =',err_max

         ! save fine variables to calculate G
         allocate(u_hf(0:p_space, 1:n_space, 1:nc, 0:m_time,1:n_time))
         allocate(r_f(0:p_space, 1:n_space, 1:nc, 0:m_time,1:n_time))
         u_hf = u_h
         ! get residual
         call mlsdc % level(l) % GetResidual(dt_slab, t_0, u_hf, r_f)

         if (corrector_switch > 1) then ! visualization

           allocate(t_ges(0:m_time, 1:n_time))
           call mlsdc%level(l)%GetTimeMesh(t_0, t_1, t_ges)
           call Pack_SpacetimeData(u=u_h, s=s, uname=['u_h'], sname=sname)
           call Pack_SpacetimeData(u=u_x, s=s, uname=['u_x'], sname=sname)
           call Pack_SpacetimeData(u=u_h-u_x, s=s, uname=['error'], sname=sname)

           write(filename,'(A,A,I0)') trim(problem_name),'__corrector_level_',l

           call ExportVTK_SpacetimeData( x     = mlsdc%level(l)%cl_operator%x &
                                       , t     = t_ges                        &
                                       , s     = s                            &
                                       , sname = sname                        &
                                       , file  = trim(adjustl(filename))      )
           deallocate(t_ges, s, sname)

         end if

       end associate
     end do

   end if

  ! test of the v-cycle ........................................................

  if (vcycle_switch > 0) then

    ! set to exact solution
    do l = 1, n_level
      call GetExactSolution(mlsdc%level(l), t_0, t_1, u_h % level(l)%val)
      call GetExactSolution(mlsdc%level(l), t_0, t_1, u_x % level(l)%val)
    end do

    write(*,'(/,A)') 'test of the v-cycle'

    ! set u before entering the v_cycle
    do l = n_level, 1, -1
      associate( u_x     => u_x       % level  (l) % val &
               , u_h     => u_h       % level  (l) % val )

        ! set initial condition and apply predictor
        u_h = 0.
        u_h(:,:,:,0,1) = u_x(:,:,:,0,1)

        ! predictor
        call mlsdc % level(l) % ApplyPredictor(dt_slab, t_0, u_h)

        if (l == n_level) print*, "error before v cycle"
        err_max = maxval(abs(u_h - u_x))
        write(*,'(2X,A,I3,A,ES10.3)') 'error on level ',l,': err_max =',err_max

      end associate
    end do

    ! enter v cycle
    !call CL_MLSDC_V_Cycle_1D(mlsdc, t_0, dt_slab, 1, 1, n_coarse, n_cycle, u_h, u_x)
    call CL_MLSDC_V_Cycle_1D(mlsdc, t_0, dt_slab, 1, 1, n_coarse, n_cycle, u_h)

    do l = n_level, 1, -1
      if (l == n_level) print*, "error after ", n_cycle," cycles"
      err_max = maxval(abs(u_h%level(l)%val - u_x%level(l)%val))
      write(*,'(2X,A,I3,A,ES10.3)') 'v error on level ',l,': err_max =',err_max
    end do

  end if 

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

end program Test_Components_MLSDC
