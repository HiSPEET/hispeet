! PROVISORIUM
program Conservation_Law
  use Kind_Parameters
  use Constants

  use DG__Element_Operators__1D
  use DG__Elliptic_Operator__1D

  use CL__Problem__Scalar__1D
  use CL__Problem__Scalar__Burgers__Breaking_Wave__1D

  use CL__Time_Integrator__1D
  use CL__Time_Integrator__TVDRK__1D
  use CL__Time_Integrator__Euler__1D
! use CL__Time_Integrator__ISD__1D
  use CL__Time_Integrator__RK__1D

  use CL__SDC__Method__1D
  use CL__SDC__Method__Euler__1D
! use CL__SDC__Method__ISD__1D
  use CL__SDC__Method__RK__1D
  use CL__SDC__Method__TVDRK__1D

  implicit none
  class(CL_Problem_Scalar_1D), allocatable :: problem
  type(DG_ElementOptions_1D) :: dg_opt
  type(DG_SchwarzOptions_1D) :: schwarz_opt
  real(RNP), allocatable     :: u(:,:,:)     ! discrete solution
  real(RNP), allocatable     :: M_inv(:,:,:) ! transformation matrix

  real(RNP) :: start, finish !for time measurement

  integer   :: po         ! polynomial degree
  integer   :: ne         ! number of elements
  real(RNP) :: penalty    ! disc. gal. penalty
  real(RNP) :: cfl
  real(RNP) :: t_end
  namelist /discretization_prm/ po, ne, penalty, cfl, t_end

  logical   :: exists
  integer   :: prm

  real(RNP) :: dt, t
  integer   :: e, i, k, ou

  integer   :: time_method ! standalone or predictor method
  integer   :: sdc_method  !                     SDC method
  namelist /input/ time_method, sdc_method


  ! standalone time-integrator or predictor ...................................

  class(CL_TimeIntegrator_1D), allocatable :: tint
  type(CL_TimeIntegrator_Options_Euler_1D) :: tint_opt_eu    ! options for Euler
  !type(CL_TimeIntegrator_Options_ISD_1D)   :: tint_opt_sd    ! options for ISD
  type(CL_TimeIntegrator_Options_RK_1D)    :: tint_opt_rk    ! options for RK
  type(CL_TimeIntegrator_Options_TVDRK_1D) :: tint_opt_tvdrk ! options for tvdrk
  namelist /input/ tint_opt_eu, tint_opt_rk, tint_opt_tvdrk!, tint_opt_sd


  ! SDC .......................................................................

  class(CL_SDC_Method_1D), allocatable :: sdc
  type(CL_SDC_Options_Euler_1D)        :: sdc_opt_eu ! options for Euler-based SDC
  !type(CL_SDC_Options_ISD_1D)          :: sdc_opt_sd ! options for ISD-based SDC
  type(CL_SDC_Options_RK_1D)           :: sdc_opt_rk ! options for RK-based SDC
  type(CL_SDC_Options_TVDRK_1D)        :: sdc_opt_tvdrk ! options for TVDRK-based SDC
  namelist /input/ sdc_opt_eu, sdc_opt_rk, sdc_opt_tvdrk!, sdc_opt_sd

  character(len=80) :: input_file = 'time_integration'
  integer           :: io

  ! read options and parameters ................................................

  input_file = trim(input_file) // '.prm'

  open(newunit=io, file=input_file)
  read(io, nml=input)
  close(io)

  ! initialize time-integration method .........................................

  select case(sdc_method)
  case(0) ! standalone time integrator
    select case(time_method)
    case(1)
      tint = CL_TimeIntegrator_Euler_1D(tint_opt_eu)
    case(2)
      print *, ' '
      print *, ' '
      print *, 'Actually there is no usable implementation of ISD at the moment!'
      stop
!      tint = CL_TimeIntegrator_ISD_1D(tint_opt_sd)
    case(3)
      tint = CL_TimeIntegrator_RK_1D(tint_opt_rk)
    case(4)
      tint = CL_TimeIntegrator_TVDRK_1D(tint_opt_tvdrk)
    case default
      tint = CL_TimeIntegrator_RK_1D(tint_opt_rk)
    end select

  case(1) ! SDC based on Euler
    select case(time_method)
    case(1)
      sdc = CL_SDC_Method_Euler_1D(tint_opt_eu, sdc_opt_eu)
    case(2)
        print *, ' '
        print *, ' '
        print *, 'Actually there is no usable implementation of ISD at the moment!'
        stop
  !    sdc = CL_SDC_Method_Euler_1D(tint_opt_sd, sdc_opt_eu)
    case(3)
      sdc = CL_SDC_Method_Euler_1D(tint_opt_rk, sdc_opt_eu)
    case(4)
      sdc = CL_SDC_Method_Euler_1D(tint_opt_tvdrk, sdc_opt_eu)
    case default
      sdc = CL_SDC_Method_Euler_1D(tint_opt_rk, sdc_opt_eu)
    end select

  case(2) ! SDC based on ISD
    print *, ' '
    print *, ' '
    print *, 'Actually there is no usable implementation of SDC-ISD at the moment!'
    stop
    select case(time_method)
      case(1)
    !    sdc = CL_SDC_Method_ISD_1D(tint_opt_eu, sdc_opt_sd)
      case(2)
    !    sdc = CL_SDC_Method_ISD_1D(tint_opt_sd, sdc_opt_sd)
      case(3)
    !    sdc = CL_SDC_Method_ISD_1D(tint_opt_rk, sdc_opt_sd)
    case(4)
    !    sdc = CL_SDC_Method_ISD_1D(tint_opt_tvdrk, sdc_opt_sd)
    case default
    !    sdc = CL_SDC_Method_ISD_1D(tint_opt_rk, sdc_opt_sd)
    end select

  case(3) ! SDC based on RK
    select case(time_method)
    case(1)
      sdc = CL_SDC_Method_RK_1D(tint_opt_eu, sdc_opt_rk)
    case(2)
        print *, ' '
        print *, ' '
        print *, 'Actually there is no usable implementation of ISD at the moment!'
        stop
    case(3)
      sdc = CL_SDC_Method_RK_1D(tint_opt_rk, sdc_opt_rk)
    case(4)
      sdc = CL_SDC_Method_RK_1D(tint_opt_tvdrk, sdc_opt_rk)
    case default
      sdc = CL_SDC_Method_RK_1D(tint_opt_rk, sdc_opt_rk)
  end select

  case(4) ! SDC based on TVDRK
    select case(time_method)
    case(1)
      sdc = CL_SDC_Method_TVDRK_1D(tint_opt_eu, sdc_opt_tvdrk)
    case(2)
        print *, ' '
        print *, ' '
        print *, 'Actually there is no usable implementation of ISD at the moment!'
        stop
    !  sdc = CL_SDC_Method_TVDRK_1D(tint_opt_sd, sdc_opt_tvdrk)
    case(3)
      sdc = CL_SDC_Method_TVDRK_1D(tint_opt_rk, sdc_opt_tvdrk)
    case(4)
      sdc = CL_SDC_Method_TVDRK_1D(tint_opt_tvdrk, sdc_opt_tvdrk)
    case default
      sdc = CL_SDC_Method_TVDRK_1D(tint_opt_rk, sdc_opt_tvdrk)
    end select

  case default ! SDC based on Runge-Kutta
    select case(time_method)
    case(1)
      sdc = CL_SDC_Method_RK_1D(tint_opt_eu, sdc_opt_rk)
    case(2)
      print *, ' '
      print *, ' '
      print *, 'Actually there is no usable implementation of ISD at the moment!'
      stop
    !  sdc = CL_SDC_Method_RK_1D(tint_opt_sd, sdc_opt_rk)
    case(3)
      sdc = CL_SDC_Method_RK_1D(tint_opt_rk, sdc_opt_rk)
    case(4)
      sdc = CL_SDC_Method_RK_1D(tint_opt_tvdrk, sdc_opt_rk)
    case default
      sdc = CL_SDC_Method_RK_1D(tint_opt_rk, sdc_opt_rk)
    end select

  end select

  ! show settings ..............................................................
  select case(sdc_method)
  case(0)
    call tint % Show()
  case default
    call sdc  % Show()
  end select

  ! read options and parameters ................................................
  inquire(file='conservation_law.prm', exist=exists)
  if (exists) then
    open(newunit=prm, file='conservation_law.prm', action='READ')
    read(prm, nml=discretization_prm)
    close(prm)
  end if

  dg_opt = DG_ElementOptions_1D( po         =  po        &
                               , penalty    =  penalty   &
                               , hybrid     = .true.     &
                               , svv        = .false.    )

  allocate(CL_Problem_Scalar_Burgers_BreakingWave_1D :: problem)
  call problem % SetProblem('cl__problem__scalar__burgers__breaking_wave__1d')
  call problem % SetSpaceDiscretization(dg_opt, schwarz_opt, ne)
  allocate(u(0:po, ne, problem%nc))
  u(0:,:,:) = problem % InitialValues()
  dt = cfl * problem%dx / po**2

  allocate(M_inv, mold = u)
  do e = 1, problem % ne
  do k = 1, problem % nc
    M_inv(:,e,k) = (2 / problem % dx) / problem % eop % w
  end do
  end do

  k = 1
  t = 0

  call cpu_time(start) ! time measurement

  do
    select case(sdc_method)
    case(0)      ! standalone time-integrator
      call tint % TimeStep(problem, t, dt, u)
      t  = t + dt
    case default ! SDC method
      call sdc % TimeStep(problem, t, dt, u, M_inv)
      t  = t + dt
    end select

    if (10 * t >= k * t_end) then
      print '(2X,I3,"%")', 10*k
      k = k + 1
    end if
    if (t >= t_end) exit
  end do

  if(isnan(maxval(abs(u)))) then   ! check if component is nan
    print*, ""
    print*, ""
    print *, "         Solution exploding!"
    stop
  end if

  call cpu_time(finish) ! just for time-measurement
  print*, ""
  print*, ""
  print '("         Integration time = ",f20.3," seconds.")',finish-start
  print*, ""
  print*, ""

  ! save results
  open(newunit=ou, file='conservation_law.dat')
  write(ou,'(A)') '# x, u'
  do k = 1, ne
  do i = 0, po
    write(ou,'(10(ES17.10,1X))') problem % x(i,k), u(i,k,:)
  end do
  end do

  call ExactSolution( problem %x, t_end, ou ) 

  close(ou)

end program Conservation_Law
