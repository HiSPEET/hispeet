! PROVISORIUM
program Conservation_Law
  use Kind_Parameters
  use Constants
  use Array_Assignments

  use CL__Operator__1D
  use CL__Problem__1D
  use CL__Problem__Burgers__Moving_Front__1D

  use CL__Time_Integrator__1D
  use CL__Time_Integrator__Euler__1D
  use CL__Time_Integrator__ISD1__1D

  implicit none

  ! declarations: control ......................................................

  character(len=80) :: case_name = 'conservation_law'
  character(len=80) :: case_file

  ! declarations: problem ......................................................

  class(CL_Problem_1D), allocatable :: cl_problem

  ! declarations: space discretization .........................................

  type(CL_Operator_1D)        :: cl_operator
  type(CL_OperatorOptions_1D) :: cl_operator_opt

  namelist/space_discretization_prm/ cl_operator_opt

  ! declarations: time integration .............................................

  integer   :: time_method = 1

  real(RNP) :: t_end   =  0.1
  real(RNP) :: dt      =  0.001
  integer   :: nt_max  = -1

  namelist/time_integration_prm/ time_method, t_end, dt, nt_max

  class(CL_TimeIntegrator_1D), allocatable :: cl_tint
  type(CL_TimeIntegrator_Options_Euler_1D) :: cl_tint_euler_opt
  type(CL_TimeIntegrator_Options_ISD1_1D)  :: cl_tint_isd1_opt

  namelist/time_integration_prm/ cl_tint_euler_opt, &
                                 cl_tint_isd1_opt

  ! declarations: variables ....................................................

  real(RNP), dimension(:,:,:), allocatable :: u_0, u

  ! declarations: auxiliary ....................................................

  real(RNP) :: t, t_run, t_run_0
  logical   :: exists
  integer   :: io, nt
  integer   :: i, k

  ! initialization .............................................................

  ! greeting
  write(*,'(/,A)') 'DG-SEM for 1D Conservation laws'

  ! load parameters
  case_file = trim(case_name) // '.prm'
  inquire(file=case_file, exist=exists)
  if (exists) then
    open(newunit=io, file=case_file)
    read(io, nml = space_discretization_prm)
    read(io, nml = time_integration_prm)
    close(io)
  end if

  nt = nint(t_end / dt)
  if (nt_max > 0) then
    nt = min(nt, nt_max)
  end if

  ! problem (only one, so far)
  allocate(CL_Problem_Burgers_MovingFront_1D :: cl_problem)
  call cl_problem % SetProblem(case_name)

  ! operators
  cl_operator_opt % nc  = cl_problem % nc
  cl_operator_opt % xb1 = cl_problem % xb1
  cl_operator_opt % xb2 = cl_problem % xb2
  cl_operator = CL_Operator_1D(cl_operator_opt)

  ! time integration
  select case(time_method)
  case(1)
    cl_tint = CL_TimeIntegrator_Euler_1D(cl_tint_euler_opt)
  case(2)
    cl_tint = CL_TimeIntegrator_ISD1_1D(cl_tint_isd1_opt)
  end select

  ! show settings: TBD

  ! variables
  allocate(u(0:cl_operator%eop%po, cl_operator%ne, cl_problem%nc))
  allocate(u_0, mold = u)

  ! initial values
  call cl_problem % GetInitialValues(cl_operator, u)
  call SetArray(u_0, u)

  ! time integration ...........................................................

  t = 0
  k = 1
  write(*,*)

  call cpu_time(t_run_0)

  do i = 1, nt

    call cl_tint % TimeStep(cl_problem, cl_operator, dt, t, u_0, u)
    call SetArray(u_0, u)
    t = t + dt

    if (10 * t >= k * t_end) then
      write(*,'(2X,I3,"%")') 10*k
      k = k + 1
    end if

  end do

  call cpu_time(t_run)

  write(*,'(/,A,ES12.5)') 't_run =', t_run  - t_run_0

  ! save results
  open(newunit=io, file=trim(case_name)//'.dat')
  write(io,'(A)') '# x, u'
  do k = 1, cl_operator % ne
  do i = 0, cl_operator % eop % po
    write(io,'(99(ES17.10,1X))') cl_operator % x(i,k), u(i,k,:)
  end do
  end do
  close(io)

  !=============================================================================

end program Conservation_Law
