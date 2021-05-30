!> summary:  Validation of incompressible single-phase flow
!> author:   Joerg Stiller
!> date:     2017/08/11
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> If present, the first argument of the invoking command will be interpreted
!> as the base name of the control file. If omitted, the program looks for
!> `isp_flow__sdc_test.prm`.
!===============================================================================

program ISP_Flow__SDC_Test

  use, intrinsic :: IEEE_Arithmetic, only: ieee_is_nan

  use Kind_Parameters, only: RNP, RDP, IXL
  use Constants,       only: ONE, ZERO, HALF
  use Array_Assignments
  use Array_Reductions
  use Standard_Operators__1D
  use Execution_Control
  use TPO__Grad_3d_R
  use Export_Volume_Data_To_VTK
  use OpenMP_Binding
  use XMPI

  use ISP_Flow_Problem__Test_Suite

  use CART__Mesh_Partition
  use CART__Generate_Structured_Mesh
  use CART__DG_Weak_Divergence
  use CART__DG_Weak_Gradient
  use CART__Elliptic_PMG
  use CART__ISP_Flow__Operators
  use CART__ISP_Flow__Time_Derivative
  use CART__ISP_Flow__Pressure

  ! time integration
  use CART__ISP_Flow__Euler                        ! old standalone Euler
  use CART__ISP_Flow__SDC                          ! current SDC
  use CART__ISP_Flow__Time_Integrator              ! TI base type
  use CART__ISP_Flow__Time_Integrator__Euler       ! standalone Euler
  use CART__ISP_Flow__Time_Integrator__BDF2        ! standalone BDF2
  use CART__ISP_Flow__Time_Integrator__TR          ! standalone trapezodial rule
  use CART__ISP_Flow__Time_Integrator__Runge_Kutta ! standalone Runge-Kutta
  use CART__ISP_Flow__SDC_Method                   ! SDC base type
  use CART__ISP_Flow__SDC_Method__Euler            ! SDC Euler method

  implicit none

  !-----------------------------------------------------------------------------
  ! Declarations

  ! MPI and OpenMP .............................................................

  type(MPI_Comm) :: comm         ! MPI communicator
  integer        :: rank         ! local MPI rank
  integer        :: n_proc       ! number of MPI processes
  integer        :: n_thread     ! number of OpenMP threads

  ! control ....................................................................

  character(len=80) :: control_file

  character(len=80) :: flow_type  = 'Vortex_TG'
  character(len=80) :: flow_case  = ''

  logical :: restart    = .false. ! continue computation from snapshot
  logical :: compute_p  = .false. ! recompute p at the end of time step
  logical :: write_stat = .false. ! evaluate+print solution metrics each step
  logical :: write_vtk  = .false. ! export results to VTK file
  logical :: eval_rms   = .true.  ! evaluate rms errors and divergence
  logical :: eval_max   = .false. ! evaluate max errors and divergence
  logical :: eval_eps   = .false. ! evaluate dissipation, ε = ∫ν∇v:∇v dΩ
  logical :: average    = .false. ! compute temporal averages
  integer :: i_avg      =  10     ! sampling rate for temporal averaging
  type(FlowOpControl) :: flow_op_control

  namelist /control/ flow_type, flow_case,         &
                     restart, compute_p,           &
                     write_stat, write_vtk,        &
                     eval_rms, eval_max, eval_eps, &
                     average, i_avg,               &
                     flow_op_control

  ! discretization parameters ..................................................

  integer   :: np(3)     =  1        ! number of partitions in directions 1:3
  integer   :: ep(3)     =  2        ! elements per partition and direction
  integer   :: po_u      =  2        ! order of u
  integer   :: po_p      = -1        ! order of pressure
  integer   :: po_q      = -1        ! order for quadrature of nonlinear terms
  real(RNP) :: penalty   =  2        ! penalty parameter (> 1)
  logical   :: adjust_dx = .false.   ! adjust mesh spacing: Δx₃ = max(Δx₁,Δx₂)

  namelist /discretization/ np, ep, po_u, po_p, po_q, penalty, adjust_dx

  ! elliptic solver parameters .................................................

  type(PMG_Options3D) :: pmg_u_opt
  type(PMG_Options3D) :: pmg_p_opt

  namelist /elliptic_solver/ pmg_u_opt, pmg_p_opt

  ! time integration ...........................................................

  real(RNP) :: t                     ! problem time
  real(RNP) :: dt     = -1           ! time step width
  real(RNP) :: t_end  =  0           ! final problem time
  real(RNP) :: c_conv = -1           ! max Courant number   (< 0 if unlimited)
  real(RNP) :: c_diff = -1           ! max diffusion number (< 0 if unlimited)
  integer   :: nt_max = huge(1)      ! max number of time steps

  namelist /time_integration/ t_end, dt, c_conv, c_diff, nt_max

  integer :: time_method = 1
  ! 0  SDC - original
  ! 1  Euler
  ! 2  BDF2
  ! 3  TR
  ! 4  Runge-Kutta
  ! 5  SDC(Eu,Eu)
  ! 6  SDC(RK,Eu)

  namelist /time_integration/ time_method


  ! flow problem ...............................................................

  class(FlowProblem), allocatable :: problem

  ! mesh .......................................................................

  real(RNP)           :: xo(3)       ! corner nearest to origin
  real(RNP)           :: dx(3)       ! spacing in directions 1:3
  logical             :: periodic(3) ! periodic directions set true
  type(MeshPartition) :: mesh        ! mesh partition

  ! operators ..................................................................

  type(FlowOperators)   :: flow_op

  ! standalone time integrator
  class(TimeIntegrator), allocatable      :: time_integrator
  type(TimeIntegrator_Euler_Options)      :: eu_opt
  type(TimeIntegrator_BDF2_Options)       :: bdf2_opt
  type(TimeIntegrator_TR_Options)         :: tr_opt
  type(TimeIntegrator_RungeKutta_Options) :: rk_opt
  namelist /time_integration/ eu_opt, bdf2_opt, tr_opt, rk_opt

  ! original SDC
  type(SDC_Method3D)  :: sdc_orig
  type(SDC_Options3D) :: sdc_orig_opt
  namelist /time_integration/ sdc_orig_opt

  ! SDC
  class(Flow_SDC_Method), allocatable :: sdc
  type(Flow_SDC_Options_Euler) :: sdc_opt_eu
  namelist /time_integration/ sdc_opt_eu

  ! variables ..................................................................

  real(RNP), allocatable, target :: var(:,:,:,:,:)  ! storage for variables
  character(len=:),  allocatable :: name_var(:)     ! variable names

  ! solution variables !!! those listed in one line share storage !!!
  real(RNP), dimension(:,:,:,:,:), pointer, contiguous, save :: &
    u,            & ! approximate solution
    u_e, u_0,     & ! exact / initial solution
    err_u, w,     & ! error / workspace
    F,            & ! time derivative for SDC
    nu,           & ! diffusivity
    q_avg           ! temporally averaged quantities

  ! additional scalar variables
  real(RNP), dimension(:,:,:,:), pointer, contiguous, save :: &
    div_v,        & ! divergence of approximate velocity
    lmb_2           ! lambda_2 vortex sensor

  ! auxiliary ..................................................................

  logical   :: exists, first, last, failed
  integer   :: prm, stat, nt, nt_10, n_avg
  real(RDP) :: time, time0

  !-----------------------------------------------------------------------------
  ! initialization

  ! MPI ........................................................................

  call XMPI_Init()

  comm = MPI_COMM_WORLD
  call MPI_Comm_rank(comm, rank)
  call MPI_Comm_size(comm, n_proc)

  !$omp parallel
  n_thread = OMP_Num_Threads()
  !$omp end parallel

  ! control and discretization parameters ......................................

  if (rank == 0) then

    write(*,'(/,A)') repeat('=',80)
    write(*,'(A)') 'Cartesian IP/DG Incompressible Flow Solver'
    write(*,*)

    call get_command_argument(1, control_file, status=stat)
    if (stat /= 0 .or. len_trim(control_file) == 0) then
      control_file = 'isp_flow__sdc_test'
    end if

    inquire(file=trim(control_file)//'.prm', exist=exists)
    if (exists) then
      write(*,'(2X,A,/)') 'reading ' // trim(control_file)//'.prm'
      open(newunit=prm, file=trim(control_file)//'.prm')
      read(prm, nml=control)
      read(prm, nml=discretization)
      read(prm, nml=elliptic_solver)
      read(prm, nml=time_integration)
      close(prm)
    end if

    if (po_p < 1) po_p = po_u
    if (po_q < 1) po_q = po_u

    pmg_u_opt % po_top = po_u
    pmg_p_opt % po_top = po_p

    write(*,'(2X,A,1X,A)')  'flow type = ', trim(flow_type)
    write(*,'(2X,A,1X,A)')  'flow case = ', trim(flow_case)
    write(*,'(2X,A,1X,I0)') 'n_proc    = ', n_proc
    write(*,'(2X,A,1X,I0)') 'n_thread  = ', n_thread
    write(*,*)

  end if

  ! control parameters
  call XMPI_Bcast(flow_type   , 0, comm)
  call XMPI_Bcast(flow_case   , 0, comm)
  call XMPI_Bcast(restart     , 0, comm)
  call XMPI_Bcast(compute_p   , 0, comm)
  call XMPI_Bcast(write_stat  , 0, comm)
  call XMPI_Bcast(write_vtk   , 0, comm)
  call XMPI_Bcast(eval_rms    , 0, comm)
  call XMPI_Bcast(eval_max    , 0, comm)
  call XMPI_Bcast(eval_eps    , 0, comm)
  call XMPI_Bcast(average     , 0, comm)
  call XMPI_Bcast(i_avg       , 0, comm)
  call flow_op_control % Bcast( 0, comm)

  ! discretization and time integration parameters
  call XMPI_Bcast(np            , 0, comm)
  call XMPI_Bcast(ep            , 0, comm)
  call XMPI_Bcast(po_u          , 0, comm)
  call XMPI_Bcast(po_p          , 0, comm)
  call XMPI_Bcast(po_q          , 0, comm)
  call XMPI_Bcast(penalty       , 0, comm)
  call XMPI_Bcast(adjust_dx     , 0, comm)

  ! elliptic solver parameters
  call pmg_u_opt % Bcast(0, comm)
  call pmg_p_opt % Bcast(0, comm)

  ! time integration parameters
  call XMPI_Bcast(time_method, 0, comm)
  call eu_opt       % Bcast(0, comm)
  call bdf2_opt     % Bcast(0, comm)
  call tr_opt       % Bcast(0, comm)
  call rk_opt       % Bcast(0, comm)
  call sdc_orig_opt % Bcast(0, comm)
  call sdc_opt_eu   % Bcast(0, comm)

  ! flow problem ...............................................................

  select case(flow_type)
  case('Channel')
    allocate(FlowProblem_Channel           :: problem)
  case('No_Flow')
    allocate(FlowProblem_No_Flow           :: problem)
  case('Poiseuille')
    allocate(FlowProblem_Poiseuille        :: problem)
  case('Stokes_DKM')
    allocate(FlowProblem_Stokes_DKM        :: problem)
  case('Stokes_GMS')
    allocate(FlowProblem_Stokes_GMS        :: problem)
  case('Transition_TG')
    allocate(FlowProblem_Transition_TG     :: problem)
  case('Vortex_HW')
    allocate(FlowProblem_Vortex_HW         :: problem)
  case('Vortex_TG')
    allocate(FlowProblem_Vortex_TG         :: problem)
  case('VortexSheet')
    allocate(FlowProblem_VortexSheet       :: problem)
  case('VariableViscosity')
    allocate(FlowProblem_VariableViscosity :: problem)
  case default
    allocate(FlowProblem_Vortex_HW         :: problem)
  end select

  call problem % SetProblem(flow_case, comm)

  ! mesh .......................................................................

  ! origin and spacing
  xo = problem % x0
  dx = (problem % x1 - problem % x0) / (np * ep)

  if (adjust_dx) then
    dx(3) = maxval(dx(1:2))
  end if

  ! periodicity
  periodic(1) = all(problem % bc(1:2,:) == 'P')
  periodic(2) = all(problem % bc(3:4,:) == 'P')
  periodic(3) = all(problem % bc(5:6,:) == 'P')

  ! mesh partition and points
  call GenerateStructuredMesh(mesh, np, ep, xo, dx, periodic, comm)

  ! operators ..................................................................

  flow_op = FlowOperators( problem, mesh              &
                         , po_u, po_p, po_q, penalty  &
                         , pmg_u_opt                  &
                         , pmg_p_opt                  &
                         , flow_op_control            &
                         )

  select case(time_method)
  case(0)
    sdc_orig = SDC_Method3D(EulerVC, EulerVC, sdc_orig_opt)
  case(1)
    time_integrator = TimeIntegrator_Euler(problem, flow_op, eu_opt)
  case(2)
    time_integrator = TimeIntegrator_BDF2(problem, flow_op, bdf2_opt)
  case(3)
    time_integrator = TimeIntegrator_TR(problem, flow_op, tr_opt)
  case(4)
    time_integrator = TimeIntegrator_RungeKutta(problem, flow_op, rk_opt)
  case(5)
    sdc = Flow_SDC_Method_Euler(problem, flow_op, eu_opt, sdc_opt_eu)
  case(6)
    sdc = Flow_SDC_Method_Euler(problem, flow_op, rk_opt, sdc_opt_eu)
  end select

  ! print time-integrator settings
  if (rank == 0) then
    select case(time_method)
    case(1:4)
      call time_integrator % Show()
    case(5:7)
      call sdc % Show()
    end select
  end if

  ! variables and initial values ...............................................

  !$omp parallel

  call InitializeMeshVariables()

  if (restart) then
    call LoadFlowVariables(t, u, n_avg, q_avg)
  else
    t = 0
    call problem % GetInitialValues(flow_op%x, u)
    if (.not. problem % HasExactSolution()) then
      ! using p = u(*,4) and F_v = u_e
      call TimeDerivative(problem, flow_op, t, u, F = u_e)
      call ComputePressure(problem, flow_op, t, u_e, u(:,:,:,:,4), w) ! p = u(*,4)
    end if
  end if

  ! time stepping ..............................................................

  call SetTimeStep(problem, flow_op, u, c_conv, c_diff, dt)
  !$omp master
  nt_max = min(nint((t_end - t)/ dt, IXL), nt_max)
  call XMPI_Bcast(nt_max, 0, comm)
  nt_10 = int((nt_max + 9)/10)
  !$omp end master

  !-----------------------------------------------------------------------------
  ! Time integration

  failed = .false.

  if (rank == 0) then
    !$omp master
    time0 = MPI_Wtime()
    !$omp end master
  end if

  !$omp barrier

  do nt = 1, nt_max

    first = nt == 1
    last  = nt == nt_max

    if (write_stat) then
      call Evaluation(failed)
      if (failed) exit
    else
      !$omp master
      if (mod(nt, nt_10) == 0 .and. mesh%part == 0) then
        print '(I5,A)', nint(100.*nt/nt_max), ' %'
      end if
      !$omp end master
    end if

    select case(time_method)
    case(0)
      call sdc_orig % TimeStep( problem, flow_op, t, dt, u, F, first, last)
    case(1:4)
      call time_integrator % TimeStep(t, dt, u, standby = .not. last)
    case(5:7)
      call sdc % TimeStep(t, dt, u, standby = .not. last )
    case default
      ! using old Euler as the fall-back
      call SetArray(u_0, u, multi=.true.)
      call EulerVC(problem, flow_op, t, dt, u_0, u)
    end select

    if (compute_p) then
      associate(p => u(:,:,:,:,4), F_v => u_e)
        call TimeDerivative(problem, flow_op, t, u, F = F_v)
        call ComputePressure(problem, flow_op, t, F_v, p, w)
      end associate
    end if

  end do

  if (rank == 0) then
    !$omp master
    time = MPI_Wtime()
    time = time - time0
    !$omp end master
  end if

  if (.not. failed) then
    call Evaluation(last=.true.)
  else if (rank == 0) then
    !$omp master
    write(*,'(A)') 'failed'
    !$omp end master
  end if

  if (rank == 0) then
    !$omp master
    nt = max(min(nt, nt_max),1)
    write(*,'(A)') 'Performance'
    write(*,'(A,ES10.3)') ' time               =', time
    write(*,'(A,ES10.3)') ' time       / step  =', time / nt
    write(*,'(A,ES10.3)') ' throughput / step  =', nt / time * size(u) * n_proc
    !$omp end master
  end if

  !-----------------------------------------------------------------------------
  ! Save current solution

  call SaveFlowVariables(t, u, n_avg, q_avg)

  !-----------------------------------------------------------------------------
  ! Export results

  if (write_vtk) then

    call problem % GetDiffusivity(flow_op%x, t, u, nu)
    call GetDerivedVariables(flow_op, u, div_v, lmb_2)

    !$omp barrier
    !$omp master
    if (mesh%part >= 0) then
      call ExportVolumeDataToVTK( po_u                 & ! polynomial order
                                , mesh%ne              & ! number of elements
                                , size(name_var)       & ! number of scalars
                                , 0                    & ! no vectors
                                , flow_op%x            & ! mesh points
                                , var                  & ! variables
                                , name_var             & ! variable names
                                , file   = flow_case   & ! VTK file base name
                                , part   = mesh%part   & ! partition ID
                                , n_part = mesh%n_part ) ! number of partitions
    end if
    !$omp end master

  end if

  !$omp end parallel

  !-----------------------------------------------------------------------------
  ! Finalization

  call MPI_Finalize()

contains

!===============================================================================
! Internal procedures
!
! The following procedures perform individual tasks, mostly relying on host
! association. Arguments are used if appropriate. When mature, the procedures
! should be moved into a separate module or attached to a matching type.

!-------------------------------------------------------------------------------
!> Initialization of workspace

subroutine InitializeMeshVariables()

  ! local variables ............................................................

  character(len=:), allocatable, save :: name_u(:)
  integer :: np, nc, ne, n_var
  integer :: i, j

  ! prerequisites ..............................................................

  np = po_u + 1
  ne = mesh % ne
  nc = problem % nc

  ! number of variables
  n_var = 5 * nc + 2
  if (average) then
    n_var = n_var + 9
  end if

  !$omp single

  ! names of solution variables
  call problem % GetVariableNames(name_u)

  allocate(character(len=len(name_u) + 10) :: name_var(n_var))
  do i = 1, size(name_var)
    write(name_var(i), '(A,1X,I0)') 'var', i
  end do

  allocate(var(np, np, np, ne, n_var))

  !$omp end single

  call SetArray(var, ZERO, multi=.true.)

  ! approximate solution .......................................................

  !$omp single

  i = 0
  u(0:, 0:, 0:, 1:, 1:) => var(:,:,:,:,i+1:i+nc)
  do j = 1, nc
    name_var(i + j) = trim(name_u(j))
  end do

  ! exact / initial solution ...................................................

  i = i + nc
  u_e(0:, 0:, 0:, 1:, 1:) => var(:,:,:,:,i+1:i+nc)
  do j = 1, nc
    name_var(i + j) = trim(name_u(j)) // '__exact'
  end do

  u_0(0:, 0:, 0:, 1:, 1:) => u_e

  ! error ......................................................................

  i = i + nc
  err_u(0:, 0:, 0:, 1:, 1:) => var(:,:,:,:,i+1:i+nc)
  do j = 1, nc
    name_var(i + j) = 'error(' // trim(name_u(j)) // ')'
  end do

  ! workspace
  w(0:, 0:, 0:, 1:, 1:) => var(:,:,:,:,i+1:i+nc)

  ! diffusivity ................................................................

  i = i + nc
  nu(0:, 0:, 0:, 1:, 1:) => var(:,:,:,:,i+1:i+nc)
  do j = 1, nc
    write(name_var(i + j), '(A,I0)') 'nu_', j
  end do

  ! time derivative / diffusivity ..............................................

  i = i + nc
  F(0:, 0:, 0:, 1:, 1:) => var(:,:,:,:,i+1:i+nc)
  do j = 1, nc
    write(name_var(i + j), '(A,I0)') 'du/dt_', j
  end do

  ! vorticity and divergence ...................................................

  i = i + nc

  ! divergence
  div_v(0:, 0:, 0:, 1:) => var(:,:,:,:,i+1)
  name_var(i+1) = 'div(v)'

  ! lambda_2
  lmb_2(0:, 0:, 0:, 1:) => var(:,:,:,:,i+2)
  name_var(i+2) = 'lambda_2'

  ! temporal averages ..........................................................

  i = i + 2
  if (average) then
    q_avg(0:, 0:, 0:, 1:, 1:) => var(:,:,:,:,i+1:i+9)
    name_var(i+1) = '< v_1 >'
    name_var(i+2) = '< v_2 >'
    name_var(i+3) = '< v_3 >'
    name_var(i+4) = '< v_1 v_1 >'
    name_var(i+5) = '< v_1 v_2 >'
    name_var(i+6) = '< v_1 v_3 >'
    name_var(i+7) = '< v_2 v_2 >'
    name_var(i+8) = '< v_2 v_3 >'
    name_var(i+9) = '< v_3 v_3 >'
  else
    q_avg => null()
  end if

  !$omp end single

end subroutine InitializeMeshVariables

!-------------------------------------------------------------------------------
!> Evaluation of the approximate solution

subroutine Evaluation(failed, last)
  logical, optional, intent(out) :: failed
  logical, optional, intent(in)  :: last

  ! local variables ............................................................

  ! shared with OpenMP
  real(RNP), save :: err_v_max = -1, err_v_loc
  real(RNP), save :: err_p_max = -1, err_p_loc
  real(RNP), save :: div_v_max = -1, div_v_loc
  real(RNP), save :: e_kin     = -1, e_kin_loc
  real(RNP), save :: eps       = -1, eps_loc
  real(RNP), save :: eps_svv   = -1, eps_svv_loc
  integer,   save :: n         = -1, n_loc
  logical,   save :: head = .true.

  real(RNP), allocatable, save :: grad_v(:,:,:,:,:,:)

  real(RNP) :: err_v_rms, err_p_rms, div_v_rms

  logical   :: failed_, failed_nan, failed_max
  integer   :: e, i, j, k

  err_v_rms = -1
  err_p_rms = -1
  div_v_rms = -1

  associate( eop   => flow_op % eop_u     &
           , x     => flow_op % x         &
           , err_v => err_u(:,:,:,:,1:3)  &
           , err_p => err_u(:,:,:,:,4)    &
           )

    ! preliminaries ............................................................

    ! total number of mesh points
    !$omp master
    n_loc = size(x) / 3
    call XMPI_Allreduce(n_loc, n, MPI_SUM, mesh%comm)
    !$omp end master

    ! exact solution and error .................................................

    if (problem % HasExactSolution()) then

      call problem % GetExactSolution(x, t, u_e)

      ! error
      call SetArray(err_u, u, multi=.true.)                 ! err_u = u
      call MergeArrays(ONE, err_u, -ONE, u_e, multi=.true.) ! err_u = err_u - u_e

      ! remove constant from pressure error
      call CalibrateArray(err_p, mesh%comm)

      ! rms error
      if (eval_rms) then
        err_v_rms = sqrt(ScalarProduct(err_v, err_v, mesh%comm) / n)
        err_p_rms = sqrt(ScalarProduct(err_p, err_p, mesh%comm) / n)
      else
        err_v_rms = -1
        err_p_rms = -1
      end if

      ! max error
      !$omp master
      if (eval_max) then
        err_v_loc = maxval(abs(err_v))
        err_p_loc = maxval(abs(err_p))
        call XMPI_Reduce(err_v_loc, err_v_max, MPI_MAX, 0, mesh%comm)
        call XMPI_Reduce(err_p_loc, err_p_max, MPI_MAX, 0, mesh%comm)
      else
        err_v_max = -1
        err_p_max = -1
      end if
      !$omp end master

    end if

    ! divergence ...............................................................

    call WeakDivergence(mesh, eop%w, eop%D, u, div_v)

    ! rms divergence
    if (eval_rms) then
      div_v_rms = sqrt(ScalarProduct(div_v, div_v, mesh%comm) / n)
    end if

    ! max divergence
    !$omp master
    if (eval_max) then
      div_v_loc = maxval(abs(div_v))
      call XMPI_Reduce(div_v_loc, div_v_max, MPI_MAX, 0, mesh%comm)
    else
      div_v_max = -1
    end if
    !$omp end master

    ! kinetic energy ...........................................................

    !$omp master
    e_kin_loc = 0
    !$omp end master

    associate(po => eop%po, Ms => eop%w)
      !$omp do reduction(+:e_kin_loc)
      do e = 1, size(u,4)
        do k = 0, po
        do j = 0, po
        do i = 0, po
          e_kin_loc = e_kin_loc  &
                    + Ms(i)*Ms(j)*Ms(k) * ( u(i,j,k,e,1) ** 2 &
                                          + u(i,j,k,e,2) ** 2 &
                                          + u(i,j,k,e,3) ** 2 )
        end do
        end do
        end do
      end do
    end associate

    !$omp master
    e_kin_loc = product(mesh%dx) / 16 * e_kin_loc
    call XMPI_Reduce(e_kin_loc, e_kin, MPI_SUM, 0, mesh%comm)
    !$omp end master

    ! dissipation -- for constant viscosity only ...............................

    if (eval_eps) then
      associate(Ms => eop%w)
        block
          real(RNP), allocatable :: Ds_svv(:,:)
          real(RNP) :: nu, nu_svv
          integer   :: po, np, ne

          po = eop%po
          np = po + 1
          ne = size(u,4)
          nu = problem % nu_ref(1)

          !$omp single
          allocate(grad_v(0:po,0:po,0:po,ne,3,3))
          !$omp end single

          ! physical dissipation . . . . . . . . . . . . . . . . . . . . . . . .

          ! element velocity gradient
          call TPO_Grad_R( eop%D, mesh%dx, u(:,:,:,:,1), grad_v(:,:,:,:,1:3,1) )
          call TPO_Grad_R( eop%D, mesh%dx, u(:,:,:,:,2), grad_v(:,:,:,:,1:3,2) )
          call TPO_Grad_R( eop%D, mesh%dx, u(:,:,:,:,3), grad_v(:,:,:,:,1:3,3) )

          ! local element contributions to physical dissipation
          !$omp master
          eps_loc = 0
          !$omp end master
          !$omp do reduction(+:eps_loc)
          do e = 1, ne
            do k = 0, po
            do j = 0, po
            do i = 0, po
              eps_loc = eps_loc &
                      + Ms(i)*Ms(j)*Ms(k) * nu * sum(grad_v(i,j,k,e,:,:)**2)
            end do
            end do
            end do
          end do

          !$omp master
          eps_loc = product(mesh%dx)/8 * eps_loc
          call XMPI_Reduce(eps_loc, eps, MPI_SUM, 0, mesh%comm)
          !$omp end master

          ! SVV dissipation  . . . . . . . . . . . . . . . . . . . . . . . . . .

          if (flow_op_control%svv) then
            allocate(Ds_svv(0:po,0:po))

            ! get modified standard diff matrix and SVV coefficient
            if (eop % Has_SVV()) then
              call eop % Get_SVV_StandardRootDiffMatrix(Ds_svv) ! Ds_svv = √Q D
              nu_svv = problem % nu_svv_ref(1)
            else
              Ds_svv = ZERO
              nu_svv = ZERO
            end if

            ! element velocity gradient based on modified diff matrix
            call TPO_Grad_R(Ds_svv, mesh%dx, u(:,:,:,:,1), grad_v(:,:,:,:,1:3,1))
            call TPO_Grad_R(Ds_svv, mesh%dx, u(:,:,:,:,2), grad_v(:,:,:,:,1:3,2))
            call TPO_Grad_R(Ds_svv, mesh%dx, u(:,:,:,:,3), grad_v(:,:,:,:,1:3,3))

            ! local element contributions to SVV dissipation
            eps_svv_loc = 0
            !$omp do reduction(+:eps_svv_loc)
            do e = 1, ne
              do k = 0, po
              do j = 0, po
              do i = 0, po
                eps_svv_loc = eps_svv_loc + Ms(i)*Ms(j)*Ms(k) * nu_svv *   &
                                            sum(grad_v(i,j,k,e,:,:)**2)
              end do
              end do
              end do
            end do

            !$omp master
            eps_svv_loc = product(mesh%dx)/8 * eps_svv_loc
            call XMPI_Reduce(eps_svv_loc, eps_svv, MPI_SUM, 0, mesh%comm)
            !$omp end master
          end if
          !$omp master
          deallocate(grad_v)
          !$omp end master
        end block
      end associate

    else
      !$omp master
      eps     = -1
      eps_svv = -1
      !$omp end master
    end if

  end associate

  ! print statistics ...........................................................

  !$omp master
  if (rank == 0) then

    if (head) then

      write(*,'(A,  6X)',advance='NO') '#'
      write(*,'(A, 11X)',advance='NO') 't'
      write(*,'(A, 11X)',advance='NO') 'dt'
      write(*,'(A, 11X)',advance='NO') 'dx'
      write(*,'(A, 11X)',advance='NO') 'dy'
      write(*,'(A,  7X)',advance='NO') 'dz'
      write(*,'(A,  3X)',advance='NO') 'po'
      if (err_v_rms >= 0) write(*,'(A, 3X)',advance='NO') 'err_v_rms '
      if (err_p_rms >= 0) write(*,'(A, 3X)',advance='NO') 'err_p_rms '
      if (div_v_rms >= 0) write(*,'(A, 3X)',advance='NO') 'div_v_rms '
      if (err_v_max >= 0) write(*,'(A, 3X)',advance='NO') 'err_v_max '
      if (err_p_max >= 0) write(*,'(A, 3X)',advance='NO') 'err_p_max '
      if (div_v_max >= 0) write(*,'(A, 3X)',advance='NO') 'div_v_max '
      if (e_kin     >= 0) write(*,'(A, 3X)',advance='NO') '  e_kin   '
      if (eps       >= 0) write(*,'(A, 3X)',advance='NO') '   eps    '
      if (eps_svv   >= 0) write(*,'(A, 3X)',advance='NO') ' eps_svv  '
      write(*,*)

      head = .false.
    end if

    write(*,'(ES12.5,1X)',advance='NO') t
    write(*,'(ES12.5,1X)',advance='NO') dt
    write(*,'(ES12.5,1X)',advance='NO') dx(1)
    write(*,'(ES12.5,1X)',advance='NO') dx(2)
    write(*,'(ES12.5,1X)',advance='NO') dx(3)
    write(*,'(I4    ,1X)',advance='NO') po_u
    if (err_v_rms >= 0) write(*,'(ES12.5,1X)',advance='NO') err_v_rms
    if (err_p_rms >= 0) write(*,'(ES12.5,1X)',advance='NO') err_p_rms
    if (div_v_rms >= 0) write(*,'(ES12.5,1X)',advance='NO') div_v_rms
    if (err_v_max >= 0) write(*,'(ES12.5,1X)',advance='NO') err_v_max
    if (err_p_max >= 0) write(*,'(ES12.5,1X)',advance='NO') err_p_max
    if (div_v_max >= 0) write(*,'(ES12.5,1X)',advance='NO') div_v_max
    if (e_kin     >= 0) write(*,'(ES12.5,1X)',advance='NO') e_kin
    if (eps       >= 0) write(*,'(ES12.5,1X)',advance='NO') eps
    if (eps_svv   >= 0) write(*,'(ES12.5,1X)',advance='NO') eps_svv
    if (present(last)) then
      if (last) write(*,'(A)') ' #last#'
    end if
    write(*,*)

  end if
  !$omp end master

  ! check for fatal errors .....................................................

  !$omp master
  failed_nan = ieee_is_nan(err_v_rms)
  failed_max = abs(err_v_max) > huge(ONE)

  if (failed_nan) then
    call Warning('*** NaN detected ***','ISP_Flow__SDC_Test')
  else if (failed_max) then
    call Warning('*** err_v_max = INFINITY ***','ISP_Flow__SDC_Test')
  end if

  if (present(failed)) then
    failed_ = failed_nan .or. failed_max .or. err_v_max > problem%v_ref
    call XMPI_Allreduce(failed_, failed, MPI_LOR, mesh%comm)
  end if
  !$omp end master

end subroutine Evaluation

!-------------------------------------------------------------------------------
!> Saves the flow variables to disk

subroutine SaveFlowVariables(t, u, n_avg, q_avg)
  real(RNP),           intent(in) :: t                !< time
  real(RNP),           intent(in) :: u(:,:,:,:,:)     !< flow variables
  integer  ,           intent(in) :: n_avg            !< current num averages
  real(RNP), optional, intent(in) :: q_avg(:,:,:,:,:) !< averaged quantities

  integer :: unit
  character(len=100) :: file

  !$omp single

  write(file,'(2A,I0,A)') trim(flow_case), '_p', rank, '.fld'
  open(newunit=unit, file=file, status='REPLACE', form='UNFORMATTED')
  write(unit) t
  write(unit) u
  if (present(q_avg)) then
    write(unit) n_avg
    write(unit) q_avg
  else
    write(unit) 0
  end if
  close(unit)

  !$omp end single

end subroutine SaveFlowVariables

!-------------------------------------------------------------------------------
!> Loads the flow variables to disk

subroutine LoadFlowVariables(t, u, n_avg, q_avg)
  real(RNP),           intent(out) :: t                !< time
  real(RNP),           intent(out) :: u(:,:,:,:,:)     !< flow variables
  integer  ,           intent(out) :: n_avg            !< current num averages
  real(RNP), optional, intent(out) :: q_avg(:,:,:,:,:) !< averaged quantities

  integer :: unit, ios
  character(len=100) :: file

  !$omp single

  write(file,'(2A,I0,A)') trim(flow_case), '_p', rank, '.fld'
  open(newunit=unit, file=file, status='OLD', form='UNFORMATTED', iostat=ios)
  if (ios == 0) then
    read(unit) t
    read(unit) u
    read(unit) n_avg
    if (present(q_avg)) then
      if (n_avg > 0) then
        read(unit) q_avg
      else
        call SetArray(q_avg, ZERO, multi = .true.)
      end if
    end if
    close(unit)
  else
    call Error('SaveFlowVariables', 'found no matching input file')
  end if

  !$omp end single

end subroutine LoadFlowVariables

!===============================================================================
! Standalone procedures

!-------------------------------------------------------------------------------
!> Determination of the time step size
!>
!> To allow for explicit prescription, any given dt > 0 will not be changed.

subroutine SetTimeStep(problem, flow_op, u, c_conv, c_diff, dt)
  use Eigenproblems

  class(FlowProblem),   intent(in)    :: problem      !< flow problem
  class(FlowOperators), intent(in)    :: flow_op      !< flow operators
  real(RNP),            intent(in)    :: u(:,:,:,:,:) !< flow variables
  real(RNP),            intent(in)    :: c_conv       !< max Courant number
  real(RNP),            intent(in)    :: c_diff       !< max diffusion number
  real(RNP),            intent(inout) :: dt           !< time step size

  complex(RDP), allocatable :: Lambda_A(:)
  real(RDP)   , allocatable :: Lambda_L(:), A(:,:), L(:,:)

  real(RNP) :: dt_conv_lmb, dt_conv_po
  real(RNP) :: dt_diff_lmb, dt_diff_po
  real(RNP) :: cn, dn, v_max, v_max_loc, nu_max, dx_min, h

  !$omp master

  associate( mesh => flow_op % mesh  &
           , eop  => flow_op % eop_u &
           , po   => flow_op % po_u  &
           , v1   => u(:,:,:,:,1)    &
           , v2   => u(:,:,:,:,2)    &
           , v3   => u(:,:,:,:,3)    )

    dx_min = minval(mesh % dx)

    ! convective time scale.....................................................

    if (problem % stokes) then

      dt_conv_lmb = huge(dt_conv_lmb)
      dt_conv_po  = huge(dt_conv_po)

    else

      if (mesh % ne > 0) then
        v_max_loc = sqrt(maxval(v1**2 + v2**2 + v3**2))
      else
        v_max_loc = 0
      end if
      v_max_loc = max(v_max_loc, problem % v_ref)

      call XMPI_Reduce(v_max_loc, v_max, MPI_MAX, 0, mesh%comm)

      if (mesh %part == 0) then

        ! time scale based on eigenvalues of the 1D convection operator
        allocate(A(po,po), Lambda_A(po))
        A = eop % D(1:,1:)
        call SolveNonsymmetricEigenproblem(A, Lambda_A)
        h = dx_min / (2 * maxval(abs(Lambda_A)))
        dt_conv_lmb = h / v_max

        ! time scale based on average point spacing
        h = dx_min / po
        dt_conv_po = h / v_max

      end if

    end if

    ! diffusive time scale .....................................................

    nu_max = maxval(problem % nu_ref)
    if (allocated(problem % nu_svv_ref)) then
      nu_max = nu_max + maxval(problem % nu_svv_ref)
    end if

    if (mesh %part == 0) then

      ! time scale based on eigenvalues of the 1D laplacian
      allocate(L(po-1,po-1), Lambda_L(po-1))
      L = eop % L(1:po-1,1:po-1)
      call SolveSymmetricEigenproblem(L, Lambda_L)
      dt_diff_lmb = 1 / (4 * nu_max * maxval(abs(Lambda_L)) * sum(1/(dx*dx)))

      ! time scale based on average point spacing
      dt_diff_po = 1 / (4 * nu_max * po**2 * sum(1/(dx*dx)))

    end if

    ! evaluation ...............................................................

    if (mesh % part == 0) then

      if (dt <= 0) then
        dt = huge(dt)
        if (c_conv > 0 .and. .not. problem % stokes) then
          dt = min(dt, c_conv * dt_conv_lmb)
        end if
        if (c_diff > 0) then
          dt = min(dt, c_diff * dt_diff_lmb)
        end if
      end if

      if (.not. problem % stokes) then
        cn = dt / dt_conv_lmb
        write(*,'(A)', advance='NO') &
          '# Courant   number  based on eigenvalues / mean spacing  C ='
        if (cn > 0.01 .and. cn < 100) then
          write(*,'(F7.3,A,F7.3)') cn, '  / ', dt / dt_conv_po
        else
          write(*,'(ES10.3,A,ES10.3)') cn ,'  / ', dt / dt_conv_po
        end if
      end if

      if (dt_diff_lmb > 0) then
        dn = dt / dt_diff_lmb
        write(*,'(A)',advance='NO') &
          '# Diffusion number  based on eigenvalues / mean spacing  D ='
        if (dn > 0.01 .and. dn < 100) then
          write(*,'(F7.3,A,F7.3)') dn, '  / ', dt / dt_diff_po
        else
          write(*,'(ES10.3,A,ES10.3)') dn, dt / dt_diff_po
        end if
      end if
      write(*,*)

    end if

    call XMPI_Bcast(dt, 0, comm)

  end associate

  !$omp end master

end subroutine SetTimeStep

!-------------------------------------------------------------------------------
!>

subroutine GetDerivedVariables(flow_op, u, div_v, lmb_2)
  class(FlowOperators), intent(in)  :: flow_op        !< flow operators
  real(RNP),            intent(in)  :: u(:,:,:,:,:)   !< flow variables
  real(RNP),            intent(out) :: div_v(:,:,:,:) !< divergence of velocity
  real(RNP),            intent(out) :: lmb_2(:,:,:,:) !< lambda_2 vortex sensor

  ! local variables ............................................................

  real(RNP), allocatable, save :: grad_v(:,:,:,:,:,:)
  real(RNP) :: grad_ve(3,3)
  integer   :: po, np, ne
  integer   :: e, i, j, k

  associate( mesh => flow_op % mesh          &
           , Ms   => flow_op % eop_u % w     &
           , Ds   => flow_op % eop_u % D     )

    ! intialization ............................................................

    po = ubound(Ms,1)
    np = po + 1
    ne = mesh % ne

    !$omp single
    allocate(grad_v(np,np,np,ne,3,3))
    !$omp end single

    ! grad_v(:,i) = ∇vᵢ
    do i = 1, 3
      call WeakGradient(mesh, Ms, Ds, u(:,:,:,:,i), grad_v(:,:,:,:,1:3,i))
    end do

    ! divegence and lambda2 ....................................................

    !$omp do
    do e = 1, ne
      do k = 1, np
      do j = 1, np
      do i = 1, np
        grad_ve = grad_v(i,j,k,e,1:3,1:3)
        div_v(i,j,k,e) = grad_ve(1,1) + grad_ve(2,2) + grad_ve(3,3)
        lmb_2(i,j,k,e) = Lambda2(grad_ve)
      end do
      end do
      end do
    end do

    ! clean-up .................................................................

    !$omp barrier
    !$omp master
    deallocate(grad_v)
    !$omp end master

  end associate

end subroutine GetDerivedVariables

!-------------------------------------------------------------------------------
!> Evaluates the λ₂ vortex sensor based on the velocity gradient tensor
!>
!> author:  Immo Huismann, Joerg Stiller

real(RNP) function Lambda2(grad_v)
  use Eigenproblems, only: SolveSymmetricEigenproblem
  real(RNP), intent(in) :: grad_v(3,3) !< gradient of velocity ∇v

  real(RNP) :: S(3,3)     ! symmetric part of grad_v
  real(RNP) :: O(3,3)     ! antisymmetric part of grad_v
  real(RDP) :: SS_OO(3,3) ! S² + O²
  real(RDP) :: lambda(3)  ! eigenvalues of S² + O²

  S = HALF * (grad_v + transpose(grad_v))
  O = HALF * (grad_v - transpose(grad_v))

  SS_OO = matmul(S,S) + matmul(O,O)

  call SolveSymmetricEigenproblem(SS_OO, lambda)

  Lambda2 = lambda(2)

end function Lambda2

!===============================================================================
! Temporal averaging

!-------------------------------------------------------------------------------
!> Perform temporal averaging

subroutine TemporalAveraging

  real(RNP) :: wa, ws
  integer   :: e, i, j, k

  if (.not. average .or. mod(nt, i_avg) > 0) return

  n_avg = n_avg + 1
  if (n_avg > 1) then
    wa = real(n_avg - 1, RNP) / n_avg   ! weight of current average
    ws = ONE - wa                       ! weight of sample
  else
    wa = 0
    ws = 1
  end if

  !$omp do
  do e = 1, mesh%ne
    do concurrent(i = 0:po_u, j = 0:po_u, k = 0:po_u)
      q_avg(i,j,k,e,1) = wa * q_avg(i,j,k,e,1) + ws * u(i,j,k,e,1)
      q_avg(i,j,k,e,2) = wa * q_avg(i,j,k,e,2) + ws * u(i,j,k,e,2)
      q_avg(i,j,k,e,3) = wa * q_avg(i,j,k,e,3) + ws * u(i,j,k,e,3)
      q_avg(i,j,k,e,4) = wa * q_avg(i,j,k,e,4) + ws * u(i,j,k,e,1) * u(i,j,k,e,1)
      q_avg(i,j,k,e,5) = wa * q_avg(i,j,k,e,5) + ws * u(i,j,k,e,1) * u(i,j,k,e,2)
      q_avg(i,j,k,e,6) = wa * q_avg(i,j,k,e,6) + ws * u(i,j,k,e,1) * u(i,j,k,e,3)
      q_avg(i,j,k,e,7) = wa * q_avg(i,j,k,e,7) + ws * u(i,j,k,e,2) * u(i,j,k,e,2)
      q_avg(i,j,k,e,8) = wa * q_avg(i,j,k,e,8) + ws * u(i,j,k,e,2) * u(i,j,k,e,3)
      q_avg(i,j,k,e,9) = wa * q_avg(i,j,k,e,9) + ws * u(i,j,k,e,3) * u(i,j,k,e,3)
    end do
  end do

end subroutine TemporalAveraging

!===============================================================================

end program ISP_Flow__SDC_Test
