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

  use Kind_Parameters, only: RNP, RDP
  use Constants,       only: ONE, ZERO, HALF
  use Array_Assignments
  use Array_Reductions
  use Standard_Operators_1D
  use Execution_Control
  use Export_Volume_Data_To_VTK
  use XMPI

  use ISP_Flow_Problem__Test_Suite

  use CART__TPO_Grad
  use CART__Mesh_Partition
  use CART__Generate_Structured_Mesh
  use CART__DG_Weak_Divergence
  use CART__DG_Weak_Gradient
  use CART__Elliptic_PMG
  use CART__ISP_Flow__Operators
  use CART__ISP_Flow__Time_Derivative
  use CART__ISP_Flow__Pressure

  ! time integration
  use CART__ISP_Flow__Euler                           ! old standalone Euler
  use CART__ISP_Flow__SDC                             ! current SDC
  use CART__ISP_Flow__Time_Integrator                 ! TI base type
  use CART__ISP_Flow__Time_Integrator__Euler_DS       ! new standalone Euler
  use CART__ISP_Flow__Time_Integrator__Runge_Kutta_DS ! new standalone Runge Kutta

  implicit none

  !-----------------------------------------------------------------------------
  ! Declarations

  ! MPI ........................................................................

  type(MPI_Comm) :: comm      ! communicator
  integer        :: comm_size ! number of processes
  integer        :: rank      ! local rank

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
  type(FlowOpControl) :: flow_op_control

  namelist /control/ flow_type, flow_case,         &
                     restart, compute_p,           &
                     write_stat, write_vtk,        &
                     eval_rms, eval_max, eval_eps, &
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

  type(SDC_Options3D) :: sdc_opt

  namelist /time_integration/ t_end, dt, c_conv, c_diff, nt_max, sdc_opt

  ! flow problem ...............................................................

  class(FlowProblem), allocatable :: problem

  ! mesh .......................................................................

  real(RNP)           :: xo(3)       ! corner nearest to origin
  real(RNP)           :: dx(3)       ! spacing in directions 1:3
  logical             :: periodic(3) ! periodic directions set true
  type(MeshPartition) :: mesh        ! mesh partition

  ! operators ..................................................................

  type(FlowOperators)   :: flow_op
  type(SDC_Method3D)    :: sdc

  ! standalone time integrator, switched on with sdc_opt%n_sub = 0
  class(TimeIntegrator), allocatable :: time_integrator

  ! variables ..................................................................

  real(RNP), allocatable, target :: var(:,:,:,:,:)  ! storage for variables
  character(len=:),  allocatable :: name_var(:)     ! variable names

  ! solution variables !!! those listed in one line share storage !!!
  real(RNP), dimension(:,:,:,:,:), pointer, contiguous :: &
    u,            & ! approximate solution
    u_e, u_0,     & ! exact / initial solution
    err_u, w,     & ! error / workspace
    F,            & ! time derivative for SDC
    nu              ! diffusivity

  ! additional scalar variables
  real(RNP), dimension(:,:,:,:), pointer, contiguous :: &
    div_v,        & ! divergence of approximate velocity
    lmb_2           ! lambda_2 vortex sensor

  ! auxiliary ..................................................................

  logical :: exists, first, last, failed
  integer :: prm, stat, nt, nt_10

  !-----------------------------------------------------------------------------
  ! initialization

  ! MPI ........................................................................

  call XMPI_Init()

  comm = MPI_COMM_WORLD
  call MPI_Comm_size(comm, comm_size)
  call MPI_Comm_rank(comm, rank)

  ! control and discretization parameters ......................................

  if (rank == 0) then

    call get_command_argument(1, control_file, status=stat)
    if (stat /= 0 .or. len_trim(control_file) == 0) then
      control_file = 'isp_flow__sdc_test'
    end if

    inquire(file=trim(control_file)//'.prm', exist=exists)
    if (exists) then
      write(*,'(/,2X,A,/)') 'reading ' // trim(control_file)//'.prm'
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

  ! SDC parameters
  call sdc_opt % Bcast(0, comm)

  ! flow problem ...............................................................

  select case(flow_type)
  case('Stokes_DKM')
    allocate(FlowProblem_Stokes_DKM        :: problem)
  case('Stokes_GMS')
    allocate(FlowProblem_Stokes_GMS        :: problem)
  case('Vortex_HW')
    allocate(FlowProblem_Vortex_HW         :: problem)
  case('Vortex_TG')
    allocate(FlowProblem_Vortex_TG         :: problem)
  case('Transition_TG')
    allocate(FlowProblem_Transition_TG     :: problem)
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

  select case (sdc_opt % n_sub)
  case(1:)
    sdc = SDC_Method3D(EulerVC, EulerVC, sdc_opt)
  case(0)
    time_integrator = TimeIntegrator_EulerDS(problem, flow_op)
  case(-1)      ! case time integration using RK- method
    time_integrator = TimeIntegrator_RungeKuttaDS(problem, flow_op)
!?  time_integrator = TimeIntegrator_RungeKuttaDS(problem, flow_op, ns, method)
  end select

  ! variables and initial values ...............................................

  call InitializeMeshVariables()

  if (restart) then
    call LoadFlowVariables(t, u)
  else
    t = 0
    call problem % GetInitialValues(flow_op%x, u)
  end if

  ! time stepping ..............................................................

  call SetTimeStep(problem, flow_op, u, c_conv, c_diff, dt)
  nt_max = min(nint((t_end - t)/ dt), nt_max)
  call XMPI_Bcast(nt_max, 0, comm)
  nt_10 = int((nt_max + 9)/10)

  !-----------------------------------------------------------------------------
  ! Time integration

  failed = .false.

  do nt = 1, nt_max

    first = nt == 1
    last  = nt == nt_max

    if (write_stat) then
      call Evaluation(failed)
      if (failed) exit
    else
      if (mod(nt, nt_10) == 0 .and. mesh%part == 0) then
        print '(I5,A)', nint(100.*nt/nt_max), ' %'
      end if
    end if

    select case (sdc % n_sub)
    case(1:)
      call sdc % TimeStep( problem, flow_op, t, dt, u, F, first, last)
    case(0)
      if (problem % HasVariableProperties()) then
        call problem % GetDiffusivity(flow_op%x, t, u, nu)
        call time_integrator % TimeStep(t, dt, u, nu)
      else
        call time_integrator % TimeStep(t, dt, u)
      end if
    case(-1)
      if (problem % HasVariableProperties()) then
        call problem % GetDiffusivity(flow_op%x, t, u, nu)
        call time_integrator % TimeStep(t, dt, u, nu)
      else
        call time_integrator % TimeStep(t, dt, u)
      end if
    case default
      call SetArray(u_0, u, multi=.true.)
      call EulerVC(problem, flow_op, t, dt, u_0, u)
    end select

    if (compute_p) then
      associate(p => u(:,:,:,:,4), F_v => u_e)
        call TimeDerivative(problem, flow_op, t, u_c = u, u_d = u, F = F_v)
        call PressureSolver(problem, flow_op, F_v, t, p, w)
      end associate
    end if

    if (last) exit

  end do

  if (.not. failed) then
    call Evaluation(last=.true.)
  else if (rank == 0) then
    !$omp master
    write(*,'(A)') 'failed'
    !$omp end master
  end if

  !-----------------------------------------------------------------------------
  ! Save current solution

  call SaveFlowVariables(t, u)

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

  character(len=:), allocatable :: name_u(:)
  integer :: np, nc, ne, n_var
  integer :: i, j

  ! prerequisites ..............................................................

  np = po_u + 1
  ne = mesh % ne
  nc = problem % nc

  ! number of variables
  n_var = 5 * nc + 2

  ! names of solution variables
  call problem % GetVariableNames(name_u)

  allocate(var(np, np, np, ne, n_var))
  call SetArray(var, ZERO, multi=.true.)

  allocate(character(len=len(name_u) + 10) :: name_var(n_var))
  do i = 1, size(name_var)
    write(name_var(i), '(A,1X,I0)') 'var', i
  end do

  ! approximate solution .......................................................

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

end subroutine InitializeMeshVariables

!-------------------------------------------------------------------------------
!> Evaluation of the approximate solution

subroutine Evaluation(failed, last)
  logical, optional, intent(out) :: failed
  logical, optional, intent(in)  :: last

  ! local variables ............................................................

  real(RNP) :: err_v_max = -1, err_v_rms = -1, err_v_loc
  real(RNP) :: err_p_max = -1, err_p_rms = -1, err_p_loc
  real(RNP) :: div_v_max = -1, div_v_rms = -1, div_v_loc
  real(RNP) :: e_kin, e_kin_loc
  real(RNP) :: eps = -1, eps_loc

  real(RNP), allocatable, save :: grad_v(:,:,:,:,:,:)
  logical   :: head = .true.
  logical   :: failed_
  integer   :: n = 0, n_loc = 0
  integer   :: e, i, j, k

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
      end if

      ! max error
      if (eval_max) then
        err_v_loc = maxval(abs(err_v))
        err_p_loc = maxval(abs(err_p))
        !$omp master
        call XMPI_Reduce(err_v_loc, err_v_max, MPI_MAX, 0, mesh%comm)
        call XMPI_Reduce(err_p_loc, err_p_max, MPI_MAX, 0, mesh%comm)
        !$omp end master
      end if

    end if

    ! divergence ...............................................................

    call WeakDivergence(mesh, eop%w, eop%D, u, div_v)

    ! rms divergence
    if (eval_rms) then
      div_v_rms = sqrt(ScalarProduct(div_v, div_v, mesh%comm) / n)
    end if

    ! max divergence
    if (eval_max) then
      div_v_loc = maxval(abs(div_v))
      !$omp master
      call XMPI_Reduce(div_v_loc, div_v_max, MPI_MAX, 0, mesh%comm)
      !$omp end master
    end if

    ! kinetic energy ...........................................................

    e_kin_loc = 0

    associate(po => eop%po, Ms => eop%w)
      !!$omp loop reduction
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

    e_kin_loc = product(mesh%dx) / 16 * e_kin_loc

    !$omp master
    call XMPI_Reduce(e_kin_loc, e_kin, MPI_SUM, 0, mesh%comm)
    !$omp end master

    ! dissipation -- for constant viscosity only ...............................

    if (eval_eps) then
      associate(Ms => eop%w)
        block
          real(RNP) :: nu
          integer   :: po, np, ne

          po = eop%po
          np = po + 1
          ne = size(u,4)
          nu = problem % nu_ref(1)

          !$omp single
          allocate(grad_v(0:po,0:po,0:po,ne,3,3))
          !$omp end single

          call TPO_Grad_Eval( np, ne, eop%D, mesh%dx, u(:,:,:,:,1), &
                              grad_v(:,:,:,:,1:3,1)                 )
          call TPO_Grad_Eval( np, ne, eop%D, mesh%dx, u(:,:,:,:,2), &
                              grad_v(:,:,:,:,1:3,2)                 )
          call TPO_Grad_Eval( np, ne, eop%D, mesh%dx, u(:,:,:,:,3), &
                              grad_v(:,:,:,:,1:3,3)                 )

          eps_loc = 0
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
          eps_loc = product(mesh%dx)/8 * eps_loc

          !$omp barrier
          !$omp master
          call XMPI_Reduce(eps_loc, eps, MPI_SUM, 0, mesh%comm)
          deallocate(grad_v)
          !$omp end master

        end block
      end associate
    end if

  end associate

  ! print statistics ...........................................................

  !$omp master
  if (rank == 0) then

    if (head) then

      write(*,'(A)')    '#'
      write(*,'(A,I0)') '# sdc_opt % n_sub   = ', sdc_opt  % n_sub
      write(*,'(A,I0)') '# n_sub   = ', sdc % n_sub
      write(*,'(A,I0)') '# n_sweep = ', sdc % n_sweep
      write(*,'(A)')    '#'

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
      if (err_p_max >= 0) write(*,'(A, 3X)',advance='NO') 'err_p_max '
      if (err_v_max >= 0) write(*,'(A, 3X)',advance='NO') 'err_v_max '
      if (div_v_max >= 0) write(*,'(A, 3X)',advance='NO') 'div_v_max '
      if (e_kin     >= 0) write(*,'(A, 3X)',advance='NO') '  e_kin   '
      if (eps       >= 0) write(*,'(A, 3X)',advance='NO') '   eps    '
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
    if (err_p_max >= 0) write(*,'(ES12.5,1X)',advance='NO') err_v_max
    if (err_v_max >= 0) write(*,'(ES12.5,1X)',advance='NO') err_p_max
    if (div_v_max >= 0) write(*,'(ES12.5,1X)',advance='NO') div_v_max
    if (e_kin     >= 0) write(*,'(ES12.5,1X)',advance='NO') e_kin
    if (eps       >= 0) write(*,'(ES12.5,1X)',advance='NO') eps
    if (present(last)) then
      if (last) write(*,'(A)') ' #last#'
    end if
    write(*,*)

  end if
  !$omp end master

  ! check for fatal errors .....................................................

  !$omp master
  if (present(failed)) then
    failed_ = err_v_max > problem%v_ref .or. ieee_is_nan(err_v_rms)
    call XMPI_Allreduce(failed_, failed, MPI_LOR, mesh%comm)
  end if
  !$omp end single

end subroutine Evaluation

!-------------------------------------------------------------------------------
!> Saves the flow variables to disk

subroutine SaveFlowVariables(t, u)
  real(RNP), intent(in) :: t            !< time
  real(RNP), intent(in) :: u(:,:,:,:,:) !< flow variables

  integer :: unit
  character(len=100) :: file

  write(file,'(2A,I0,A)') trim(flow_case), '_p', rank, '.fld'
  open(newunit=unit, file=file, status='REPLACE', form='UNFORMATTED')
  write(unit) t
  write(unit) u
  close(unit)

end subroutine SaveFlowVariables

!-------------------------------------------------------------------------------
!> Loads the flow variables to disk

subroutine LoadFlowVariables(t, u)
  real(RNP), intent(out) :: t            !< time
  real(RNP), intent(out) :: u(:,:,:,:,:) !< flow variables

  integer :: unit, ios
  character(len=100) :: file

  write(file,'(2A,I0,A)') trim(flow_case), '_p', rank, '.fld'
  open(newunit=unit, file=file, status='OLD', form='UNFORMATTED', iostat=ios)
  if (ios == 0) then
    read(unit) t
    read(unit) u
    close(unit)
  else
    call Error('SaveFlowVariables', 'found no matching input file')
  end if

end subroutine LoadFlowVariables

!===============================================================================
! Standalone procedures

!-------------------------------------------------------------------------------
!> Determination of the time step size
!>
!> To allow for explicit prescription, any given dt > 0 will not be changed.

subroutine SetTimeStep(problem, flow_op, u, c_conv, c_diff, dt)
  class(FlowProblem),   intent(in)    :: problem      !< flow problem
  class(FlowOperators), intent(in)    :: flow_op      !< flow operators
  real(RNP),            intent(in)    :: u(:,:,:,:,:) !< flow variables
  real(RNP),            intent(in)    :: c_conv       !< max Courant number
  real(RNP),            intent(in)    :: c_diff       !< max diffusion number
  real(RNP),            intent(inout) :: dt           !< time step size

  real(RNP) :: dt_conv, dt_diff, cn, dn
  real(RNP) :: v_max, v_max_loc
  real(RNP) :: nu_max

  associate( mesh => flow_op % mesh &
           , po   => flow_op % po_u &
           , v1   => u(:,:,:,:,1)   &
           , v2   => u(:,:,:,:,2)   &
           , v3   => u(:,:,:,:,3)   )

    ! convective time scale.....................................................

    if (problem % stokes) then

      dt_conv = huge(dt_conv)

    else

      if (mesh % ne < 0) then
        v_max_loc = sqrt(maxval(v1**2 + v2**2 + v3**2))
      else
        v_max_loc = 0
      end if
      v_max_loc = max(v_max_loc, problem % v_ref)

      call XMPI_Reduce(v_max_loc, v_max, MPI_MAX, 0, mesh%comm)

      if (mesh %part == 0) then
        dt_conv = minval(mesh % dx)/ (po * v_max)
      end if

    end if

    ! diffusive time scale .....................................................

    nu_max = maxval(problem % nu_ref)

    if (mesh %part == 0) then
      dt_diff = 1 / (2 * nu_max * po**2 * sum(1/(dx*dx)))
    end if

    ! evaluation ...............................................................

    if (mesh % part == 0) then

      if (dt <= 0) then
        dt = huge(dt)
        if (c_conv > 0 .and. .not. problem % stokes) then
          dt = min(dt, c_conv * dt_conv)
        end if
        if (c_diff > 0) then
          dt = min(dt, c_diff * dt_diff)
        end if
      end if

      if (.not. problem % stokes) then
        cn = dt / dt_conv
        write(*,'(A)',advance='NO') '# Courant number    C ='
        if (cn > 0.01 .and. cn < 100) then
          write(*,'(F7.3)') cn
        else
          write(*,'(ES10.3)') cn
        end if
      end if

      if (dt_diff > 0) then
        dn = dt / dt_diff
        write(*,'(A)',advance='NO') '# Diffusion number  D ='
        if (dn > 0.01 .and. dn < 100) then
          write(*,'(F7.3)') dn
        else
          write(*,'(ES10.3)') dn
        end if
      end if
      write(*,*)

    end if

    call XMPI_Bcast(dt, 0, comm)

  end associate

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

end program ISP_Flow__SDC_Test
