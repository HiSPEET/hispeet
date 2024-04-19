!> summary:  Validation of incompressible Navier-Stokes time integrators
!> author:   Joerg Stiller
!> date:     2022/10/19
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> If present, the first argument of the invoking command will be interpreted
!> as the base name of the control file. If omitted, the program looks for
!> `ins_timeintegrator_3d_test.prm`.
!===============================================================================

program INS_TimeIntegrator_3D_Test
  use Kind_Parameters
  use Constants
  use OpenMP_Binding
  use XMPI
  use Execution_Control
  use Array_Assignments
  use Array_Reductions
  use Logging_Levels

  use Mesh__3D
  use Boundary_Variable__3D
  use Trace_Operators__3D
  use Volume_Integrals__3D
  use Surface_Integrals__3D
  use Export_VTK_Volume_Data__3D

  use INS__Problem__3D
  use INS__Problem__Test_Suite__3D

  use INS__Operator__3D
  use INS__Time_Scales__3D
  use INS__Time_Integrator__3D
  use INS__Time_Integrator__Euler__3D
  use INS__Time_Integrator__BDF2__3D
  use INS__Time_Integrator__Runge_Kutta__3D
  use INS__Flow_Characteristics__3D

  use Create_Cuboid_Cartesian
  use Create_Cuboid_Diamonds
  use Create_Cylinder
  use Create_Annulus

  use Generic_Mesh__3D
  use Verify_Mesh__3D
  use Import_GMSH__3D

  use Root_Mesh_Partitioning__3D

  implicit none

  !-----------------------------------------------------------------------------
  ! Declarations

  ! MPI and OpenMP variables ...................................................

  type(MPI_Comm) :: comm      ! MPI communicator
  integer        :: rank      ! local MPI rank
  integer        :: n_proc    ! number of MPI processes
  integer        :: n_thread  ! number of OpenMP threads

  ! control parameters .........................................................

  character(len=*), parameter :: default_case = 'ins_timeintegrator_3d_test'
  character(len=80) :: flow_case ! flow case name
  character(len=80) :: case_file ! flow case input file: trim(flow_case).prm
  ! input file (*.prm)

  character(len=80) :: flow_problem = 'Vortex_TG' ! problem name
  character(len=80) :: problem_file = 'vortex_tg' ! problem parameters file

  integer :: flow_domain = 1
  ! computational flow domain (u/s = un/structured, r = regular, d = deformed)
  !   1  cuboidal domain with Cartesian mesh                               (s+r)
  !   2  cuboidal domain with unstructured "diamond" mesh                  (u+d)
  !   3  cylindrical domain                                                (u+d)
  !   4  annular domain                                                    (u+d)
  !  10  import from GMSH

  character(len=80) :: raw_mesh_file = ''

  namelist/control_prm/ flow_problem, problem_file, flow_domain, raw_mesh_file

  integer :: time_method = 1
  ! 1  Euler
  ! 2  BDF2
  ! 3  Runge-Kutta

  namelist/control_prm/ time_method

  type(INS_OperatorOptions_3D) :: ins_op_opts
  type(INS_TimeIntegrator_Euler_Options_3D)      :: ins_ti_euler_opts
  type(INS_TimeIntegrator_BDF2_Options_3D)       :: ins_ti_bdf2_opts
  type(INS_TimeIntegrator_RungeKutta_Options_3D) :: ins_ti_runge_kutta_opts

  namelist/control_prm/ ins_op_opts,            &
                        ins_ti_euler_opts,      &
                        ins_ti_bdf2_opts,       &
                        ins_ti_runge_kutta_opts

  real(RNP) :: t_end    = 1  ! final time
  real(RNP) :: dt       = 1  ! time step size
  integer   :: nt_max   = 0  ! max num time steps

  namelist/control_prm/ t_end, dt, nt_max

  logical :: export_vtk = .false.  ! generate VTK files
  integer :: char_freq  = 1        ! characteristics output frequency

  namelist/control_prm/ export_vtk, char_freq

  ! control of logging levels
  namelist/control_prm/ log_level
  namelist/control_prm/ log_level_inner_iteration
  namelist/control_prm/ log_level_outer_iteration

  ! operators and variables ....................................................

  type(GenericMesh_3D) :: generic_mesh
  ! intermediate generic mesh for importing raw meshes

  type(Mesh_3D), allocatable, save :: initial_mesh
  ! mesh before repartitioning

  type(PartitioningOptions_3D) :: part_opt

  class(INS_Problem_3D), allocatable, save :: problem
  ! flow problem

  type(INS_Operator_3D), save :: ins_op
  ! incompressible Navier-Stokes operator

  class(INS_TimeIntegrator_3D), allocatable, save :: ins_ti
  ! incompressible Navier-Stokes time integrator

  type(INS_TimeScales_3D)          :: time_scales
  type(INS_FlowCharacteristics_3D) :: flow_char

  real(RNP) :: t = 0 ! problem time

  real(RNP), allocatable, target, save :: var(:,:,:,:,:)
  character(len=20), allocatable, save :: var_name(:)

  real(RNP), pointer, contiguous, save :: u(:,:,:,:,:)    ! u = [v, p]
  real(RNP), pointer, contiguous, save :: v(:,:,:,:,:)    ! velocity
  real(RNP), pointer, contiguous, save :: p(:,:,:,:)      ! pressure

  real(RNP), pointer, contiguous, save :: u_ex(:,:,:,:,:) ! u_ex = [v_ex, p_ex]
  real(RNP), pointer, contiguous, save :: v_ex(:,:,:,:,:) ! exact velocity
  real(RNP), pointer, contiguous, save :: p_ex(:,:,:,:)   ! exact pressure

  real(RNP), pointer, contiguous, save :: err_u(:,:,:,:,:) ! error, u - u_ex
  real(RNP), pointer, contiguous, save :: err_v(:,:,:,:,:) ! velocity error
  real(RNP), pointer, contiguous, save :: err_p(:,:,:,:)   ! pressure error

  real(RNP), allocatable, save :: w(:,:,:,:,:)   ! workspace

  type(BoundaryVariable_3D), allocatable, save :: bv_vn(:) ! n⋅v on Γ=∂Ω
  real(RNP), allocatable, save :: int_vn(:,:) ! ∫n⋅v dΓ

  ! auxiliaries ................................................................

  character(len=80) :: domain_name = ''
! real(RDP) :: time, time0
  real(RNP) :: domain_volume
  logical   :: exists, last, passed
  integer   :: io, stat
  integer   :: n_bound, n_elem, n_elem_tot, n_ghost, n_point, n_var, po
  integer   :: i, nt

  !-----------------------------------------------------------------------------
  ! Initialization

  ! MPI and OpenMP .............................................................

  call XMPI_Init()

  comm = MPI_COMM_WORLD
  call MPI_Comm_rank(comm, rank)
  call MPI_Comm_size(comm, n_proc)

  !$omp parallel
  n_thread = OMP_Num_Threads()
  !$omp end parallel

  ! parameters .................................................................

  ! read control parameters
  if (rank == 0) then

    write(*,'(/,A)') repeat('=',80)
    write(*,'(A)') 'Validation of incompressible Navier-Stokes time integrators'
    write(*,*)
    write(*,'(T3,A,T30,9(G0,X))') 'number of processes:'  , n_proc
    write(*,'(T3,A,T30,9(G0,X))') 'number of threads:'    , n_thread
    write(*,*)

    call get_command_argument(1, flow_case, status=stat)
    if (stat /= 0 .or. len_trim(flow_case) == 0) then
      flow_case = default_case
    end if
    case_file = trim(flow_case) // '.prm'

    inquire(file=case_file, exist=exists)
    if (exists) then
      write(*,'(2X,A)') 'reading ' // trim(case_file)
      open(newunit = io, file = case_file)
      read(io, nml = control_prm)
      close(io)
    else
       call Error( 'INS_TimeIntegrator_3D_Test', &
                   'input file "' // trim(case_file) // '" not found' )
    end if

  end if

  ! globalize control parameters
  call XMPI_Bcast(flow_case   , 0, comm)
  call XMPI_Bcast(case_file   , 0, comm)
  call XMPI_Bcast(flow_problem, 0, comm)
  call XMPI_Bcast(problem_file, 0, comm)
  call XMPI_Bcast(flow_domain , 0, comm)
  call XMPI_Bcast(time_method , 0, comm)
  call XMPI_Bcast(t_end       , 0, comm)
  call XMPI_Bcast(dt          , 0, comm)
  call XMPI_Bcast(nt_max      , 0, comm)
  call XMPI_Bcast(export_vtk  , 0, comm)
  call XMPI_Bcast(char_freq   , 0, comm)

  ! globalize logging levels
  call XMPI_Bcast_LoggingLevels(0, comm)

  ! globalize options
  call ins_op_opts             % Bcast(0, comm)
  call ins_ti_euler_opts       % Bcast(0, comm)
  call ins_ti_bdf2_opts        % Bcast(0, comm)
  call ins_ti_runge_kutta_opts % Bcast(0, comm)

  ! mesh .......................................................................

  allocate(initial_mesh)

  select case(flow_domain)
  case(1)
    call CreateCuboidCartesian(comm, case_file, initial_mesh)
    domain_name = 'Cuboidal domain with Cartesian mesh'
  case(2)
    call CreateCuboidDiamonds(comm, case_file, initial_mesh)
    domain_name = 'Cuboidal domain with unstructured "diamond" mesh'
  case(3)
    call CreateCylinder(comm, case_file, initial_mesh)
    domain_name = 'Cylindrical domain with unstructured mesh'
  case(4)
    call CreateAnnulus(comm, case_file, initial_mesh)
    domain_name = 'Annular domain with unstructured mesh'
  case(10)
    if (rank == 0) then
      call ImportGMSH_3D(raw_mesh_file, generic_mesh)
    end if
    call initial_mesh % ImportGenericMesh(generic_mesh, comm)
    domain_name = raw_mesh_file
  end select

  ! root mesh partitioning .....................................................

  if (initial_mesh % n_parts /= n_proc) then
    part_opt = PartitioningOptions_3D(n_parts = n_proc, w_comp = [1,1,0,0])
    call RootMeshPartitioning_3D(part_opt, initial_mesh, ins_op % mesh)
  else
    ins_op % mesh = initial_mesh
  end if

  deallocate(initial_mesh)

  call VerifyMesh_3D(ins_op % mesh, passed)
  if (rank == 0) then
    write(*,'(/,A)') 'Verification of computational mesh'
  end if
  if (ins_op % mesh % part >= 0) then
    write(*,'(2X,A,I4,A,L1)') 'part:',ins_op%mesh%part,': passed = ',passed
  end if

  n_elem  = ins_op % mesh % n_elem
  n_ghost = ins_op % mesh % n_ghost
  n_bound = ins_op % mesh % n_bound

  ! problem ....................................................................

  call Set_INS_TestProblem_3D(problem, flow_problem, problem_file, n_bound, comm)

  ! check & fix boundary conditions
  do i = 1, n_bound
    if (ins_op % mesh % boundary(i) % coupled > 0) then
      if (problem % bc_v(i) /= 'P') then
        if (rank == 0) then
          write(*,'(A,I0)') '  *** enforcing periodic BC on coupled boundary ',i
        end if
        problem % bc_v(i) = 'P'
      end if
    end if
  end do

  ! operators ..................................................................

  ! Navier-Stokes operator
  call ins_op % Init(ins_op_opts, problem)

  ! time integrator
  select case(time_method)
  case(1)
    ins_ti = INS_TimeIntegrator_Euler_3D(problem, ins_op, ins_ti_euler_opts)
  case(2)
    ins_ti = INS_TimeIntegrator_BDF2_3D(problem, ins_op, ins_ti_bdf2_opts)
  case(3)
    ins_ti = INS_TimeIntegrator_RungeKutta_3D &
                 (problem, ins_op, ins_ti_runge_kutta_opts)
  end select

  ! variables ..................................................................

  po = ins_op % eop_v % po
  n_var = 12

  allocate(var(0:po,0:po,0:po,1:n_elem,1:n_var), source = ZERO)
  allocate(var_name(1:n_var))

  u(0:,0:,0:,1:,1:)  =>  var(:,:,:,:,1:4)
  v(0:,0:,0:,1:,1:)  =>  var(:,:,:,:,1:3)
  p(0:,0:,0:,1:)     =>  var(:,:,:,:,4)

  var_name(1:4) = [ 'v_x', 'v_y', 'v_z', 'p  ']

  u_ex(0:,0:,0:,1:,1:)  =>  var(:,:,:,:,5:8)
  v_ex(0:,0:,0:,1:,1:)  =>  var(:,:,:,:,5:7)
  p_ex(0:,0:,0:,1:)     =>  var(:,:,:,:,8)

  var_name(5:8) = [ 'v_x__exact', 'v_y__exact', 'v_z__exact', 'p__exact  ']

  err_u(0:,0:,0:,1:,1:)  =>  var(:,:,:,:,9:12)
  err_v(0:,0:,0:,1:,1:)  =>  var(:,:,:,:,9:11)
  err_p(0:,0:,0:,1:)     =>  var(:,:,:,:,12)

  var_name(9:12) = [ 'error(v_x)', 'error(v_y)', 'error(v_z)', 'error(p)  ']

  allocate(w(0:po,0:po,0:po,1:n_elem,1:4) )

  ! initial conditions .........................................................

  call problem % GetInitialValues(ins_op % sem_v % metrics % x, u)

  ! time scales
  call time_scales % Evaluate(problem, ins_op, u)

  ! info .......................................................................

  call ins_op % sem_v % Get_Volume(domain_volume)

  call XMPI_Reduce(n_elem, n_elem_tot, MPI_SUM, 0, comm)
  n_point = n_elem_tot * (po+1)**3

  if (rank == 0) then
    write(*,'(/,A)') 'problem and discretization parameters'
    write(*,'(T3,A,T30,9(G0,X))') 'flow problem:', trim(flow_problem)
    write(*,'(T3,A,T30,9(G0,X))') 'domain:',  trim(domain_name)
    write(*,'(T3,A,T30,9(G0,X))') 'domain volume:', domain_volume
    write(*,'(T3,A,T30,9(G0,X))') 'boundary conditions:'  , problem % bc_v
    write(*,'(T3,A,T30,9(G0,X))') 'polynomial order of v:', ins_op % eop_v % po
    write(*,'(T3,A,T30,9(G0,X))') 'polynomial order of p:', ins_op % eop_p % po
    write(*,'(T3,A,T30,9(G0,X))') 'conv quadrature order:', ins_op % sop_q % po
    write(*,'(T3,A,T30,9(G0,X))') 'conv quadrature type:' , ins_op % sop_q % basis
    write(*,'(T3,A,T30,9(G0,X))') 'number of mesh points:', n_point
    write(*,'(T3,A,T30,9(G0,X))') 'time integrator:'      , trim(ins_ti % name)
    write(*,'(T3,A,T29,ES18.11)') 'time step size:'       , dt
    write(*,'((T7,A,T29,ES12.5,2X,A,ES12.5,X,A))')  &
            'dt / tau_c(v_0)'   , dt / time_scales % tau_conv_ve     , &
                            '(' , dt / time_scales % tau_conv_vm,')' , &
            'dt / tau_c(v_ref)' , dt / time_scales % tau_conv_re     , &
                            '(' , dt / time_scales % tau_conv_rm,')' , &
            'dt / tau_d(nu_ref)', dt / time_scales % tau_diff_re     , &
                             '(', dt / time_scales % tau_diff_rm,')'
    write(*,*)
  end if

  !-----------------------------------------------------------------------------
  ! Time integration

  if (char_freq > 0) then
    call flow_char % Evaluate(problem, ins_op, t, u, domain_volume)
    call flow_char % PrintHeader()
    call flow_char % PrintValues()
  end if

  do nt = 1, nt_max
    last = t + dt >= t_end .or. nt == nt_max
    call ins_ti % TimeStep(t, dt, u, standby = .not. last)
    if (mod(nt, char_freq) == 0) then
      call flow_char % Evaluate(problem, ins_op, t, u, domain_volume)
      call flow_char % PrintValues()
    end if
    if (last) exit
  end do

  !-----------------------------------------------------------------------------
  ! evaluation

  ! errors .....................................................................

  if (ins_op % mesh % part >= 0) then

    call problem % GetExactSolution(ins_op % sem_v % metrics % x, t, u_ex)

    call SetArray(err_u, u, multi=.true.)
    call MergeArrays(ONE, err_u, -ONE, u_ex, multi=.true.)
    call CalibrateArray(err_p, comm = ins_op % mesh % comm_parts)

  end if

  ! boundary fluxes ............................................................

  if (ins_op % mesh % part >= 0) then

    !$omp master
    allocate(bv_vn(n_bound), int_vn(1,n_bound))
    !$omp end master
    !$omp barrier

    do i = 1, n_bound
      call bv_vn(i) % Create(ins_op % mesh % boundary(i), po, nc=1)
      call bv_vn(i) % ExtractNormalComponent(ins_op % sem_v, v)
    end do

    call GetSurfaceIntegrals(ins_op % sem_v, bv_vn, int_vn)

  end if

  if (ins_op % mesh % part == 0) then
    !$omp master
    write(*,'(/,A)') 'boundary fluxes'
    do i = 1, n_bound
       write(*,'(T3,A,I3,A,T29,ES12.5)') &
           'boundary',i,',   int(vn)  =', int_vn(1,i)
    end do
    write(*,'(T3,A,T29,ES12.5)') 'total:     sum(int(vn)) =', sum(int_vn)
    write(*,*)
    !$omp end master
  end if

  !-----------------------------------------------------------------------------
  ! Write plot files

  if (export_vtk .and. ins_op % mesh%part >= 0) then
    call ExportVTK_VolumeData( x       = ins_op % sem_v % metrics % x &
                             , s       = var                          &
                             , sname   = var_name                     &
                             , file    = flow_case                    &
                             , part    = ins_op % mesh % part         &
                             , n_parts = ins_op % mesh % n_parts      )
  end if

  !-----------------------------------------------------------------------------
  ! Finalization

  call MPI_Finalize()

end program INS_TimeIntegrator_3D_Test
