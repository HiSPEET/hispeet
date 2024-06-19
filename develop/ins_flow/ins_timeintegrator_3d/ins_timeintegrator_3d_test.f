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

  ! restart options
  character(len=80) :: restart_tag_in  = ''  ! tag for restart input files
  character(len=80) :: restart_tag_out = ''  ! tag for restart output files
  namelist/control_prm/ restart_tag_in, restart_tag_out
  !
  ! restart input is read from
  !   - trim(flow_case)_trim(restart_tag_in)_mesh_<rank>.h5  for the mesh
  !   - trim(flow_case)_trim(restart_tag_in)_data_<rank>.h5  for flow data
  ! where <rank> is the process rank in mesh%comm_world
  !
  ! output written to
  !   - trim(flow_case)_trim(restart_tag_out)_mesh_<rank>.h5  for the mesh
  !   - trim(flow_case)_trim(restart_tag_out)_data_<rank>.h5  for flow data
  !
  ! no restart data is read or written if the corresponding tag is empty
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

  character(:), allocatable :: domain_name
  character(:), allocatable :: mesh_file
  character(:), allocatable :: data_file

  real(RNP) :: domain_volume
  logical   :: exists, last, passed, restart_in, restart_out
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
  call XMPI_Bcast(flow_case      , 0, comm)
  call XMPI_Bcast(case_file      , 0, comm)
  call XMPI_Bcast(flow_problem   , 0, comm)
  call XMPI_Bcast(problem_file   , 0, comm)
  call XMPI_Bcast(flow_domain    , 0, comm)
  call XMPI_Bcast(time_method    , 0, comm)
  call XMPI_Bcast(t_end          , 0, comm)
  call XMPI_Bcast(dt             , 0, comm)
  call XMPI_Bcast(nt_max         , 0, comm)
  call XMPI_Bcast(export_vtk     , 0, comm)
  call XMPI_Bcast(char_freq      , 0, comm)
  call XMPI_Bcast(restart_tag_in , 0, comm)
  call XMPI_Bcast(restart_tag_out, 0, comm)

  ! globalize logging levels
  call XMPI_Bcast_LoggingLevels(0, comm)

  ! globalize options
  call ins_op_opts             % Bcast(0, comm)
  call ins_ti_euler_opts       % Bcast(0, comm)
  call ins_ti_bdf2_opts        % Bcast(0, comm)
  call ins_ti_runge_kutta_opts % Bcast(0, comm)

  ! restart switches
  restart_in  = len_trim(restart_tag_in)  > 0
  restart_out = len_trim(restart_tag_out) > 0

  ! domain name ................................................................

  select case(flow_domain)
  case(1)
    domain_name = 'Cuboidal domain with Cartesian mesh'
  case(2)
    domain_name = 'Cuboidal domain with unstructured "diamond" mesh'
  case(3)
    domain_name = 'Cylindrical domain with unstructured mesh'
  case(4)
    domain_name = 'Annular domain with unstructured mesh'
  case(10)
    domain_name = raw_mesh_file
  end select

  ! create mesh ................................................................

  if (restart_in) then

    mesh_file = trim(flow_case) // '_' // trim(restart_tag_in) // '_mesh'
    call ins_op % mesh % ReadHDF5(mesh_file, comm)

  else

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
      domain_name = trim(raw_mesh_file)
    end select

    ! mesh partitioning ........................................................

    if (initial_mesh % n_parts /= n_proc) then
      part_opt = PartitioningOptions_3D(n_parts = n_proc, w_comp = [1,1,0,0])
      call RootMeshPartitioning_3D(part_opt, initial_mesh, ins_op % mesh)
    else
      ins_op % mesh = initial_mesh
    end if

    deallocate(initial_mesh)

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

  if (restart_in) then
    data_file = trim(flow_case) // '_' // trim(restart_tag_in) // '_data'
    call ReadRestartData(data_file, rank, t, u)
    call XMPI_Bcast(t, 0, comm)
  else
    call problem % GetInitialValues(ins_op % sem_v % metrics % x, u)
    t = 0
  end if

  ! time scales
  call time_scales % Evaluate(problem, ins_op, u)

  ! info .......................................................................

  call ins_op % sem_v % Get_Volume(domain_volume)

  call XMPI_Reduce(n_elem, n_elem_tot, MPI_SUM, 0, comm)
  n_point = n_elem_tot * (po+1)**3

  if (rank == 0) then
    write(*,'(/,A)') 'problem and discretization parameters'
    write(*,'(T3,A,T30,9(G0,X))') 'flow problem:', trim(flow_problem)
    write(*,'(T3,A,T30,9(G0,X))') 'domain:',  domain_name
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
  ! Write restart data

  if (restart_out) then
    mesh_file = trim(flow_case) // '_' // trim(restart_tag_out) // '_mesh'
    data_file = trim(flow_case) // '_' // trim(restart_tag_out) // '_data'
    ! mesh
    call ins_op % mesh % WriteHDF5(mesh_file)
    ! flow data
    call WriteRestartData(data_file, rank, t, u)
  end if

  !-----------------------------------------------------------------------------
  ! Finalization

  call MPI_Finalize()

contains

  !-----------------------------------------------------------------------------
  !> Write restart data

  subroutine WriteRestartData(file, rank, t, u)
    use, intrinsic ::  ISO_C_Binding
    use Kind_Parameters
    use HDF5_Binding
    implicit none

    character(len=*),  intent(in) :: file              !< file base name
    integer,           intent(in) :: rank              !< process rank
    real(RNP), target, intent(in) :: t                 !< problem time
    real(RNP), target, intent(in) :: u(0:,0:,0:,:,:)   !< variables

    integer(hid_t)    :: data_id, file_id, group_id, space_id
    integer(hsize_t)  :: dims_t(1), dims_var(5)
    integer           :: err
    character(len=80) :: tag
    character(len=:), allocatable :: file_pr

    call Init_HDF5_Binding()

    ! append process rank to file name
    write(tag,'(I0)') rank
    file_pr = trim(file)//'_'//trim(tag)//'.h5'

    ! create HDF5 file and group
    call H5Fcreate_f(file_pr, H5F_ACC_TRUNC_F, file_id, err)
    call H5Gcreate_f(file_id, 'data', group_id, err)

    ! time
    dims_t = 1
    call H5Screate_simple_f(size(dims_t), dims_t, space_id, err)
    call H5Dcreate_f(group_id, 't', H5T_REAL_RNP, space_id, data_id, err)
    call H5Dwrite_f(data_id, H5T_REAL_RNP, C_Loc(t), err)
    call H5Sclose_f(space_id, err)
    call H5Dclose_f(data_id, err)

    ! u
    dims_var = shape(u)
    call H5Screate_simple_f(size(dims_var), dims_var, space_id, err)
    call H5Dcreate_f(group_id, 'u', H5T_REAL_RNP, space_id, data_id, err)
    call H5Dwrite_f(data_id, H5T_REAL_RNP, C_Loc(u), err)
    call H5Sclose_f(space_id, err)
    call H5Dclose_f(data_id, err)

    ! release resources
    call H5Gclose_f(group_id, err)
    call H5Fclose_f(file_id, err)

  end subroutine WriteRestartData

  !-----------------------------------------------------------------------------
  !> Read restart data

  subroutine ReadRestartData(file, rank, t, u)
    use, intrinsic ::  ISO_C_Binding
    use Kind_Parameters
    use Execution_Control
    use HDF5_Binding
    implicit none

    character(len=*),  intent(in)    :: file              !< file base name
    integer,           intent(in)    :: rank              !< process rank
    real(RNP), target, intent(inout) :: t                 !< problem time
    real(RNP), target, intent(inout) :: u(0:,0:,0:,:,:)   !< variables

    integer(hid_t)    :: data_id, file_id, group_id, space_id, type_id
    integer(hsize_t)  :: dims_var(5), maxdims_var(5)
    integer           :: err
    logical           :: exists
    type(C_Ptr)       :: buf
    character(len=80) :: tag
    character(len=:), allocatable :: file_pr

    call Init_HDF5_Binding()

    ! append process rank to file name and check if it exists
    write(tag,'(I0)') rank
    file_pr = trim(file)//'_'//trim(tag)//'.h5'
    inquire(file=file_pr, exist=exists)

    if (.not. exists) return

    ! open HDF5 file and group for reading
    call H5Fopen_f(file_pr, H5F_ACC_RDWR_F, file_id, err)
    call H5Gopen_f(file_id, 'data', group_id, err)

    ! time
    buf = C_Loc(t)
    call H5Dopen_f(group_id, 't', data_id, err)
    call H5Dget_type_f(data_id, type_id, err)
    call H5Dread_f(data_id, type_id, buf, err)
    call H5Tclose_f(type_id, err)
    call H5Dclose_f(data_id, err)

    ! u
    buf = C_Loc(u)
    call H5Dopen_f(group_id, 'u', data_id, err)
    call H5Dget_space_f(data_id, space_id, err)
    call H5Sget_simple_extent_dims_f(space_id, dims_var, maxdims_var, err)
    if (any(dims_var /= shape(u))) then
      call Error('ReadRestartData','shape(u) not matching')
    end if
    call H5Dget_type_f(data_id, type_id, err)
    call H5Dread_f(data_id, type_id, buf, err)
    call H5Tclose_f(type_id, err)
    call H5Dclose_f(data_id, err)

    ! release resources
    call H5Gclose_f(group_id, err)
    call H5Fclose_f(file_id, err)

  end subroutine ReadRestartData

  !=============================================================================

end program INS_TimeIntegrator_3D_Test
