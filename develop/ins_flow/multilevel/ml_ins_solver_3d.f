!> summary:  3D multilevel solver for incompressible Navier-Stokes problems
!> author:   Joerg Stiller
!> date:     2025/03/28
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> If present, the first argument of the invoking command will be interpreted
!> as the base name of the control file. If omitted, the program looks for
!> `ml_ins_solver_3d.prm`.
!===============================================================================

program ML_INS_Solver_3D
  use Kind_Parameters
  use Constants
  use OpenMP_Binding
  use XMPI
  use Execution_Control
  use Logging_Levels

  use Create_Cuboid_Cartesian
  use Create_Cuboid_Diamonds
  use Create_Cylinder
  use Create_Annulus

  use Mesh__3D
  use Generic_Mesh__3D
  use Verify_Mesh__3D
  use Import_GMSH__3D

  use INS__Problem__3D
  use INS__Problem__Test_Suite__3D

  use ML__Mesh__3D
  use ML__Mesh_Variable__3D
  use ML__INS__Operator__3D
  use ML__INS__Integrator__BDF2__3D

  implicit none

  !-----------------------------------------------------------------------------
  ! Declarations

  ! MPI and OpenMP variables ...................................................

  type(MPI_Comm) :: comm      ! MPI communicator
  integer        :: rank      ! local MPI rank
  integer        :: n_proc    ! number of MPI processes
  integer        :: n_thread  ! number of OpenMP threads

  ! control parameters .........................................................

  ! control of logging levels
  namelist/control_prm/ log_level
  namelist/control_prm/ log_level_inner_iteration
  namelist/control_prm/ log_level_outer_iteration
  namelist/control_prm/ log_level_multigrid_cycle

  character(len=*), parameter :: default_case = 'ml_ins_solver_3d'
  character(len=80) :: flow_case ! flow case name
  character(len=80) :: case_file ! flow case input file: trim(flow_case).prm

  character(len=80) :: flow_problem = 'VariableViscosity'  ! problem name
  character(len=80) :: problem_file = 'variable_viscosity' ! parameters file

  integer :: flow_domain = 1
    ! computational flow domain (u/s = un/structured, r = regular, d = deformed)
    !   1  cuboidal domain with Cartesian mesh                             (s+r)
    !   2  cuboidal domain with unstructured "diamond" mesh                (u+d)
    !   3  cylindrical domain                                              (u+d)
    !   4  annular domain                                                  (u+d)
    !  10  import from GMSH

  character(len=80) :: raw_mesh_file = ''

  namelist/control_prm/ flow_problem, problem_file, flow_domain, raw_mesh_file

  integer :: char_freq  = 1 ! characteristics output frequency
  integer :: vtk_mode   = 0 ! VTK export mode, 0/1/2/3: none/all/active/leafs

  namelist/control_prm/ char_freq, vtk_mode

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

  ! mesh .......................................................................

  type(GenericMesh_3D),       save :: generic_mesh
  type(Mesh_3D), allocatable, save :: base_mesh
  type(ML_Mesh_3D),           save :: ml_mesh

  type(ML_Mesh_Options_3D), save :: ml_mesh_opt
    ! multilevel mesh options, defining top level, partitioning and refinement,
    ! stacic components defined in namelist `ml_mesh_options_3d__static` and
    ! dynamic components in `ml_mesh_options_3d__dynamic`

  ! problem ....................................................................

  class(INS_Problem_3D), allocatable, save :: problem

  ! spatial operators ..........................................................

  integer, allocatable, save :: po(:)
  type(ML_INS_OperatorOptions_3D), save :: ml_ins_opt
  type(ML_INS_Operator_3D), save :: ml_ins

  namelist/spatial_prm/ po, ml_ins_opt

  ! time integration ...........................................................

  type(ML_INS_Integrator_BDF2_Options_3D), save :: ml_bdf2_opt
  type(ML_INS_Integrator_BDF2_3D), save :: ml_bdf2

  namelist/temporal_prm/ ml_bdf2_opt

  real(RNP) :: t_end  = 0.25
  real(RNP) :: dt     = 1E-3
  integer   :: nt_max = 1

  namelist/temporal_prm/ t_end, dt, nt_max

  ! variables ..................................................................

  real(RNP) :: t = 0 ! problem time

  type(ML_MeshVariable_3D), save :: var   ! container for essential variables
  type(ML_MeshVariable_3D), save :: u     ! handle for numerical solution
  type(ML_MeshVariable_3D), save :: u_ex  ! handle for exact solution
  type(ML_MeshVariable_3D), save :: err_u ! handle for error, err_u = u - u_ex

  ! auxiliaries ................................................................

  character(:), allocatable :: domain_name
  character(:), allocatable :: var_name(:)

  logical :: exists, restart_in, restart_out
  logical :: first, last
  integer :: io, stat
  integer :: l_top, n_bound, n_var
  integer :: i, l, nt

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
    write(*,'(A)') 'Multilevel Navier-Stokes solver for incompressible flow'
    write(*,*)
    write(*,'(T3,A,T30,9(G0,X))') 'number of processes:', n_proc
    write(*,'(T3,A,T30,9(G0,X))') 'number of threads:'  , n_thread
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
      if (char_freq < 1) then
        ! disable intermediate control output
        char_freq = huge(1)
      end if
    else
       call Error( 'ML_INS_Solver_3D', &
                   'input file "' // trim(case_file) // '" not found' )
    end if

  end if

  ! globalize control parameters
  call XMPI_Bcast(flow_case      , 0, comm)
  call XMPI_Bcast(case_file      , 0, comm)
  call XMPI_Bcast(flow_problem   , 0, comm)
  call XMPI_Bcast(problem_file   , 0, comm)
  call XMPI_Bcast(flow_domain    , 0, comm)
  call XMPI_Bcast(raw_mesh_file  , 0, comm)
  call XMPI_Bcast(char_freq      , 0, comm)
  call XMPI_Bcast(restart_tag_in , 0, comm)
  call XMPI_Bcast(restart_tag_out, 0, comm)

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

  if (restart_in) then

    ! TBD ......................................................................

  else

    ! create base mesh .........................................................

    allocate(base_mesh)

    select case(flow_domain)
    case(1)
      call CreateCuboidCartesian(comm, case_file, base_mesh)
      domain_name = 'Cuboidal domain with Cartesian mesh'
    case(2)
      call CreateCuboidDiamonds(comm, case_file, base_mesh)
      domain_name = 'Cuboidal domain with unstructured "diamond" mesh'
    case(3)
      call CreateCylinder(comm, case_file, base_mesh)
      domain_name = 'Cylindrical domain with unstructured mesh'
    case(4)
      call CreateAnnulus(comm, case_file, base_mesh)
      domain_name = 'Annular domain with unstructured mesh'
    case(10)
      if (rank == 0) then
        call ImportGMSH_3D(raw_mesh_file, generic_mesh)
      end if
      call base_mesh % ImportGenericMesh(generic_mesh, comm)
      domain_name = trim(raw_mesh_file)
    end select

    ! create multilevel mesh ...................................................

    if (rank == 0) then
      write(*,'(2X,A)') 'creating multilevel mesh'
      ! read multilevel mesh options
      open(newunit = io, file = case_file)
      ml_mesh_opt = ML_Mesh_Options_3D(io, n_proc)
      close(io)
    end if

    call ml_mesh_opt % Bcast(0, comm)
    ml_mesh = ML_Mesh_3D(base_mesh, ml_mesh_opt)
    deallocate(base_mesh)

  end if

  l_top   = size(ml_mesh % mesh)
  n_bound = ml_mesh % mesh(1) % n_bound

  ! problem ....................................................................

  call Set_INS_TestProblem_3D(problem, flow_problem, problem_file, n_bound, comm)

  ! check & fix boundary conditions
  do i = 1, n_bound
    if (ml_mesh % mesh(1) % boundary(i) % coupled > 0) then
      if (problem % bc_v(i) /= 'P') then
        if (rank == 0) then
          write(*,'(2X,A,X,I0)') '*** enforcing periodic BC on boundary', i
        end if
        problem % bc_v(i) = 'P'
      end if
    end if
  end do

  ! spatial ....................................................................

  allocate(po(l_top), source = -1)

  ! read
  if (rank == 0) then
    write(*,'(2X,A)') 'creating spatial operators'
    open(newunit = io, file = case_file)
    read(io, nml = spatial_prm)
    close(io)
  end if

  ! globalize
  call XMPI_Bcast(po, 0, comm)
  call ml_ins_opt % Bcast(0, comm)

  ml_ins = ML_INS_Operator_3D(ml_mesh, po, problem, ml_ins_opt)

  ! temporal ...................................................................

  ! read
  if (rank == 0) then
    write(*,'(2X,A)') 'creating multilevel integrator'
    open(newunit = io, file = case_file)
    read(io, nml = temporal_prm)
    close(io)
  end if

  ! globalize
  call ml_bdf2_opt % Bcast(0, comm)

  ml_bdf2 = ML_INS_Integrator_BDF2_3D(problem, ml_ins, ml_bdf2_opt)

  ! variables ..................................................................

  associate(nc => problem % nc)

    n_var = nc
    if (problem % HasExactSolution()) then
      n_var = n_var + 2 * nc
    end if

    allocate(character(len=20) :: var_name(n_var))

    ! names of solution components
    var_name(1:4) = [ 'v_x', 'v_y', 'v_z', 'p  ']
    do i = 5, nc
      write(var_name(i),'(A,I0)') 'u_', i
    end do

    ! names of exact solution and error components
    if (problem % HasExactSolution()) then
      do i = 1, nc
        write(var_name(i +   nc),'(2A)') trim(var_name(i)), '__exact'
        write(var_name(i + 2*nc),'(2A)') trim(var_name(i)), '__error'
      end do
    end if

    ! variable container
    call var % Init(ml_ins%ml_op_u, n_var, var_name)

    ! handles
    call var % GetSlice(u, first = 1, last = nc)
    if (problem % HasExactSolution()) then
      call var % GetSlice(u_ex , first = 1 +   nc, last = 2*nc)
      call var % GetSlice(err_u, first = 1 + 2*nc, last = 3*nc)
    end if

  end associate

  !-----------------------------------------------------------------------------
  ! Initial conditions

  if (restart_in) then

    ! TBD ......................................................................

  else

    do l = 1, l_top
      associate(ins_l => ml_ins % ins_op(l), u_l => u % level(l) % val)
        call problem % GetInitialValues(ins_l % sem_u % metrics % x, u_l)
      end associate
    end do

  end if

  !-----------------------------------------------------------------------------
  ! Time integration

  do nt = 1, nt_max
    first = nt == 1
    last  = t + dt >= t_end .or. nt == nt_max
    call ml_bdf2 % TimeStep(t, dt, u, first, last)
  end do

  !-----------------------------------------------------------------------------
  ! Write plot files

  if (vtk_mode > 0) then
    call u % ExportVTK(ml_ins%ml_op_u, file = flow_case, mode = vtk_mode)
  end if

  !-----------------------------------------------------------------------------
  ! Finalization

  call MPI_Finalize()

  !=============================================================================

end program ML_INS_Solver_3D
