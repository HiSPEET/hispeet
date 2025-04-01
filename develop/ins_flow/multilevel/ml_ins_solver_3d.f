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
  use ML__INS__Operator__3D

  implicit none

  !-----------------------------------------------------------------------------
  ! Declarations

  ! MPI and OpenMP variables ...................................................

  type(MPI_Comm) :: comm      ! MPI communicator
  integer        :: rank      ! local MPI rank
  integer        :: n_proc    ! number of MPI processes
  integer        :: n_thread  ! number of OpenMP threads

  ! control parameters .........................................................

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

  integer :: char_freq  = 1        ! characteristics output frequency

  namelist/control_prm/ char_freq

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

  ! spatial parameters .........................................................

  type(ML_Mesh_Options_3D), save :: ml_mesh_opt
    ! multilevel mesh options, defining top level, partitioning and refinement,
    ! stacic components defined in namelist `ml_mesh_options_3d__static` and
    ! dynamic components in `ml_mesh_options_3d__dynamic`

  type(ML_INS_OperatorOptions_3D), save :: ml_ins_opt

  namelist/spatial_prm/ ml_ins_opt

  ! mesh .......................................................................

  type(GenericMesh_3D),       save :: generic_mesh
  type(Mesh_3D), allocatable, save :: base_mesh
  type(ML_Mesh_3D),           save :: ml_mesh

  ! problem, operators and solvers .............................................

  class(INS_Problem_3D), allocatable, save :: problem
  type(ML_INS_Operator_3D), save :: ml_ins

  ! auxiliaries ................................................................

  character(:), allocatable :: domain_name

  logical   :: exists, restart_in, restart_out
  integer   :: io, stat
  integer   :: n_bound
  integer   :: l_top

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

  end if

  l_top = size(ml_mesh % mesh)

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

  allocate(ml_ins_opt % po(l_top))

  ! read
  if (rank == 0) then
    write(*,'(2X,A)') 'creating spatial operators'
    open(newunit = io, file = case_file)
    read(io, nml = spatial_prm)
    close(io)
  end if

  ! globalize
  call ml_ins_opt % Bcast(0, comm)

  ml_ins = ML_INS_Operator_3D(ml_ins_opt, ml_mesh, problem)

  !-----------------------------------------------------------------------------
  ! Finalization

  call MPI_Finalize()

  !=============================================================================

end program ML_INS_Solver_3D
