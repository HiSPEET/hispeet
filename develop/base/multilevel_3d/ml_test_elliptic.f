!> summary:  Program for testing the elliptic multilevel solvers
!> author:   Joerg Stiller
!> date:     2024/09/16
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

program ML_Test_Elliptic
  use Kind_Parameters
  use Constants
  use OpenMP_Binding
  use XMPI
  use Execution_Control
  use Logging_Levels

  use Elliptic_Problem__3D
  use Elliptic_Problem__Simple__3D
  use Elliptic_Problem__Knotty__3D

  use Create_Cuboid_Cartesian
  use Create_Cuboid_Diamonds
  use Create_Cylinder
  use Create_Annulus

  use Import_GMSH__3D
  use Generic_Mesh__3D
  use Mesh__3D
  use Verify_Mesh__3D

  use ML__Mesh__3D

  implicit none

  !-----------------------------------------------------------------------------
  ! Declarations

  ! MPI and OpenMP variables ...................................................

  type(MPI_Comm) :: comm      ! MPI communicator
  integer        :: rank      ! local MPI rank
  integer        :: n_proc    ! number of MPI processes
  integer        :: n_thread  ! number of OpenMP threads

  ! control parameters .........................................................

  ! NOTE
  ! The case name is defined by the first command line argument.
  ! If no argument is given, the default case is assumed.

  character(len=*), parameter :: default_case = 'ml_test_elliptic'
  character(len=80) :: case_name ! case name
  character(len=80) :: case_file ! case input file: trim(test_case).prm

  namelist/control_prm/ log_level
  namelist/control_prm/ log_level_inner_iteration
  namelist/control_prm/ log_level_outer_iteration
  namelist/control_prm/ log_level_multigrid_cycle

  logical :: export_vtk = .false.  ! switch for VTK export

  namelist/control_prm/ export_vtk

  ! domain parameters ..........................................................

  integer :: test_domain = 1 ! computational domain
                             !   0  import from GMSH
                             !   1  cuboidal with Cartesian mesh
                             !   2  cuboidal with unstructured "diamond" mesh
                             !   3  cylindrical domain

  character(len=80) :: gmsh_file = '../gmsh_3d/cylinder_2d'
  integer :: pg = 1 ! polynomial degree of element geometry

  namelist/domain_prm/ test_domain, gmsh_file, pg

  ! problem parameters .........................................................

  integer :: test_problem    =  3      ! 1/2/3/4: simple_{1/2/3}d / knotty
  logical :: has_variable_nu = .false. ! T/F: variable/constant ν

  namelist/problem_prm/ test_problem, has_variable_nu

  ! NOTE
  ! The fluctuation amplitude ν₁ is ignored in case of constant ν

  real(RNP) :: lambda = 0         ! Helmholtz parameter
  real(RNP) :: nu_0   = 1         ! diffusivity mean value ν₀
  real(RNP) :: nu_1   = 0         ! diffusivity fluctuation amplitude ν₁
  real(RNP) :: d_nu   = 0         ! diffusivity fluctuation phase shift
  integer   :: k_nu   = 1         ! diffusivity fluctuation wave number
  integer   :: k_u    = 1         ! solution wave number
  character, allocatable :: bc(:) ! boundary conditions {'D','N','P'} ['D']

  namelist/problem_prm/ lambda, nu_0, nu_1, d_nu, k_nu, k_u, bc

  ! mesh, operator, variables ..................................................

  type(GenericMesh_3D)     , save :: generic_mesh
  type(Mesh_3D)            , save :: base_mesh
  type(ML_Mesh_Options_3D) , save :: ml_mesh_opt
  type(ML_Mesh_3D)         , save :: ml_mesh
! type(ML_MeshOperators_3D), save :: ml_op
! type(ML_MeshVariable_3D) , save :: ml_var

  integer, allocatable :: po(:)    ! sequence of polynomial orders

  namelist/operator_prm/ po

  ! auxiliaries ................................................................

  character(:), allocatable :: domain_name

  logical   :: exists, passed, all_passed
  integer   :: io, stat
  integer   :: l_top, ne_max, ne_min, ne_tot
  integer   :: l

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

  ! control and domain parameters ..............................................

  if (rank == 0) then

    write(*,'(/,A)') repeat('=',80)
    write(*,'(A)') 'Validation of elliptic multilevel solvers'
    write(*,*)
    write(*,'(T3,A,T30,9(G0,X))') 'number of processes:', n_proc
    write(*,'(T3,A,T30,9(G0,X))') 'number of threads:'  , n_thread
    write(*,*)

    call get_command_argument(1, case_name, status=stat)
    if (stat /= 0 .or. len_trim(case_name) == 0) then
      case_name = default_case
    end if
    case_file = trim(case_name) // '.prm'

    inquire(file=case_file, exist=exists)
    if (exists) then
      write(*,'(2X,A)') 'reading ' // trim(case_file)
      open(newunit = io, file = case_file)
      read(io, nml = control_prm)
      read(io, nml = domain_prm)
      ml_mesh_opt = ML_Mesh_Options_3D(io, n_proc)
      allocate(po(ml_mesh_opt%l_top), source = -1)
      read(io, nml = operator_prm)
      close(io)
    else
      call Error( 'ML_Test_Elliptic', &
                  'input file "' // trim(case_file) // '" not found' )
      call Error('ML_Test_Elliptic', 'file "'// trim(case_file) //'" not found')
    end if

  end if

  ! globalize multilevel mesh options
  call ml_mesh_opt % Bcast(0, comm)

  ! globalize logging levels
  call XMPI_Bcast_LoggingLevels(0, comm)

  ! globalize remaining parameters
  call XMPI_Bcast(export_vtk , 0, comm)
  call XMPI_Bcast(test_domain, 0, comm)
  call XMPI_Bcast(gmsh_file  , 0, comm)
  call XMPI_Bcast(pg         , 0, comm)
  call XMPI_Bcast(po         , 0, comm)

  ! base mesh ..................................................................

  select case(test_domain)
  case(0)
    if (rank == 0) then
      call ImportGMSH_3D(gmsh_file, generic_mesh)
    end if
    call base_mesh % ImportGenericMesh(generic_mesh, comm)
    domain_name = gmsh_file
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
  end select

  if (rank == 0) then
    write(*,'(/,A)') 'verifying base mesh'
  end if

  call VerifyMesh_3D(base_mesh, passed)
  call XMPI_Reduce(passed, all_passed, MPI_LAND, 0, comm)

  if (rank == 0) then
    write(*,'(2X,A,G0)') 'passed = ', all_passed
  end if

  call MPI_Barrier(comm)

  ! multilevel mesh ............................................................

  ml_mesh = ML_Mesh_3D(base_mesh, ml_mesh_opt)
  l_top   = size(ml_mesh%mesh)

  associate(mesh => ml_mesh%mesh)

    ! verification
    do l = 1, l_top
      call VerifyMesh_3D(mesh(l), passed)
      call XMPI_Allreduce(passed, all_passed, MPI_LAND, comm)
      if (.not. all_passed) exit
    end do
    call MPI_Barrier(comm)
    if (rank == 0) then
      if (all_passed) then
        write(*,'(2X,9G0)') 'verification: all levels passed'
      else
        write(*,'(2X,9G0)') 'verification of level ',l,' failed'
      end if
    end if

    ! print info
    do l = 1, l_top
      if (mesh(l)%part >= 0) then
        call XMPI_Reduce(mesh(l)%n_elem, ne_min, MPI_MIN, 0, mesh(l)%comm_parts)
        call XMPI_Reduce(mesh(l)%n_elem, ne_max, MPI_MAX, 0, mesh(l)%comm_parts)
        call XMPI_Reduce(mesh(l)%n_elem, ne_tot, MPI_SUM, 0, mesh(l)%comm_parts)
      end if
      if (mesh(l)%part == 0) then
        write(*,'(2X,A,I4,A,I5,A,3(A,I6))') &
          'level ',l,': n_parts =',mesh(l)%n_parts,',  ', &
          'min/max/sum(n_elem) = ',ne_min,' / ',ne_max,' / ',ne_tot
      end if
    end do

  end associate

  !-----------------------------------------------------------------------------
  ! Finalization

  call MPI_Finalize()

  !=============================================================================

end program ML_Test_Elliptic
