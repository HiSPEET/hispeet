!> summary:  Test of mesh adaptation
!> author:   Joerg Stiller
!> date:     2023/03/26
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> If present, the first argument of the invoking command will be interpreted
!> as the base name of the control file. If omitted, the program looks for
!> `mesh3d_adapt.prm`.
!===============================================================================

program Mesh3d_Adapt
  use Kind_Parameters
  use Constants
  use OpenMP_Binding
  use Execution_Control
  use XMPI

  use Create_Cuboid_Cartesian
  use Create_Cuboid_Diamonds
  use Create_Cuboid_OneRotated
  use Create_Cylinder
  use Create_Annulus

  use Mesh__3D
  use Root_Mesh_Partitioning__3D

  implicit none

  !-----------------------------------------------------------------------------
  ! Declarations

  ! control parameters .........................................................

  ! input file (*.prm)
  character(len=*), parameter :: input_default = 'dg_elliptic_3d_test'
  character(len=80) :: input_file = ''

  integer :: config = 1
  ! configuration (u/s = un/structured, r = regular, d = deformed)
  !   1  cuboidal domain with Cartesian mesh                               (s+r)
  !   2  cuboidal domain with unstructured "diamond" mesh                  (u+d)
  !   3  cuboidal domain with  3x3x3 elements and rotated center           (u+r)
  !   4  cylindrical domain                                                (u+d)
  !   5  annular domain                                                    (u+d)

  ! plotting
  character(len=80) :: plot_file  = ''
  logical :: plot_subdiv = .true.

  namelist/control_prm/ config, plot_file, plot_subdiv

  ! MPI and OpenMP variables ...................................................

  type(MPI_Comm) :: comm     ! MPI communicator
  integer        :: rank     ! local MPI rank
  integer        :: n_proc   ! number of MPI processes
  integer        :: n_thread ! number of OpenMP threads

  ! mesh and variables .........................................................

  integer :: n_level       = 1
  integer :: n_parts_base   = 1
  integer :: n_parts_growth = 1

  type(PartitioningOptions_3D), allocatable, save :: part_opt(:)
  type(Mesh_3D), allocatable, save :: orig_mesh(:), mesh(:)

  namelist/adaptation_prm/ n_level, n_parts_base, n_parts_growth

  ! auxiliary variables ........................................................

  character(len=80) :: config_name = ''
  logical :: exists
  integer :: io, stat
  integer :: i

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
    write(*,'(A)') 'Test of mesh adaptation'
    write(*,*)

    call get_command_argument(1, input_file, status=stat)
    if (stat /= 0 .or. len_trim(input_file) == 0) then
      input_file = input_default
    end if
    input_file = trim(input_file) // '.prm'

    inquire(file=trim(input_file), exist=exists)
    if (exists) then
      write(*,'(2X,A)') 'reading ' // trim(input_file)
      read(io, nml = control_prm)
      read(io, nml = adaptation_prm)
      close(io)
    else
       call Warning( 'Mesh3d_Adapt', 'input file "'//trim(input_file)// &
                     '" not found, using defaults' )
    end if

    n_level        = max(1, n_level)
    n_parts_base   = max(1, min(n_proc, n_parts_base))
    n_parts_growth = max(1, n_parts_growth)

    allocate(part_opt(n_level))

    part_opt(1) % n_parts = n_parts_base
    do i = 2, n_level
      part_opt(i)%n_parts = min(n_parts_growth * part_opt(i-1)%n_parts, n_proc)
    end do

  end if

  ! globalize parameters
  call XMPI_Bcast( plot_file     , 0, comm )
  call XMPI_Bcast( plot_subdiv   , 0, comm )
  call XMPI_Bcast( config        , 0, comm )
  call XMPI_Bcast( n_level       , 0, comm )
  call XMPI_Bcast( n_parts_base  , 0, comm )
  call XMPI_Bcast( n_parts_growth, 0, comm )

  if (rank > 0) then
    allocate(part_opt(n_level))
  end if

  ! globalize partitioning parameters
  do i = 1, n_level
    call part_opt(i) % Bcast( 0, comm )
  end do

  ! mesh generation ............................................................

  allocate(orig_mesh(n_level))

  select case(config)
  case(2)
    call CreateCuboidDiamonds(comm, input_file, orig_mesh(1))
    config_name = 'Cuboidal domain with unstructured "diamond" mesh'
  case(3)
    call CreateCuboidOneRotated(comm, input_file, orig_mesh(1))
    config_name = 'Cuboidal domain with 3x3x3 elements and rotated center'
  case(4)
    call CreateCylinder(comm, input_file, orig_mesh(1))
    config_name = 'Cylindrical domain with unstructured mesh'
  case(5)
    call CreateAnnulus(comm, input_file, orig_mesh(1))
    config_name = 'Annular domain with unstructured mesh'
  case default
    call CreateCuboidCartesian(comm, input_file, orig_mesh(1))
    config_name = 'Cuboidal domain with Cartesian mesh'
  end select

  ! root mesh partitioning .....................................................

  allocate(mesh(n_level))

  if (orig_mesh(1) % n_parts /= part_opt(1) % n_parts) then
    call RootMeshPartitioning_3D(part_opt(1), orig_mesh(1), mesh(1))
  else
    mesh(1) = orig_mesh(1)
  end if

  !-----------------------------------------------------------------------------
  ! Finalization

  call MPI_Finalize()

  !=============================================================================

end program Mesh3d_Adapt
