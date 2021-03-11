program Mesh3d_Import_Generic
  use Kind_Parameters
  use Constants
  use XMPI
  use Generic_Mesh_3d
  use Mesh_3d__Partition
  use Verify_Mesh_3d
  use Export_VTK_3d__Volume_Data
  implicit none

  type(GenericMesh3d) :: generic_mesh
  type(Mesh3d_Partition) :: mesh
  real(RNP), allocatable :: x(:,:,:,:,:)
  real(RNP), allocatable :: s(:,:,:,:,:)

  character(len=*), parameter :: input_file = 'mesh3d_import_generic.prm'
  integer   :: conf = 1   ! configuration (1 cylinder, 2 annular gap)
  real(RNP) :: r0   = 0.5 ! inner radius  (annular gap only)
  real(RNP) :: r1   = 1   ! outer radius
  real(RNP) :: h    = 2   ! height
  integer   :: nr   = 2   ! num elements in radial    direction
  integer   :: np   = 4   ! num elements in azimuthal direction (gap only)
  integer   :: nz   = 2   ! num elements in axial     direction
  integer   :: po   = 3   ! polynomial order of mesh elements
  logical   :: periodic = .true. ! switch for axial periodicity

  namelist/input/ conf, r0, r1, h, nr, np, nz, po, periodic

  type(MPI_Comm) :: comm = MPI_COMM_WORLD
  integer :: rank
  integer :: io     ! IO unit number
  integer :: ns = 1 ! number of scalar variables
  logical :: passed

  call MPI_Init()
  call MPI_Comm_rank(comm, rank)

  if (rank == 0) then

    ! read parameters
    open(newunit = io, file = input_file)
    read(io, nml = input)
    close(io)

    select case(conf)
    case(2)
      call generic_mesh % CreateAnnularGap(r0, r1, h, nr, np, nz, po, periodic)
    case default
      call generic_mesh % CreateCylinder(nr, nz, po, periodic)
    end select

    call mesh % ImportGenericMesh(generic_mesh, comm = comm)

    call VerifyMesh3d(mesh, passed)
    write(*,'(/,A,G0,/)') 'VerifyMesh3d: passed = ', passed

    call mesh % GetPoints(po, 'GLL', x)

    allocate(s(0:po, 0:po, 0:po, mesh%n_elem, ns), source = ZERO)

    call ExportVTK_VolumeData( x, s                    &
                             , sname  = ['s']          &
                             , file   = 'generic_mesh' &
                             , part   = mesh % part    &
                             , n_part = mesh % n_part  )

  end if

  call MPI_Finalize()

end program Mesh3d_Import_Generic
