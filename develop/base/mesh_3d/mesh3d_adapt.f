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
  use Spectral_Element_Mesh__3D
  use Root_Mesh_Partitioning__3D
  use Child_Mesh_Adaptation__3D
  use Export_VTK_Volume_Data__3D

  use Smiling_Face

  implicit none

  !-----------------------------------------------------------------------------
  ! Declarations

  ! input parameters ...........................................................

  ! input file (*.prm)
  character(len=*), parameter :: input_default = 'mesh3d_adapt'
  character(len=80) :: input_file = ''

  integer :: config = 1
  ! configuration (u/s = un/structured, r = regular, d = deformed)
  !   1  cuboidal domain with Cartesian mesh                               (s+r)
  !   2  cuboidal domain with unstructured "diamond" mesh                  (u+d)
  !   3  cuboidal domain with  3x3x3 elements and rotated center           (u+r)
  !   4  cylindrical domain                                                (u+d)
  !   5  annular domain                                                    (u+d)

  ! SE mesh and plotting
  logical :: export_vtk = .true.
  integer :: po = 3

  namelist/control_prm/ config, export_vtk, po

  ! adaptation
  integer :: n_level        = 1
  integer :: n_parts_base   = 1
  integer :: n_parts_growth = 1

  namelist/adaptation_prm/ n_level, n_parts_base, n_parts_growth

  ! MPI and OpenMP variables ...................................................

  type(MPI_Comm) :: comm     ! MPI communicator
  integer        :: rank     ! local MPI rank
  integer        :: n_proc   ! number of MPI processes
  integer        :: n_thread ! number of OpenMP threads

  ! mesh and variables .........................................................

  type(PartitioningOptions_3D), allocatable, save :: part_opt(:)
  type(Mesh_3D),                allocatable, save :: orig_mesh(:), mesh(:)
  type(SpectralElementMesh_3D), allocatable, save :: sem(:)

  ! auxiliary variables ........................................................

  type(SmilingFace), save :: smiley

  real(RNP), allocatable, save :: s(:,:,:,:,:)
  logical,   allocatable, save :: mask(:)

  character(len=80) :: config_name = ''
  character(len=80) :: plot_file   = ''
  character(len=80) :: tag         = ''
  logical :: exists
  integer :: io, stat
  integer :: l

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
      open(newunit = io, file = input_file)
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
    do l = 2, n_level
      part_opt(l) % mode    = 2 ! switch child partitioning
      part_opt(l) % n_parts = min(n_parts_growth * part_opt(l-1)%n_parts, n_proc)
    end do

  end if

  ! globalize parameters
  call XMPI_Bcast( config        , 0, comm )
  call XMPI_Bcast( export_vtk    , 0, comm )
  call XMPI_Bcast( po            , 0, comm )
  call XMPI_Bcast( n_level       , 0, comm )
  call XMPI_Bcast( n_parts_base  , 0, comm )
  call XMPI_Bcast( n_parts_growth, 0, comm )

  if (rank > 0) then
    allocate(part_opt(n_level))
  end if

  ! globalize partitioning parameters
  do l = 1, n_level
    call part_opt(l) % Bcast( 0, comm )
  end do

  allocate(orig_mesh(n_level))
  allocate(mesh(n_level))
  allocate(sem(n_level))

  ! mesh generation ............................................................

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

  if (orig_mesh(1) % n_parts /= part_opt(1) % n_parts) then
    call RootMeshPartitioning_3D(part_opt(1), orig_mesh(1), mesh(1))
  else
    mesh(1) = orig_mesh(1)
  end if

  sem(1) = SpectralElementMesh_3D(mesh(1), po)

  !-----------------------------------------------------------------------------
  ! Adaptation

!### CHECK
print '(99(G0,1X))', '# 0'
!### CHECK END
  do l = 2, n_level

    mesh(l-1) % element % adaptation % mark = 0
    mesh(l-1) % element(1) % adaptation % mark = 100
    call ChildMeshAdaptation_3D(part_opt(l), mesh(l-1), mesh(l))
!### CHECK
print '(99(G0,1X))', '# 1a, l =',l
!### CHECK END
    sem(l) = SpectralElementMesh_3D(mesh(l), po)
!### CHECK
print '(99(G0,1X))', '# 1b, l =',l
!### CHECK END

  end do
!### CHECK
print '(99(G0,1X))', '# 1'
!### CHECK END

  !-----------------------------------------------------------------------------
  ! Plotting

  if (export_vtk) then
    do l = 1, n_level

allocate(s(0:po,0:po,0:po,mesh(l)%n_elem,1))
s(:,:,:,:,1) = smiley % Density(sem(l)%metrics%x(:,:,:,:,1), sem(l)%metrics%x(:,:,:,:,2))

      write(tag, fmt='(A2,I0)') '_l', l
      plot_file = trim(input_file) // trim(tag)

      allocate(mask(mesh(l)%n_elem))
      if (l < n_level) then
        mask = mesh(l) % element % adaptation % mark == 0
      else
        mask = .true.
      end if

      call ExportVTK_VolumeData( sem(l) % metrics % x        &
                               , s, sname = ['f']            &
                               , file    = plot_file         &
                               , part    = mesh(l) % part    &
                               , n_parts = mesh(l) % n_parts &
                               , mask    = mask              )
      deallocate(mask)
deallocate(s)
!### CHECK
print '(99(G0,1X))', '# 1, l =', l
!### CHECK END
    end do
  end if
!### CHECK
print '(99(G0,1X))', '# X'
!### CHECK END

  !-----------------------------------------------------------------------------
  ! Finalization

  call MPI_Finalize()

  !=============================================================================

end program Mesh3d_Adapt
