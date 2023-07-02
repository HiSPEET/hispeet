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
  use Data_Exchange__3D
  use Spectral_Element_Mesh__3D

  use Child_Mesh_Adaptation__3D
  use Globalize_Adaptation_Pattern__3D
  use Process_Adaptation_Pattern__3D
  use Restrict_Adaptation_Pattern__3D
  use Root_Mesh_Partitioning__3D

  use Export_VTK_Volume_Data__3D

  use Smiling_Face

  implicit none

  !-----------------------------------------------------------------------------
  ! Declarations

  ! input parameters ...........................................................

  ! case file
  character(len=*), parameter :: case_default = 'mesh3d_adapt'
  character(len=80) :: case_file = 'mesh3d_adapt'

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
  integer :: scenario       = 1 ! 1 global
  integer :: n_level        = 1
  integer :: n_parts_base   = 1
  integer :: n_parts_growth = 1

  namelist/adaptation_prm/ scenario, n_level, n_parts_base, n_parts_growth

  ! MPI and OpenMP variables ...................................................

  type(MPI_Comm) :: comm     ! MPI communicator
  integer        :: rank     ! local MPI rank
  integer        :: n_proc   ! number of MPI processes
  integer        :: n_thread ! number of OpenMP threads

  ! mesh and variables .........................................................

  type(PartitioningOptions_3D), allocatable, save :: part_opt(:)
  type(Mesh_3D),                allocatable, save :: old_mesh(:), mesh(:)
  type(DataExchangePlan_3D),    allocatable, save :: exch_plan(:)
  type(SpectralElementMesh_3D), allocatable, save :: sem(:)

  ! auxiliary variables ........................................................

  type(SmilingFace), save :: smiley

  real(RNP), allocatable, save :: s(:,:,:,:,:)
  logical,   allocatable, save :: mask(:)

  character(len=80) :: config_name = ''
  character(len=80) :: input_file  = ''
  character(len=80) :: plot_file   = ''
  character(len=80) :: tag         = ''
  logical :: exists
  integer :: io, stat
  integer :: l, m

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

    call get_command_argument(1, case_file, status=stat)
    if (stat /= 0 .or. len_trim(case_file) == 0) then
      case_file = case_default
    end if
    input_file = trim(case_file) // '.prm'

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
  call XMPI_Bcast( scenario      , 0, comm )
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

  allocate(old_mesh(n_level))
  allocate(mesh(n_level))
  allocate(sem(n_level))

  ! mesh generation ............................................................

  select case(config)
  case(2)
    call CreateCuboidDiamonds(comm, input_file, old_mesh(1))
    config_name = 'Cuboidal domain with unstructured "diamond" mesh'
  case(3)
    call CreateCuboidOneRotated(comm, input_file, old_mesh(1))
    config_name = 'Cuboidal domain with 3x3x3 elements and rotated center'
  case(4)
    call CreateCylinder(comm, input_file, old_mesh(1))
    config_name = 'Cylindrical domain with unstructured mesh'
  case(5)
    call CreateAnnulus(comm, input_file, old_mesh(1))
    config_name = 'Annular domain with unstructured mesh'
  case default
    call CreateCuboidCartesian(comm, input_file, old_mesh(1))
    config_name = 'Cuboidal domain with Cartesian mesh'
  end select

  ! root mesh partitioning .....................................................

  if (old_mesh(1) % n_parts /= part_opt(1) % n_parts) then
    call RootMeshPartitioning_3D(part_opt(1), old_mesh(1), mesh(1))
  else
    mesh(1) = old_mesh(1)
  end if

  sem(1) = SpectralElementMesh_3D(mesh(1), po)

  !-----------------------------------------------------------------------------
  ! Adaptation

  do m = 1, n_level-1
!### CHECK
write(*,'(99(G0,1X))') 'starting cycle m =', m
!### CHECK END

    ! set adaptation marks
    mesh(m) % element % adaptation % mark = 1

!### CHECK
write(*,'(99(G0,1X))') 'make adaptation pattern consistent'
!### CHECK END
    ! make adaptation pattern consistent
    do l = m, 2, -1
!### CHECK
write(*,'(99(G0,1X))') '... level l =',l
!### CHECK END
      call GlobalizeAdaptationPattern_3D(mesh(l))
!### CHECK
write(*,'(99(G0,1X))') '...... globalized'
!### CHECK END
      call RestrictAdaptationPattern_3D(mesh(l), mesh(l-1))
!### CHECK
write(*,'(99(G0,1X))') '...... restricted'
!### CHECK END
    end do
!### CHECK
write(*,'(99(G0,1X))') '... level l =',1
!### CHECK END
    call GlobalizeAdaptationPattern_3D(mesh(1))

!### CHECK
write(*,'(99(G0,1X))') 'make present mesh the original one'
!### CHECK END
    ! make present mesh the original one
    call move_alloc(mesh, old_mesh)
    allocate(mesh(n_level))

    ! create data exchange plan for redistribution of retained data
    allocate(exch_plan(n_level))

!### CHECK
write(*,'(99(G0,1X))') 'root level'
!### CHECK END
    ! root level
    call ProcessAdaptationPattern_3D(old_mesh(1))
    if (max(old_mesh(1) % n_parts, part_opt(1) % n_parts) == 1) then
      mesh(1) = old_mesh(1)
    else if (old_mesh(1) % is_top) then
      call RootMeshPartitioning_3D( opt       = part_opt(1)  &
                                  , old_mesh  = old_mesh(1)  &
                                  , new_mesh  = mesh(1)      &
                                  , exch_plan = exch_plan(1) )
    else
      call RootMeshPartitioning_3D( opt       = part_opt(1)  &
                                  , old_mesh  = old_mesh(1)  &
                                  , new_mesh  = mesh(1)      &
                                  , child     = old_mesh(2)  &
                                  , exch_plan = exch_plan(1) )
    end if

    do l = 1, m
!### CHECK
write(*,'(99(G0,1X))') 'refining level l =', l
!### CHECK END
      if (l > 1) then
        call ProcessAdaptationPattern_3D(mesh(l))
      end if
!### CHECK
write(*,'(99(G0,1X))') 'creating level', l+1
write(*,'(99(G0,1X))') 'mesh(',l,')%mark =', mesh(l)%element%adaptation%mark
write(*,'(99(G0,1X))') 'mesh(',l,')%subl =', mesh(l)%element%adaptation%sublevels
!### CHECK END
      select case(m-l)
      case(0)
        call ChildMeshAdaptation_3D( opt        = part_opt(l+1) &
                                   , parent     = mesh(l)       &
                                   , new_child  = mesh(l+1)     )
!### CHECK
write(*,'(99(G0,1X))') 'mesh(',l,')%n_child =',mesh(l)%n_child
write(*,'(99(G0,1X))') 'mesh(',l+1,')%n_elem =',mesh(l+1)%n_elem
!### CHECK END
      case(1)
        call ChildMeshAdaptation_3D( opt       = part_opt(l+1)  &
                                   , parent    = mesh(l)        &
                                   , new_child = mesh(l+1)      &
                                   , old_child = old_mesh(l+1)  &
                                   , exch_plan = exch_plan(l+1) )
      case(2:)
        call ChildMeshAdaptation_3D( opt        = part_opt(l+1)  &
                                   , parent     = mesh(l)        &
                                   , new_child  = mesh(l+1)      &
                                   , old_child  = old_mesh(l+1)  &
                                   , grandchild = old_mesh(l+2)  &
                                   , exch_plan  = exch_plan(l+1) )
      end select
!### CHECK
write(*,'(99(G0,1X))') 'finished child mesh adaptation'
!### CHECK END
     sem(l+1) = SpectralElementMesh_3D(mesh(l+1), po)

    end do

    deallocate(exch_plan)

  end do

  !-----------------------------------------------------------------------------
  ! Plotting

  if (export_vtk) then
    do l = 1, n_level

allocate(s(0:po,0:po,0:po,mesh(l)%n_elem,1))
s(:,:,:,:,1) = smiley % Density(sem(l)%metrics%x(:,:,:,:,1), sem(l)%metrics%x(:,:,:,:,2))

      write(tag, fmt='(A2,I0)') '_l', l
      plot_file = trim(case_file) // trim(tag)

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
