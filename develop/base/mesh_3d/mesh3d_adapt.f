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
  use Verify_Mesh__3D

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
  character(len=80) :: case_file = case_default

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
  integer :: n_level        = 1 ! maximum top level L ≤ L_max ≡ n_level
  integer :: n_parts_base   = 1
  integer :: n_parts_growth = 1

  namelist/adaptation_prm/ scenario, n_level, n_parts_base, n_parts_growth

  ! adaptation criterion: elements on level l will be refined if the quantity
  ! of interest equals or exceeds c(l) = c₁ + ∆c⋅(l-1)/(L_max-1) in any point
  real(RNP) :: c1 = 0  ! c₁
  real(RNP) :: dc = 0  ! ∆c
  real(RNP), allocatable :: c(:)

  namelist/adaptation_prm/ c1, dc

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

  real(RNP), allocatable :: s_e(:,:,:)
  character(len=80) :: config_name = ''
  character(len=80) :: input_file  = ''
  character(len=80) :: plot_file   = ''
  character(len=80) :: tag         = ''
  logical :: exists, passed, passed_loc
  integer :: io, stat
  integer :: e, l, m

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
  call XMPI_Bcast( c1            , 0, comm )
  call XMPI_Bcast( dc            , 0, comm )

  if (rank > 0) then
    allocate(part_opt(n_level))
  end if

  ! globalize partitioning parameters
  do l = 1, n_level
    call part_opt(l) % Bcast( 0, comm )
  end do

  ! adaptation criterion
  allocate(c(n_level-1))
  do l = 1, n_level-1
    c(l) = c1 + dc * (l - 1) / (n_level - 1)
  end do

  allocate(s_e(0:po,0:po,0:po))

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
  call VerifyMesh_3D(mesh(1), passed)

  if (passed) then
    write(*,'(2X,A,/)') 'Root mesh successfully partitioned'
  end if

  sem(1) = SpectralElementMesh_3D(mesh(1), po)

  !-----------------------------------------------------------------------------
  ! Adaptation

  do m = 1, n_level-1

    if (rank == 0) then
      write(*,'(2X,99G0,/)') 'Adaptation cycle ', m
    end if

    ! set adaptation marks
    do l = 1, m
      associate(x => sem(l) % metrics % x)
        do e = 1, mesh(l) % n_elem
          associate(element => mesh(l) % element(e))
            if (element % frozen) then
              call element % MarkForRemoval()
            else
              s_e = smiley % Density(x = x(:,:,:,e,1), y = x(:,:,:,e,2))
              if (any(s_e >= c(l))) then
                call element % MarkForRefinement()
              else
                call element % MarkForRemoval()
              end if
            end if
          end associate
        end do
      end associate
    end do

    ! make adaptation pattern consistent
    do l = m, 2, -1
      call GlobalizeAdaptationPattern_3D(mesh(l))
      call RestrictAdaptationPattern_3D(mesh(l), mesh(l-1))
    end do
    call GlobalizeAdaptationPattern_3D(mesh(1))

    ! make present mesh the original one
    call move_alloc(mesh, old_mesh)
    allocate(mesh(n_level))

    ! create data exchange plan for redistribution of retained data
    allocate(exch_plan(n_level))

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
    sem(1) = SpectralElementMesh_3D(mesh(1), po)

    do l = 1, m
      if (l > 1) then
        call ProcessAdaptationPattern_3D(mesh(l))
      end if
      select case(m-l)
      case(0)
        call ChildMeshAdaptation_3D( opt        = part_opt(l+1) &
                                   , parent     = mesh(l)       &
                                   , new_child  = mesh(l+1)     )
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
     sem(l+1) = SpectralElementMesh_3D(mesh(l+1), po)

    end do


    if (rank == 0) then
      write(*,'(4X,99G0)') 'verify mesh'
    end if
    do l = 1, min(m+1, n_level)
      call VerifyMesh_3D(mesh(l), passed_loc)
      call XMPI_Reduce(passed_loc, passed, MPI_LAND, 0, comm)
      if (rank == 0) then
        write(*,'(6X,A,I3,2X,A,I7,2X,A,2X,L1)') &
          'level',l,', n_elem =',mesh(l)%n_elem,', passed:',passed
      end if
    end do

    deallocate(exch_plan)

  end do

  !-----------------------------------------------------------------------------
  ! Plotting

  if (export_vtk) then
    do l = 1, n_level

      allocate(s(0:po,0:po,0:po,mesh(l)%n_elem,3))
      s(:,:,:,:,1) = smiley % Density( sem(l)%metrics%x(:,:,:,:,1) &
                                     , sem(l)%metrics%x(:,:,:,:,2) )
      do e = 1, mesh(l)%n_elem
        s(:,:,:,e,2) = e
        s(:,:,:,e,3) = mesh(l) % part
      end do

      write(tag, fmt='(A2,I0)') '_l', l
      plot_file = trim(case_file) // trim(tag)

      allocate(mask(mesh(l)%n_elem))
      if (l < n_level) then
        mask = mesh(l) % element % adaptation % refinement < 100 &
               .and. .not. mesh(l) % element % frozen
      else
        mask = .not. mesh(l) % element % frozen
      end if

      call ExportVTK_VolumeData( sem(l) % metrics % x        &
                               , s, sname = ['f','e','p']    &
                               , file    = plot_file         &
                               , part    = mesh(l) % part    &
                               , n_parts = mesh(l) % n_parts &
                               , mask    = mask              )
      deallocate(mask)
      deallocate(s)
    end do
  end if

  !-----------------------------------------------------------------------------
  ! Finalization

  call MPI_Finalize()

  !=============================================================================

end program Mesh3d_Adapt
