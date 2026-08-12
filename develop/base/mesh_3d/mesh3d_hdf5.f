!------------------------------------------------------------------------------!
! This file is part of HiSPEET: High-order Spectral Element Techniques         !
!                                                                              !
! Copyright (C) 2026 by the HiSPEET authors and the                            !
! Chair of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany                 !
!                                                                              !
! HiSPEET is free software: you can redistribute it and/or modify              !
! it under the terms of the GNU General Public License as published by         !
! the Free Software Foundation, either version 3 of the License, or            !
! (at your option) any later version.                                          !
!                                                                              !
! HiSPEET is distributed in the hope that it will be useful,                   !
! but WITHOUT ANY WARRANTY; without even the implied warranty of               !
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                         !
! See the GNU General Public License for more details.                         !
!                                                                              !
! You should have received a copy of the GNU General Public License            !
! along with HiSPEET. If not, see <http://www.gnu.org/licenses/>.              !
!------------------------------------------------------------------------------!

!> summary:  Test of mesh HDF5 read/write
!> author:   Joerg Stiller
!> date:     2024/06/11
!===============================================================================

program Mesh3d_HDF5
  use Kind_Parameters
  use Constants
  use HDF5_Binding
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

  use Root_Mesh_Partitioning__3D

  use Export_VTK_Volume_Data__3D

  use QOI__Smiley__3D

  implicit none

  !-----------------------------------------------------------------------------
  ! Declarations

  ! input parameters ...........................................................

  ! case file
  character(len=*), parameter :: case_default = 'mesh3d_hdf5'
  character(len=80) :: case_file = case_default

  integer :: config = 1
  ! configuration (u/s = un/structured, r = regular, d = deformed)
  !   1  cuboidal domain with Cartesian mesh                               (s+r)
  !   2  cuboidal domain with unstructured "diamond" mesh                  (u+d)
  !   3  cuboidal domain with  3x3x3 elements and rotated center           (u+r)
  !   4  cylindrical domain                                                (u+d)
  !   5  annular domain                                                    (u+d)

  ! SE mesh and plotting
  integer :: n_parts = 0
  integer :: po = 3
  logical :: export_vtk = .true.

  namelist/control_prm/ config, n_parts, po, export_vtk

  ! MPI and OpenMP variables ...................................................

  type(MPI_Comm) :: comm     ! MPI communicator
  integer        :: rank     ! local MPI rank
  integer        :: n_proc   ! number of MPI processes
  integer        :: n_thread ! number of OpenMP threads

  ! mesh and variables .........................................................

  type(MeshPartitionerOptions_3D), save :: part_opt
  type(Mesh_3D),                   save :: old_mesh, mesh
  type(DataExchangePlan_3D),       save :: x_plan
  type(SpectralElementMesh_3D),    save :: sem

  ! HDF5 .......................................................................

  integer(hid_t) :: file_id, group_id

  ! auxiliary variables ........................................................

  type(QOI_Smiley_3D), save :: smiley

  real(RNP), allocatable, save :: s(:,:,:,:,:)

  character(len=80) :: config_name = ''
  character(len=80) :: input_file  = ''
  character(len=80) :: hdf5_file   = ''
  character(len=80) :: plot_file   = ''
  character(len=80) :: tag         = ''
  logical :: exists, passed, passed_loc
  integer :: err, io, stat
  integer :: e

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
    write(*,'(A)') 'Test of HDF5 read/write'
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
      close(io)
    else
       call Warning( 'Mesh3d_HDF5', 'input file "'//trim(input_file)// &
                     '" not found, using defaults' )
    end if

    if (n_parts < 1) then
      n_parts = n_proc
    else
      n_parts = max(1, min(n_proc, n_parts))
    end if

    part_opt % n_parts = n_parts

  end if

  ! globalize parameters
  call XMPI_Bcast( config        , 0, comm )
  call XMPI_Bcast( n_parts       , 0, comm )
  call XMPI_Bcast( po            , 0, comm )
  call XMPI_Bcast( export_vtk    , 0, comm )

  ! globalize partitioning parameters
  call part_opt % Bcast( 0, comm )

  ! mesh generation ............................................................

  select case(config)
  case(2)
    call CreateCuboidDiamonds(comm, input_file, old_mesh)
    config_name = 'cuboidal domain with unstructured "diamond" mesh'
  case(3)
    call CreateCuboidOneRotated(comm, input_file, old_mesh)
    config_name = 'cuboidal domain with 3x3x3 elements and rotated center'
  case(4)
    call CreateCylinder(comm, input_file, old_mesh)
    config_name = 'cylindrical domain with unstructured mesh'
  case(5)
    call CreateAnnulus(comm, input_file, old_mesh)
    config_name = 'annular domain with unstructured mesh'
  case default
    call CreateCuboidCartesian(comm, input_file, old_mesh)
    config_name = 'cuboidal domain with Cartesian mesh'
  end select

  if (rank == 0) then
    write(*,'(2X,A,/)') 'Created '//trim(config_name)
  end if

  ! root mesh partitioning .....................................................

  if (old_mesh % n_parts /= part_opt % n_parts) then
    call RootMeshPartitioning_3D( opt      = part_opt  &
                                , old_mesh = old_mesh  &
                                , new_mesh = mesh      &
                                , x_plan   = x_plan    )
  else
    mesh = old_mesh
  end if
  call VerifyMesh_3D(mesh, passed_loc)
  call XMPI_Reduce(passed_loc, passed, MPI_LAND, 0, comm)

  if (rank == 0) then
    if (passed) then
      write(*,'(2X,A,/)') 'Mesh successfully partitioned'
    else
      call Error('Mesh3d_HDF5', 'Mesh partitioning failed')
    end if
  end if

  !-----------------------------------------------------------------------------
  ! HDF5 test

  ! preliminaries ..............................................................

  ! make present mesh the original one
  old_mesh = mesh

  ! safeguard
  call Init_HDF5_Binding()

  ! each MPI rank gets own file name
  write(tag, fmt='(A2,I0)') '_r', rank
  hdf5_file = trim(case_file)//trim(tag)//'.h5'

  ! write mesh parts ...........................................................

  if (mesh%part >= 0) then

    ! create HDF5 file and group
    call H5Fcreate_f(trim(hdf5_file), H5F_ACC_TRUNC_F, file_id, err)
    call H5Gcreate_f(file_id, 'mesh', group_id, err)

    ! write mesh partition
    call old_mesh % WriteHDF5(group_id)

    ! close HDF5 group and file
    call H5Gclose_f(group_id, err)
    call H5Fclose_f(file_id, err)

  end if

  if (rank == 0) then
    if (err == 0) then
      write(*,'(2X,A,/)') 'Wrote mesh to HDF5'
    else
      call Error('Mesh3d_HDF5', 'Failed to write HDF5')
    end if
  end if

  ! read mesh parts ............................................................

  inquire(file=trim(hdf5_file), exist=exists)

  if (exists) then
    ! open HDF5 file and group for reading
    call H5Fopen_f(trim(hdf5_file), H5F_ACC_RDWR_F, file_id, err)
    call H5Gopen_f(file_id, 'mesh', group_id, err)
  else
    group_id = H5I_INVALID_HID_F
  end if

  ! read mesh partition
  call mesh % ReadHDF5(group_id, comm)

  if (exists) then
    ! close HDF5 group and file
    call H5Gclose_f(group_id, err)
    call H5Fclose_f(file_id, err)
  end if

  if (rank == 0) then
    if (err == 0) then
      write(*,'(2X,A,/)') 'Read mesh from HDF5'
    else
      call Error('Mesh3d_HDF5', 'Failed to read HDF')
    end if
  end if

  ! release HDF5 resources .....................................................

  call H5close_f(err)

  ! verify mesh ................................................................

  call VerifyMesh_3D(mesh, passed_loc)
  call XMPI_Reduce(passed_loc, passed, MPI_LAND, 0, comm)

  if (rank == 0) then
    if (passed) then
      write(*,'(2X,A,/)') 'Mesh successfully restored'
    else
      call Error('Mesh3d_HDF5', 'Restoration of the mesh failed')
    end if
  end if

  !-----------------------------------------------------------------------------
  ! Plotting

  sem = SpectralElementMesh_3D(mesh, po)

  if (export_vtk) then  ! .and. mesh%part >= 0 !??

    allocate(s(0:po,0:po,0:po,mesh%n_elem,3))

    do e = 1, mesh%n_elem
      s(:,:,:,e,1) = smiley % Density( sem%metrics%x(:,:,:,e,1) &
                                     , sem%metrics%x(:,:,:,e,2) &
                                     , sem%metrics%x(:,:,:,e,3) )
      s(:,:,:,e,2) = e
      s(:,:,:,e,3) = mesh % part
    end do

    plot_file = trim(case_file) !// trim(tag)

    call ExportVTK_VolumeData( sem % metrics % x        &
                             , s, sname = ['f','e','p'] &
                             , file    = plot_file      &
                             , part    = mesh % part    &
                             , n_parts = mesh % n_parts )
    deallocate(s)

  end if

  !-----------------------------------------------------------------------------
  ! Finalization

  call MPI_Finalize()

  !=============================================================================

end program Mesh3d_HDF5
