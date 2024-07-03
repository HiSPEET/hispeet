!> summary:  Program for testing basic multilevel functionality
!> author:   Joerg Stiller
!> date:     2024/07/02
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

program ML_Test_Functionality
! use Kind_Parameters
  use XMPI
  use Import_GMSH__3D
  use Generic_Mesh__3D
  use Mesh__3D
  use Verify_Mesh__3D
  use ML__Mesh__3D
  use ML__Mesh_Operators__3D
  use ML__Mesh_Variable__3D
  implicit none

  ! variables ..................................................................

  ! base name of the mesh file (gmsh/*.msh)
  character(len=80) :: file = 'cylinder_2d'

  ! spectral element mesh options
  integer, allocatable :: po(:)  ! sequence of polynomial orders

  namelist /control/   file
  namelist /operators/ po

  type(GenericMesh_3D)     , save :: generic_mesh
  type(Mesh_3D)            , save :: base_mesh
  type(ML_Mesh_Options_3D) , save :: ml_mesh_opt
  type(ML_Mesh_3D)         , save :: ml_mesh
  type(ML_MeshOperators_3D), save :: ml_op

  type(MPI_Comm) :: comm = MPI_COMM_WORLD
  character(len=:), allocatable :: gmsh_file
  logical :: passed, all_passed
  integer :: rank, n_proc
  integer :: prm

  ! initialization .............................................................

  call XMPI_Init()
  call MPI_Comm_rank(comm, rank)
  call MPI_Comm_size(comm, n_proc)

  ! read parameters
  if (rank == 0) then
    write(*,'(A)') 'Testing basic multilevel functionality'
    write(*,'(2X,A)') 'reading input parameters ..'
    open(newunit = prm, file = 'ml_test_functionality.prm')
    read(prm, nml = control)
    ml_mesh_opt = ML_Mesh_Options_3D(prm, n_proc)
    allocate(po(ml_mesh_opt%l_top), source = -1)
    read(prm, nml = operators)
    close(prm)
  end if

  ! globalize parameters
  call ml_mesh_opt % Bcast(0, comm)
  if (rank > 0) then
    allocate(po(ml_mesh_opt%l_top), source = -1)
  end if
  call XMPI_Bcast(po, 0, comm)

  ! mesh import ................................................................

  gmsh_file = 'gmsh/' // trim(file)

  if (rank == 0) then
    call ImportGMSH_3D(gmsh_file, generic_mesh)
  end if

  call base_mesh % ImportGenericMesh(generic_mesh, comm = MPI_COMM_WORLD)
  call VerifyMesh_3D(base_mesh, passed)

  call XMPI_Reduce(passed, all_passed, MPI_LAND, 0, MPI_COMM_WORLD)
  if (rank == 0) then
    write(*,'(/,2X,A,G0)') 'verifying imported mesh: passed = ', all_passed
  end if

  ! multilevel mesh ............................................................

  ml_mesh = ML_Mesh_3D(base_mesh, ml_mesh_opt)

  ! finalization ...............................................................

  call MPI_Finalize()

  !=============================================================================

end program ML_Test_Functionality
