!> summary:  Program for testing basic multilevel functionality
!> author:   Joerg Stiller
!> date:     2024/07/02
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

program ML_Test_Functionality
! use Kind_Parameters
  use Logging_Levels
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

  namelist /control/   log_level, file
  namelist /operators/ po

  type(GenericMesh_3D)     , save :: generic_mesh
  type(Mesh_3D)            , save :: base_mesh
  type(ML_Mesh_Options_3D) , save :: ml_mesh_opt
  type(ML_Mesh_3D)         , save :: ml_mesh
  type(ML_MeshOperators_3D), save :: ml_op

  type(MPI_Comm) :: comm = MPI_COMM_WORLD
  character(len=:), allocatable :: gmsh_file
  logical :: passed, all_passed
  integer :: rank, n_proc, prm
  integer :: ne_max, ne_min, ne_tot

  integer :: l

  ! initialization .............................................................

  call XMPI_Init()
  call MPI_Comm_rank(comm, rank)
  call MPI_Comm_size(comm, n_proc)

  ! read parameters
  if (rank == 0) then
    write(*,'(/,A,/)') 'Testing basic multilevel functionality'
    write(*,'(2X,A)') 'reading input parameters'
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

  call base_mesh % ImportGenericMesh(generic_mesh, comm)

  if (rank == 0) then
    write(*,'(/,A)') 'verifying imported mesh'
  end if

  call VerifyMesh_3D(base_mesh, passed)
  call XMPI_Reduce(passed, all_passed, MPI_LAND, 0, comm)

  if (rank == 0) then
    write(*,'(2X,A,G0)') 'passed = ', all_passed
  end if

  ! multilevel mesh ............................................................

  ml_mesh = ML_Mesh_3D(base_mesh, ml_mesh_opt)

  associate(mesh => ml_mesh%mesh)

    ! verification
    do l = 1, size(mesh)
      call VerifyMesh_3D(mesh(l), passed)
      call XMPI_Allreduce(passed, all_passed, MPI_LAND, comm)
      if (.not. all_passed) exit
    end do
    if (rank == 0) then
      if (all_passed) then
        write(*,'(2X,9G0)') 'verification: all levels passed'
      else
        write(*,'(2X,9G0)') 'verification of level ',l,' failed'
      end if
    end if

    ! print info
    do l = 1, size(mesh)
      if (mesh(l)%part >= 0) then
        call XMPI_Reduce(mesh(l)%n_elem, ne_min, MPI_MIN, 0, mesh(l)%comm_parts)
        call XMPI_Reduce(mesh(l)%n_elem, ne_max, MPI_MAX, 0, mesh(l)%comm_parts)
        call XMPI_Reduce(mesh(l)%n_elem, ne_tot, MPI_SUM, 0, mesh(l)%comm_parts)
      end if
      if (mesh(l)%part == 0) then
        write(*,'(2X,9G0)') &
          'level ',l,': min/max/sum(n_elem) = ',ne_min,' / ',ne_max,' / ',ne_tot
      end if
    end do

  end associate

  ! multilevel operators .......................................................

  ml_op = ML_MeshOperators_3D(ml_mesh, po)

  ! finalization ...............................................................

  call MPI_Finalize()

  !=============================================================================

end program ML_Test_Functionality
