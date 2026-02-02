!> summary:  Test program for importing Gmsh mesh files
!> author:   Matthias Frey, Joerg Stiller
!> date:     2023/05/10
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

program Mesh3d_Import_GMSH
  use Kind_Parameters
  use Constants
  use XMPI
  use Import_GMSH__3D
  use Generic_Mesh__3D
  use Mesh__3D
  use Spectral_Element_Mesh__3D
  use Verify_Mesh__3D
  use Export_VTK_Volume_Data__3D
  use Export_VTK_Mesh_SFC__3D
  implicit none

  character(len=80) :: file = '../gmsh_3d/pipe' ! mesh file base name (*.msh)
  integer           :: po   = 1                 ! degree of spectral elements

  namelist/input/ file, po

  type(GenericMesh_3D)         :: generic_mesh
  type(Mesh_3D)                :: mesh
  type(SpectralElementMesh_3D) :: sem

  real(RNP), allocatable :: var(:,:,:,:,:)
  character(len=20), allocatable :: var_name(:)
  logical :: passed, all_passed
  integer :: rank, n_proc
  integer :: n_elem, n_var, prm
  integer :: e, i

  ! initialization .............................................................

  call XMPI_Init()
  call MPI_Comm_rank(MPI_COMM_WORLD, rank)
  call MPI_Comm_size(MPI_COMM_WORLD, n_proc)

  ! read parameters
  if (rank == 0) then
    write(*,'(A)') 'Testing import of GMSH files'
    write(*,'(2X,A)') 'reading control parameters ..'
    open(newunit = prm, file = 'mesh3d_import_gmsh.prm')
    read(prm, nml = input)
    close(prm)
  end if
  call XMPI_Bcast(file, 0, comm = MPI_COMM_WORLD)
  call XMPI_Bcast(po  , 0, comm = MPI_COMM_WORLD)

  ! mesh import ................................................................

  if (rank == 0) then
    call ImportGMSH_3D(file, generic_mesh)
  end if

  call mesh % ImportGenericMesh(generic_mesh, comm = MPI_COMM_WORLD)
  call VerifyMesh_3D(mesh, passed)

  call XMPI_Reduce(passed, all_passed, MPI_LAND, 0, MPI_COMM_WORLD)
  if (rank == 0) then
    write(*,'(/,2X,A,G0)') 'verifying imported mesh: passed = ', all_passed
  end if

  ! eval mesh quality ..........................................................

  if (rank == 0) then
    write(*,'(/,2X,A)') 'evaluating mesh quality ..'
  end if

  sem = SpectralElementMesh_3D(mesh, po)

  n_elem = mesh % n_elem

  n_var = 6
  allocate(var(0:po,0:po,0:po,n_elem,n_var), source = ZERO)
  allocate(var_name(n_var))
  var_name(1) = 'elem_id'
  var_name(2) = 'xi_1'
  var_name(3) = 'xi_2'
  var_name(4) = 'xi_3'
  var_name(5) = 'J'
  var_name(6) = 'Q'

  associate(xi => sem % std_op % x)
    do e = 1, n_elem
      var(:,:,:,e,1) = e
      do i = 0, po
        var(i,:,:,e,2) = xi(i)
        var(:,i,:,e,3) = xi(i)
        var(:,:,i,e,4) = xi(i)
      end do
      var(:,:,:,e,5) = sem % metrics % Jd(:,:,:,e)
      var(:,:,:,e,6) = minval(sem % metrics % Jd(:,:,:,e)) &
                     / maxval(sem % metrics % Jd(:,:,:,e))
    end do
  end associate

  ! export mesh ................................................................

  if (rank == 0) then
    write(*,'(/,2X,A)') 'exporting mesh to VTK ..'
  end if

  call ExportVTK_VolumeData( sem%metrics%x            &
                           , s       = var            &
                           , sname   = var_name       &
                           , file    = file           &
                           , part    = mesh % part    &
                           , n_parts = mesh % n_parts &
                           , subdiv  = .false.        )


  if (mesh % has_sfc) then
    call ExportVTK_MeshSFC(mesh, file = trim(file)//'_sfc')
  end if

  ! finalization ...............................................................

  call MPI_Finalize()

  !=============================================================================

end program Mesh3d_Import_GMSH
