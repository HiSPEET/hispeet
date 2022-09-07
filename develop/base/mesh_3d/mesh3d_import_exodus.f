program Mesh3d_Import_Exodus
  use Kind_Parameters
  use Constants
  use XMPI
  use Generic_Mesh__3D
  use Import_Exodus__3D
  use Mesh__3D
  use Verify_Mesh__3D
  use Export_VTK_Volume_Data__3D
  implicit none

  character(len=80) :: exodus_file = 'exodus/naca_0012.e'

  type(GenericMesh_3D) :: generic_mesh
  type(Mesh_3D) :: mesh
  real(RNP), allocatable :: x(:,:,:,:,:)
  real(RNP), allocatable :: s(:,:,:,:,:)
  integer :: po = 1, ne, ns = 1
  logical :: passed

  call MPI_Init()
  call ImportExodus3d(exodus_file, generic_mesh)
  call generic_mesh % SwitchToLexicalNumbering()
  call mesh % ImportGenericMesh(generic_mesh, comm = MPI_COMM_WORLD)
  call VerifyMesh_3D(mesh, passed)
  write(*,'(A,G0,/)') 'VerifyMesh_3D: passed = ', passed
  call mesh % GetPoints(po, 'L', x)
  ne = mesh % n_elem
  allocate(s(0:po,0:po,0:po,ne,ns), source = ZERO)
  call ExportVTK_VolumeData( x, s                     &
                           , sname   = ['s']          &
                           , file    = 'naca_0012'    &
                           , part    = mesh % part    &
                           , n_parts = mesh % n_parts )
  call MPI_Finalize()

end program Mesh3d_Import_Exodus
