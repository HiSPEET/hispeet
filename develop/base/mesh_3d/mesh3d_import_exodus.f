program Mesh3d_Import_Exodus
  use Kind_Parameters
  use Constants
  use XMPI
  use Generic_Mesh_3d
  use Import_Exodus_3d
  use Mesh_3d__Partition
  use Verify_Mesh_3d
  use Export_VTK_3d__Volume_Data
  implicit none

  character(len=80) :: exodus_file = 'exodus/naca_0012.e'

  type(GenericMesh3d) :: generic_mesh
  type(Mesh3d_Partition) :: mesh
  real(RNP), allocatable :: x(:,:,:,:,:)
  real(RNP), allocatable :: s(:,:,:,:,:)
  integer :: po = 1, ne, ns = 1
  logical :: passed

  call MPI_Init()
  call ImportExodus3d(exodus_file, generic_mesh)
  call generic_mesh % SwitchToLexicalNumbering()
  call mesh % ImportGenericMesh(generic_mesh, comm = MPI_COMM_WORLD)
  call VerifyMesh3d(mesh, passed)
  write(*,'(A,G0,/)') 'VerifyMesh3d: passed = ', passed
  call mesh % GetPoints(po, 'L', x)
  ne = mesh % n_elem
  allocate(s(0:po,0:po,0:po,ne,ns), source = ZERO)
  call ExportVTK_VolumeData( x, s                   &
                           , sname  = ['s']         &
                           , file   = 'naca_0012'   &
                           , part   = mesh % part   &
                           , n_part = mesh % n_part )
  call MPI_Finalize()

end program Mesh3d_Import_Exodus
