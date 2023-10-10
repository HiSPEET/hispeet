!> summary:  Test program for importing Gmsh mesh files
!> author:   Matthias Frey, Joerg Stiller
!> date:     2023/05/10
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================
!-------------------------------------------------------------------------------

program Mesh3d_Import_GMSH
  use Kind_Parameters
  use Constants
  use XMPI
  use Import__GMSH__3D
  use Generic_Mesh__3D
  use Mesh__3D
  use Verify_Mesh__3D
  use Export_VTK_Volume_Data__3D
  implicit none

  character(len=80) :: gmsh_file = 'gmsh/example.msh'

  type(GenericMesh_3D) :: generic_mesh
  type(Mesh_3D) :: mesh
  real(RNP), allocatable :: x(:,:,:,:,:)
  real(RNP), allocatable :: s(:,:,:,:,:)
  integer :: po = 1, ne, ns = 1
  logical :: passed

  call MPI_Init()
  call Import_GMSH_3D(gmsh_file, generic_mesh)
! call generic_mesh % SwitchToLexicalNumbering() ! do that when importing
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

end program Mesh3d_Import_GMSH
