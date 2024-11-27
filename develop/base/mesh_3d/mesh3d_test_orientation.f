program Mesh3d_Test_Orientation
  use Kind_Parameters
  use Constants
  use XMPI
  use Standard_Element_Operators__1D
  use Generic_Mesh__3D
  use Mesh__3D
  use Mesh_Element__3D
  use Verify_Mesh__3D
  use Export_VTK_Volume_Data__3D
  implicit none

  character(len=*), parameter :: input_file = 'mesh3d_test_orientation.prm'

  real(RNP) :: xo(3)       = 0       ! corner closest to -infinity
  real(RNP) :: lx(3)       = 1       ! domain extensions
  integer   :: po          = 3       ! polynomial order of mesh elements
  integer   :: rotation(3) = 0       ! x/y/z-rotation of center element
  logical   :: periodic(3) = .false. ! F/T for non/periodic directions
  logical   :: export_vtk  = .true.  ! generate VTK file
  namelist/input/ xo, lx, po, rotation, periodic, export_vtk

  type(MPI_Comm) :: comm = MPI_COMM_WORLD
  integer :: rank

  type(GenericMesh_3D) :: generic_mesh
  type(Mesh_3D)        :: mesh
  type(StandardElementOperators_1D) :: std_op

  real(RNP), allocatable, target :: var(:,:,:,:,:)
  real(RNP), allocatable  :: x(:,:,:,:,:)  ! element point coordinates
! real(RNP), allocatable  :: v(:,:,:,:)    ! test variable
  real(RNP), pointer      :: xi(:,:,:,:,:) ! element point std coordinates
  character(len=80) :: stdin
  logical   :: passed
  integer   :: io
  integer   :: l, m, n

  call Init_MPI_Binding()
  call MPI_Comm_rank(MPI_COMM_WORLD, rank)

  if (rank == 0) then

    ! set communicator to run tests on rank 0 only .............................

    comm = MPI_COMM_SELF

    ! read parameters ..........................................................

    open(newunit = io, file = input_file)
    read(io, nml = input)
    close(io)

    std_op = StandardElementOperators_1D(po, nodes='L')

    ! create and import generic mesh ...........................................

    call generic_mesh % CreateOneRotated(xo, lx, po, rotation, periodic)
    call mesh % ImportGenericMesh(generic_mesh, comm = comm)

    ! verification .............................................................

    call VerifyMesh_3D(mesh, passed)
    write(*,'(/,A,G0,/)') 'VerifyMesh_3D: passed = ', passed

    write(*,'(A,/)') 'Review of element information'
    ELEMENT_INFO: do
      write(*,'(2X,2A)',advance='NO') 'triple index (press ENTER to exit): ', &
                                      '3 ≥ l,m,n = '
      read(*,'(A)') stdin
      if (len_trim(stdin) == 0) exit
      read(stdin,*) l, m, n
      l = max(1, min(27, l + 3*(m-1 + 3*(n-1))))
      write(*,'(/,2X,3G0)') 'element(',l,'):'
      associate(element => mesh % element(l))
        write(*,*)
        do m = 1, 6
          write(*,'(4X,A,I2,A,4I4)', advance = 'NO') 'face', m  &
              , ': normal, rotation, n_neighbor, neighbor orientation ='  &
              , element % face(m) % normal         &
              , element % face(m) % rotation       &
              , element % face(m) % n_neighbor
          call PrintNeighborOrientation( element % neighbor             &
                                       , element % face(m) % i_neighbor &
                                       , element % face(m) % n_neighbor )
        end do
        write(*,*)
        do m = 1, 12
          write(*,'(4X,A,I3,A,4I4)', advance = 'NO') 'edge', m  &
            , ': orientation, n_neighbor, neighbor orientations =' &
            , element % edge(m) % orientation    &
            , element % edge(m) % n_neighbor
          call PrintNeighborOrientation( element % neighbor             &
                                       , element % edge(m) % i_neighbor &
                                       , element % edge(m) % n_neighbor )
        end do
        write(*,*)
        do m = 1, 8
          write(*,'(4X,A,I2,A,4I4)', advance = 'NO') 'vertex', m  &
            , ': n_neighbor, neighbor orientations =' &
            , element % vertex(m) % n_neighbor
          call PrintNeighborOrientation( element % neighbor               &
                                       , element % vertex(m) % i_neighbor &
                                       , element % vertex(m) % n_neighbor )
        end do
        write(*,*)
      end associate
    end do ELEMENT_INFO

    ! set up data ..............................................................

    allocate(var(0:po, 0:po, 0:po, mesh%n_elem, 3), source = ZERO)
!   allocate(v  (0:po, 0:po, 0:po, mesh%n_elem + mesh%n_ghost))

    ! element point coordinates
    call mesh % GetPoints(po, 'L', x)

    ! element point standard coordinates
    xi(0:,0:,0:,1:,1:) => var(:,:,:,:,1:3)
    do l = 1, mesh % n_elem
      do n = 0, po
      do m = 0, po
        xi(:,m,n,l,1) = std_op % x
        xi(m,:,n,l,2) = std_op % x
        xi(m,n,:,l,3) = std_op % x
      end do
      end do
    end do

    ! export mesh and data .....................................................

    if (export_vtk) then
      call ExportVTK_VolumeData( x, var                        &
                               , sname   = ['xi1','xi2','xi3'] &
                               , file    = 'element_mesh'      &
                               , part    = mesh % part         &
                               , n_parts = mesh % n_parts      )
    end if

  end if

  call MPI_Finalize()

contains

  subroutine PrintNeighborOrientation(neighbor, i1, nn)
    class(MeshElementNeighbor_3D), intent(in) :: neighbor(:)
    integer(IXS), intent(in) :: i1 !< start index
    integer(IXS), intent(in) :: nn !< number of neighbors

    if (nn > 0) then
      write(*,'(99I4)') neighbor(i1:i1+nn-1) % orientation
    else
      write(*,*)
    end if

  end subroutine PrintNeighborOrientation

end program Mesh3d_Test_Orientation
