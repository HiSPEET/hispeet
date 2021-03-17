program Mesh3d_Import_Generic
  use Kind_Parameters
  use Constants
  use XMPI
  use Generic_Mesh_3d
  use Mesh_3d__Partition
  use Element_Transfer_Buffer_3d
  use Verify_Mesh_3d
  use Assembly_3d
  use Export_VTK_3d__Volume_Data
  implicit none

  character(len=*), parameter :: input_file = 'mesh3d_import_generic.prm'
  integer   :: conf = 1   ! configuration (1 cylinder, 2 annular gap)
  real(RNP) :: r0   = 0.5 ! inner radius  (annular gap only)
  real(RNP) :: r1   = 1   ! outer radius
  real(RNP) :: h    = 2   ! height = axial extension
  integer   :: nr   = 2   ! num elements in radial    direction
  integer   :: np   = 4   ! num elements in azimuthal direction ≥ 3 (gap only)
  integer   :: nz   = 3   ! num elements in axial     direction ≠ 2 if periodic
  integer   :: po   = 3   ! polynomial order of mesh elements
  logical   :: periodic = .true. ! switch for axial periodicity
  namelist/input/ conf, r0, r1, h, nr, np, nz, po, periodic

  logical   :: test_avg   = .false. ! perform averaging test
  logical   :: export_vtk = .false. ! generate VTK file
  namelist/control/ test_avg, export_vtk

  type(MPI_Comm) :: comm = MPI_COMM_WORLD
  integer :: rank

  type(GenericMesh3d)    :: generic_mesh
  type(Mesh3d_Partition) :: mesh
  real(RNP), allocatable :: x(:,:,:,:,:)

  real(RNP), allocatable, target :: var(:,:,:,:,:)
  real(RNP), allocatable         :: v(:,:,:,:)     ! test variable
  real(RNP), pointer             :: r(:,:,:,:)     ! reference variable
  real(RNP), pointer             :: e(:,:,:,:)     ! error

  type(ElementTransferBuffer3d), asynchronous, allocatable :: v_buf

  real(RNP) :: kappa(3), y(3), err
  logical   :: passed
  integer   :: io
  integer   :: i, j, k, l
  integer   :: i_err, j_err, k_err, l_err

  call MPI_Init()
  call MPI_Comm_rank(comm, rank)

  if (rank == 0) then

    ! read parameters ..........................................................

    open(newunit = io, file = input_file)
    read(io, nml = input)
    read(io, nml = control)
    close(io)

    ! create and import generic mesh ...........................................

    select case(conf)
    case(2)
      call generic_mesh % CreateAnnularGap(r0, r1, h, nr, np, nz, po, periodic)
    case default
      call generic_mesh % CreateCylinder(nr, nz, po, periodic)
    end select

    ! verification .............................................................

    call mesh % ImportGenericMesh(generic_mesh, comm = comm)

    call VerifyMesh3d(mesh, passed)
    write(*,'(/,A,G0,/)') 'VerifyMesh3d: passed = ', passed

    ! set up data ..............................................................

    if (test_avg .or. export_vtk) then
      call mesh % GetPoints(po, 'GLL', x)
      allocate(v  (0:po, 0:po, 0:po, mesh%n_elem + mesh%n_ghost))
      allocate(var(0:po, 0:po, 0:po, mesh%n_elem, 2), source = ZERO)
      r(0:,0:,0:,1:) => var(:,:,:,:,1)
      e(0:,0:,0:,1:) => var(:,:,:,:,2)
    end if

    ! averaging test ...........................................................

    if (test_avg) then

      kappa = 2 * PI / h

      do l = 1, mesh%n_elem
        do k = 0, po
        do j = 0, po
        do i = 0, po
          y(1) = x(i,j,k,l,1)
          y(2) = x(i,j,k,l,2)
          y(3) = x(i,j,k,l,3)
          r(i,j,k,l) = cos(sum(kappa * y))
          v(i,j,k,l) = r(i,j,k,l)
        end do
        end do
        end do
      end do

      v_buf = ElementTransferBuffer3d(mesh, v)
      call Assembly3d(mesh, v, v_buf, avg=.true.)

      err = 0
      do l = 1, mesh%n_elem
        do k = 0, po
        do j = 0, po
        do i = 0, po
          e(i,j,k,l) = v(i,j,k,l) - r(i,j,k,l)
          if (abs(e(i,j,k,l)) > err) then
            err   = abs(e(i,j,k,l))
            i_err = i
            j_err = j
            k_err = k
            l_err = l
          end if
        end do
        end do
        end do
      end do
      write(*,'(A,ES10.3,A,4I5)') 'Average over element boundaries: err = ', &
                                  err, ' at ', i_err, j_err, k_err, l_err
    end if

    ! export mesh and data .....................................................

    if (export_vtk) then
      call ExportVTK_VolumeData( x, var                  &
                               , sname  = ['r','e']      &
                               , file   = 'element_mesh' &
                               , part   = mesh % part    &
                               , n_part = mesh % n_part  )

      call mesh % GetCuboids(x)
      call ExportVTK_VolumeData( x                       &
                               , file   = 'cuboid_mesh'  &
                               , part   = mesh % part    &
                               , n_part = mesh % n_part  )

    end if

  end if

  call MPI_Finalize()

end program Mesh3d_Import_Generic
